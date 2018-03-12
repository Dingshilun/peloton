//===----------------------------------------------------------------------===//
//
//                         Peloton
//
// skiplist.h
//
// Identification: src/include/index/skiplist.h
//
// Copyright (c) 2015-17, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#ifndef _SKIPLIST_H
#define _SKIPLIST_H

#pragma once

#include <atomic>
#include <functional>
#include <thread>
#include <deque>

#include "index/index.h"

namespace peloton {
namespace index {

#define WORD(x) reinterpret_cast<std::uintptr_t>((x))
#define GET_DELETE(addr) ((WORD((addr))) & 1ll)
#define GET_FLAG(addr) ((WORD((addr))) & 2ll)
#define GET_SUCC(node) ((node)->next_.load())
#define CHECK_DELETE(node) (WORD(GET_SUCC(node)) & 1ll)
#define CHECK_FLAG(node) (WORD(GET_SUCC(node)) & 2ll)
#define SET_DELETE(addr, bit) (((WORD((addr))) & ~1ll) | (bit))
#define SET_FLAG(addr, bit) (((WORD((addr))) & ~2ll) | ((bit) << 1))
/*******************************/
#define PHYSICALLY_DELETE(addr) ((WORD(addr)) | 4ll)
#define GET_PHYSICALLY_DELETE(addr) ((WORD(addr)) & 4ll)
#define DISABLE_CLEAR_EPOCH false
#define DISABLE_VERIFY true
/*******************************/
#define GET_NEXT(node) \
  reinterpret_cast<SkipListBaseNode *>(WORD((node)->next_.load()) & ~7ll)

#define SKIP_LIST_INITIAL_MAX_LEVEL_ 32
#define MAX_THREAD_COUNT ((int)0x7FFFFFFF)

/*
 * SKIPLIST_TEMPLATE_ARGUMENTS - Save some key strokes
 */
#define SKIPLIST_TEMPLATE_ARGUMENTS                                       \
  template <typename KeyType, typename ValueType, typename KeyComparator, \
            typename KeyEqualityChecker, typename ValueEqualityChecker>
template <typename KeyType, typename ValueType, typename KeyComparator,
          typename KeyEqualityChecker, typename ValueEqualityChecker>
class SkipList {
 public:
  class NodeManager;
  class EpochManager;
  class OperationContext;
  class SkipListBaseNode;
  class SkipListInnerNode;

  using NodePair = std::pair<SkipListBaseNode *, SkipListBaseNode *>;
  using NodeList = std::vector<SkipListBaseNode *>;
  using InnerNodePair = std::pair<SkipListInnerNode *, SkipListInnerNode *>;
  using InnerNodeList = std::vector<SkipListInnerNode *>;

  using KeyValuePair = std::pair<KeyType *, ValueType *>;

 private:
  ///////////////////////////////////////////////////////////////////
  // Core components
  ///////////////////////////////////////////////////////////////////
  std::atomic<SkipListBaseNode *> skip_list_head_;
  u_int32_t max_level_;
  NodeManager node_manager_;
  EpochManager epoch_manager_;
  bool duplicate_support_;
  int GC_Interval_;

  /*
   * Get() - Search a key in the skip-list and fill in the value_list
   *
   * The return value is a indicator of the success get
   */
  bool Get(const KeyType &key, std::vector<ValueType> &value_list,
           OperationContext &ctx) {
    // LOG_INFO("Get()");
    auto pair = Search(key, ctx);
    auto node = pair.second;
    while (node != nullptr && KeyCmpEqual(node->key_, key)) {
      if (CHECK_DELETE(node)) {
        node = GET_NEXT(node);
        continue;
      }
      value_list.push_back(static_cast<SkipListInnerNode *>(node)->GetValue());
      node = GET_NEXT(node);
    }
    return true;
  }

  /*
   * GetLimit() - Search a key in the skip-list and fill in the value_list with
   * offset and limit
   *
   * The return value is a indicator of the success get
   */
  bool GetLimit(const KeyType &key, std::vector<ValueType> &value_list,
                uint64_t limit, uint64_t offset, OperationContext &ctx) {
    LOG_INFO("Get()");
    auto pair = Search(key, ctx);
    auto node = pair.second;
    uint64_t count = 0;
    while (node != nullptr && KeyCmpEqual(node->key_, key) &&
           count < offset + limit) {
      if (CHECK_DELETE(node)) {
        node = GET_NEXT(node);
        continue;
      }
      if (count >= offset) {
        value_list.push_back(
            static_cast<SkipListInnerNode *>(node)->GetValue());
      }
      node = GET_NEXT(node);
      count++;
    }
    return true;
  }

  /*
   * Search() - Search the first interval that node1.key < key and key <=
   * node2.key
   *
   * The return value is a pair of node1, node2
   * if duplicate is available, the node2 is the first node among all
   * duplicators
   *
   * NOTE: the second pointer might be nullptr!!!!!!!!
   */
  NodePair Search(const KeyType &key, OperationContext &ctx) {
    SkipListBaseNode *headNode = this->skip_list_head_.load();
    while (1) {
      auto sr = SearchFrom(key, headNode, ctx);
      PL_ASSERT(sr.first != nullptr);
      headNode = sr.first;
      if (headNode->down_.load() == nullptr) {
        return sr;
      } else {
        headNode = headNode->down_.load();
      }
    }
  }
  /*
   * SearchPrevNode() - search for the target_node in the skiplist
   *
   * will return its prev_node in the skiplist
   * might return a deleted previous node, so we need to double check outside
   * and the previous node might be a head(base) node, also need to check this
   */
  SkipListInnerNode *SearchPrevNode(SkipListInnerNode *target_node,
                                    OperationContext &ctx) {
    SkipListInnerNode *cursor = nullptr;
    KeyType target_key = target_node->key_;
    auto pair = Search(target_key, ctx);
    while (pair.second != target_node) {
      if (CHECK_DELETE(target_node)) {
        return reinterpret_cast<SkipListInnerNode *>(
            target_node->back_link_.load());
      }
      cursor = reinterpret_cast<SkipListInnerNode *>(pair.second);

      while (cursor) {
        if (!KeyCmpEqual(cursor->key_, target_node->key_)) {
          break;
        }
        if (CHECK_DELETE(cursor)) {
          cursor = reinterpret_cast<SkipListInnerNode *>(GET_NEXT(cursor));
        }
        if (CHECK_FLAG(cursor)) {
          HelpFlagged(cursor, GET_NEXT(cursor), ctx);
        }
        if (GET_NEXT(cursor) == target_node) {
          return cursor;
        } else {
          cursor = reinterpret_cast<SkipListInnerNode *>(GET_NEXT(cursor));
        }
      }
      pair = Search(target_key, ctx);
    }
    return reinterpret_cast<SkipListInnerNode *>(pair.first);
  }

  /*
   * GetPrevNode - Return previous node in the linked list
   *
   * will return null if there is not such pointer
   */
  SkipListBaseNode *GetPrevNode(SkipListBaseNode *target_node,
                                SkipListBaseNode *start) {
    auto cursor = start;
    while (cursor) {
      if (!cursor->isHead_ && KeyCmpGreater(cursor->key_, target_node->key_))
        break;
      if (CHECK_DELETE(cursor)) {
        cursor = cursor->back_link_.load();
      }
      if (GET_NEXT(cursor) == target_node) {
        return cursor;
      }
      cursor = GET_NEXT(cursor);
    }
    return nullptr;
  }
  /*
   * GetPrevNodeWithRootNode - Return previous node with the root
   *
   * will return null if there is no such node
   */
  NodePair GetPrevNodeWithRootNode(SkipListBaseNode *root,
                                   SkipListBaseNode *start) {
    auto cursor = start;
    while (cursor) {
      if (!cursor->isHead_ && KeyCmpGreater(cursor->key_, root->key_)) break;
      if (CHECK_DELETE(cursor)) {
        cursor = cursor->back_link_.load();
      }
      if (GET_NEXT(cursor)) {
        if (reinterpret_cast<SkipListInnerNode *>(GET_NEXT(cursor))
                ->GetRoot()
                .load() == root) {
          return std::make_pair(cursor, GET_NEXT(cursor));
        }
      }
      cursor = GET_NEXT(cursor);
    }
    return std::make_pair(nullptr, nullptr);
  }
  /*
   * SearchFrom() - Search the first interval that node1.key < key and
   * key <= node2.key
   * from given skip list node
   *
   * The return value is a pair of node1, node2
   * For duplicate enabled skip list, the type of return value is the same as
   * Search()
   * There is no guarantee that the nodes would be succeed after being returned
   *
   * Call this function again in insert and delete if the node pair is not
   * consistent (node1.next != node2)
   */
  NodePair SearchFrom(const KeyType &key, const SkipListBaseNode *Node,
                      UNUSED_ATTRIBUTE OperationContext &ctx) {
    // TODO: physically deletion when search in the list
    // LOG_INFO("SearchFrom key:%s", key.GetInfo().c_str());
    if (Node == nullptr) {
      return std::make_pair(nullptr, nullptr);
    }
    SkipListBaseNode *curr_node = const_cast<SkipListBaseNode *>(Node);
    while (curr_node) {
      SkipListBaseNode *tmp_pointer = curr_node->next_.load();
      if (GET_FLAG(tmp_pointer)) {
        // HelpFlagged(curr_node, GET_NEXT(curr_node), ctx);
      } else if ((GET_DELETE(tmp_pointer))) {
        curr_node = curr_node->back_link_.load();
        // LOG_INFO("GO BACK: %p", curr_node);
      } else if (tmp_pointer == nullptr) {
        return std::make_pair(curr_node, nullptr);
      } else {
        if (KeyCmpGreaterEqual(tmp_pointer->key_, key)) {
          return std::make_pair(curr_node, tmp_pointer);
        } else {
          curr_node = tmp_pointer;
        }
      }
    }
    return std::make_pair(nullptr, nullptr);
  }

  /*
   * SearchWithPath() - Search the skiplist for the key, would store the path of
   * every level
   *
   * @param:
   *  call_stack: used for storing the path
   *  key: the search key
   *  curr_node: the same as SearchFrom, but please send in a SkipListHead
   *  expected_stored_level: from which level the function starts to record the
   * path, default to start recording from
   *  curr_node's level
   *
   * The lowest level starts from 0, the search would start from the curr_node
   * The function defaults to store all nodes in the path from the head to
   * target node
   *
   * returns nothing but will store the path at call_stack
   */
  void SearchWithPath(std::vector<NodePair> &call_stack, const KeyType &key,
                      SkipListBaseNode *curr_node, OperationContext &ctx,
                      u_int32_t expected_stored_level = 0) {
    expected_stored_level = (expected_stored_level == 0)
                                ? curr_node->level_
                                : expected_stored_level;
    u_int32_t level_now = curr_node->level_;
    call_stack.resize(expected_stored_level + 1);
    // LOG_INFO("SearchWithPath %d, levelNow: %u", expected_stored_level,
    //         level_now);
    while (level_now >= 0) {
      if (level_now <= expected_stored_level) {
        call_stack[level_now] = SearchFrom(key, curr_node, ctx);
        curr_node = call_stack[level_now].first->down_.load();
      } else {
        curr_node = SearchFrom(key, curr_node, ctx).first->down_.load();
      }
      // Stop at root level
      if (level_now == 0) break;
      level_now--;
    }
  }

  /*
 * AddLevel() - add corresponding level to the SkipList
 *
 * return true if successfully added or the level is already added
 * return false if the level cannot be reached from the highest level now
 */
  bool AddLevel(u_int32_t level) {
    // LOG_INFO("AddLevel %u", level);
    SkipListBaseNode *head = this->skip_list_head_.load();
    if (head->level_ + 1 < level) {
      return false;
    } else {
      if (head->level_ + 1 == level) {
        SkipListBaseNode *new_head = node_manager_.GetSkipListHead(level);
        new_head->down_ = head;
        if (this->skip_list_head_.compare_exchange_strong(head, new_head)) {
          return true;
        } else {
          node_manager_.ReturnSkipListNode(new_head);
          head = this->skip_list_head_.load();
          return head->level_ == level;
        }
      } else {
        return true;
      }
    }
  }

  /*
   * InsertTowerIntoInterval() - this method would try to insert the tower into
   * the interval and retry due to contention
   *
   * It has a call_stack array to accelorate the process
   * NOTE: this method would retry until the world ends (or the root is
   * deleted)!!!
   *
   * Would only return true, or it would retry until succeed or the root is
   * deleted
   * Use this function only the tower can be exactly inserted
   */
  bool InsertTowerIntoInterval(
      const KeyType &key, std::vector<SkipListInnerNode *> &tower,
      std::vector<std::pair<SkipListBaseNode *, SkipListBaseNode *>> &
          call_stack,
      OperationContext &ctx, u_int32_t start_level = 0,
      bool check_multiple_key_value = false,
      std::function<bool(const void *)> predicate = nullptr,
      bool *predicate_satisfied = nullptr) {
    // LOG_INFO("InsertTower");
    u_int32_t expected_level = tower.size();
    for (u_int32_t i = start_level; i < expected_level; i++) {
      bool insert_flag = false;
      do {
        //        if (i != 0 && CHECK_DELETE(tower[i]->GetRoot().load())) {
        //          // the root has been deleted
        //          // there is no need to continue
        //          for (auto j = i; j < expected_level; j++) {
        //            epoch_manager_.AddGarbageNode(tower[j]);
        //          }
        //          return true;
        //        }

        // if the level is 0, multiple key-value pair is required
        // if the check bool is true, then it should be non-unique index
        if (i == 0) {
          if (!check_multiple_key_value) {
            // unique index, just need to check the next one
            if (!GET_DELETE(call_stack[i].second->next_.load()) &&
                ValueCmpEqual(tower[i]->GetValue(),
                              static_cast<SkipListInnerNode *>(
                                  call_stack[i].second)->GetValue())) {
              if (predicate_satisfied != nullptr) *predicate_satisfied = false;
              for (auto tmp : tower) {
                node_manager_.ReturnSkipListNode(tmp);
              }
              return false;
            }
          } else {
            // There can be node with the same key and different values, we need
            // to ensure no nodes with both key and value are the same
            // TODO: add optimization of last pointer check
            auto cursor =
                static_cast<SkipListInnerNode *>(call_stack[i].second);
            while (cursor) {
              if (!KeyCmpEqual(key, cursor->key_)) break;
              if (CHECK_DELETE(cursor)) {
                cursor = static_cast<SkipListInnerNode *>(GET_NEXT(cursor));
                continue;
              }
              if (ValueCmpEqual(tower[i]->GetValue(), cursor->GetValue())) {
                if (predicate_satisfied != nullptr)
                  *predicate_satisfied = false;
                for (auto tmp : tower) {
                  node_manager_.ReturnSkipListNode(tmp);
                }
                return false;
              }
              if (predicate != nullptr && predicate(cursor->GetValue())) {
                *predicate_satisfied = false;
                for (auto tmp : tower) {
                  node_manager_.ReturnSkipListNode(tmp);
                }
                return false;
              }
              cursor = static_cast<SkipListInnerNode *>(GET_NEXT(cursor));
            }
          }
        }

        // if multiple test passed, try to insert
        tower[i]->next_ = call_stack[i].second;
        insert_flag = call_stack[i].first->next_.compare_exchange_strong(
            call_stack[i].second, tower[i]);
        if (insert_flag)
          break;
        else
          call_stack[i] = SearchFrom(key, call_stack[i].first, ctx);
      } while (!insert_flag);
    }
    // verify the whole tower to not be deleted
    if (predicate_satisfied != nullptr) *predicate_satisfied = true;
    if (CHECK_DELETE(tower[0])) {
      // go through the whole tower to verify
      for (int i = tower.size() - 1; i >= 0; i++) {
        if (CHECK_DELETE(tower[i])) break;
        while (!CHECK_DELETE(tower[i])) {
          // try to delete this node as the deleter missed it
          DeleteSingleNode(tower[0], call_stack[i].first, ctx);
        }
      }
    }
    return true;
  }
  /*
   * InsertNode() - Insert key value tuple to the skip-list
   *
   * The return value is a indicator of success or not
   */
  bool InsertNode(const KeyType &key, const ValueType &value,
                  OperationContext &ctx,
                  std::function<bool(const void *)> predicate = nullptr,
                  bool *predicate_satisfied = nullptr) {
    // LOG_INFO("Insert node");

    u_int32_t expected_level = 0;

    while (expected_level < max_level_) {
      // TODO: try not to call too many rand()
      if (rand() & 1) {
        expected_level++;
      } else {
        break;
      }
    }

    SkipListBaseNode *curr_node = this->skip_list_head_.load();

    while (curr_node->level_ < expected_level) {
      AddLevel(curr_node->level_ + 1);
      curr_node = this->skip_list_head_.load();
    }

    // used to store the path
    std::vector<NodePair> call_stack;
    InnerNodeList tower(expected_level + 1);
    // build the tower of expected level
    SkipListInnerNode *new_node =
        node_manager_.GetSkipListInnerNode(key, value, 0);
    tower[0] = new_node;
    for (u_int32_t i = 1; i < expected_level + 1; i++) {
      tower[i] =
          node_manager_.GetSkipListInnerNode(key, tower[0], tower[i - 1], i);
    }
    SearchWithPath(call_stack, key, curr_node, ctx, expected_level);
    PL_ASSERT(curr_node != nullptr);
    // if duplicate support, then just try insert
    // else need to verify the next node
    if (this->duplicate_support_) {
      // insert the node from the lowest level
      // redo the search from stack if the insert fails
      return InsertTowerIntoInterval(key, tower, call_stack, ctx, 0, true,
                                     predicate, predicate_satisfied);
    } else {
      // unique key
      // need to compare with the second return value's key
      bool insert_flag;
      do {
        // try to insert the key in the lowest level
        // if failed then abort the insert
        if (call_stack[0].second == nullptr ||
            CHECK_DELETE(call_stack[0].second) ||
            !KeyCmpEqual(call_stack[0].second->key_, key)) {
          tower[0]->next_ = call_stack[0].second;
          insert_flag = call_stack[0].first->next_.compare_exchange_strong(
              call_stack[0].second, tower[0]);
        } else {
          // found duplicate key not deleted
          // abort the insertion
          if (predicate_satisfied != nullptr) *predicate_satisfied = false;
          for (auto tmp : tower) {
            node_manager_.ReturnSkipListNode(tmp);
          }
          return false;
        }
        if (insert_flag) break;
        call_stack[0] = SearchFrom(key, call_stack[0].first, ctx);
      } while (!insert_flag);
      // insertion at the lowest level has succeeded
      // those towers should all be inserted into the skiplist successfully
      return InsertTowerIntoInterval(key, tower, call_stack, ctx, 1);
    }
  }

  /*
   * SearchKeyValueInList() - Search specific key value pair in the current
   * level
   *
   * return the prev of the node to delete
   */
  NodePair SearchKeyValueInList(const KeyType &key, const ValueType &value,
                                SkipListBaseNode *prev, SkipListBaseNode *del) {
    // LOG_INFO("Search Key Value in list");
    PL_ASSERT(prev != nullptr);
    while (del && KeyCmpEqual(del->key_, key)) {
      auto inner_node = static_cast<SkipListInnerNode *>(del);
      if (ValueCmpEqual(inner_node->GetRootValue(), value)) {
        return std::make_pair(prev, del);
      } else {
        prev = del;
        del = GET_NEXT(del);
      }
    }
    return std::make_pair(nullptr, nullptr);
  }

  /*
   * DeleteSingleNode - delete one node with specific root node
   */
  void DeleteSingleNode(SkipListBaseNode *root, SkipListBaseNode *start,
                        OperationContext &ctx) {
    NodePair tmp_pair = GetPrevNodeWithRootNode(root, start);
    auto target = tmp_pair.second;
    auto previous = tmp_pair.first;
    bool indicator = false;
    // LOG_INFO("getting int delete single node: %p, %p", root, start);
    while (target) {
      // found something to delete
      // flag its previous node to notify other deleters
      indicator = TryFlag(previous, target, ctx);
      if (indicator) {
        // flag succeed, we can try to delete the node
        HelpFlagged(previous, target, ctx);
        // LOG_INFO("deleted previous:%p target:%p", previous, target);
        break;
      }
      // if we didn't find the previous node, then the node has been deleted
      // from skiplist
      previous = GetPrevNode(target, start);
    }
  }
  /*
   * DeleteNode() - Try delete the key value pair from the skip list
   *
   * return true if successful deleted
   *
   * Modified By Shilun: delete the node from top down
   */
  bool DeleteNode(const KeyType &key, UNUSED_ATTRIBUTE const ValueType &value,
                  NodePair &pair, OperationContext &ctx) {
    // LOG_INFO("DeleteNode()");
    SkipListBaseNode *prev_node = pair.first;
    SkipListBaseNode *del_node = pair.second;
    bool result = false;
    if (prev_node && del_node) {
      // mark the delete node first
      // in case that it belongs to a tower that is being inserted
      SkipListBaseNode *tmp = nullptr;
      // only one thread can get into the deleting stage
      // others wouldn't try to delete the tower as the back_link is not nullptr
      if (del_node->back_link_.compare_exchange_strong(tmp, prev_node)) {
        result = TryDelete(del_node, ctx);
        std::vector<NodePair> call_stack;
        auto head = this->skip_list_head_.load();
        SearchWithPath(call_stack, key, head, ctx);
        // delete the tower from top if possible
        for (int i = call_stack.size() - 1; i >= 1; i--) {
          DeleteSingleNode(del_node, call_stack[i].first, ctx);
        }
        // delete the lowest node by ourselves
        auto previous = GetPrevNode(del_node, call_stack[0].first);
        while (previous) {
          bool indicator = TryFlag(previous, del_node, ctx);
          if (indicator) {
            HelpDeleted(previous, del_node, ctx);
            result = true;
            break;
          }
          previous = GetPrevNode(del_node, call_stack[0].first);
        }
      }
    }
    return result;
  }
  /*
   * Delete() - Delete certain key from the skip-list and fill in the
   * deleted nodes to del_nodes
   *
   * Return true if delete succeeded, false otherwise
   */
  bool Delete(const KeyType &key, const ValueType &value,
              OperationContext &ctx) {
    // LOG_INFO("Delete()");
    NodePair pair = Search(key, ctx);
    SkipListBaseNode *prev_node = pair.first;
    SkipListBaseNode *del_node = pair.second;

    // TODO: Improve efficiency
    while (del_node) {
      // Keep searching in the root level
      pair = SearchFrom(key, prev_node, ctx);
      prev_node = pair.first;
      del_node = pair.second;
      if (del_node && !CHECK_FLAG(del_node) && !CHECK_DELETE(del_node)) {
        if (!KeyCmpEqual(key, del_node->key_)) {
          // No such pair
          return false;
        }
        auto root_node = static_cast<SkipListInnerNode *>(del_node);
        if (ValueCmpEqual(root_node->GetValue(), value)) {
          // KV pair found, delete it
          return DeleteNode(key, value, pair, ctx);
        }
        // Continue checking if there are duplicate keys
        if (this->duplicate_support_) {
          prev_node = GET_NEXT(prev_node);
        } else {
          return false;
        }
      }
    }
    return false;
  }

  /*
   * HelpDeleted() Attempts to physically delete the del_node and unflag
   * prev_node
   */
  void HelpDeleted(SkipListBaseNode *prev_node, SkipListBaseNode *del_node,
                   UNUSED_ATTRIBUTE OperationContext &ctx) {
    auto set_ptr = GET_NEXT(del_node);
    auto cmp_ptr = reinterpret_cast<SkipListBaseNode *>(SET_FLAG(del_node, 1));

    while (!prev_node->next_.compare_exchange_strong(cmp_ptr, set_ptr))
      ;
    // TODO: Notify EpochManager for GC
    this->epoch_manager_.AddGarbageNode(del_node);
    // else LOG_INFO("CAS failed in HelpDeleted!!!!!!!");
  }

  /*
   * HelpFlagged() - Attempts to mark and physically delete del_node
   */
  void HelpFlagged(SkipListBaseNode *prev_node, SkipListBaseNode *del_node,
                   OperationContext &ctx) {
    SkipListBaseNode *prev_assert = nullptr;
    // this would always succeed now as only one thread can delete a specific
    // tower
    bool flag =
        del_node->back_link_.compare_exchange_strong(prev_assert, prev_node);
    if (flag) {
      if (!CHECK_DELETE(del_node)) {
        TryDelete(del_node, ctx);
      }
      HelpDeleted(prev_node, del_node, ctx);
    }
  }

  /*
   * TryDelete() Attempts to mark the node del node.
   */
  bool TryDelete(SkipListBaseNode *del_node,
                 UNUSED_ATTRIBUTE OperationContext &ctx) {
    bool ret = false;
    while (!CHECK_DELETE(del_node)) {
      auto cmp_ptr = GET_NEXT(del_node);
      auto set_ptr =
          reinterpret_cast<SkipListBaseNode *> SET_DELETE(cmp_ptr, 1);
      ret = del_node->next_.compare_exchange_strong(cmp_ptr, set_ptr);
      //      cannot do this now, as nodes can only be deleted by its own
      //      deleter
      //      if (CHECK_FLAG(del_node)) {
      //        HelpFlagged(del_node, GET_NEXT(del_node), ctx);
      //      }
    }
    return ret;
  }

  /*
   * TryFlag() - Attempts to flag the prev_node, which is the last node known to
   * be the predecessor of target_node
   *
   * The return value is a tuple of deleted node and the success indicator
   */
  bool TryFlag(SkipListBaseNode *prev_node, SkipListBaseNode *target_node,
               UNUSED_ATTRIBUTE OperationContext &ctx) {
    auto flag_ptr =
        reinterpret_cast<SkipListBaseNode *> SET_FLAG(target_node, 1);
    SkipListBaseNode *cmp_ptr = nullptr;
    while (true) {
      cmp_ptr = reinterpret_cast<SkipListBaseNode *> SET_FLAG(target_node, 0);
      if (prev_node->next_.load() == flag_ptr) {
        return true;
      }
      bool result = prev_node->next_.compare_exchange_strong(cmp_ptr, flag_ptr);
      if (result) {
        return true;
      }
      return false;
    }
  }

 public:
  SkipList(KeyComparator key_cmp_obj = KeyComparator{},
           KeyEqualityChecker key_eq_obj = KeyEqualityChecker{},
           ValueEqualityChecker val_eq_obj = ValueEqualityChecker{})
      : epoch_manager_(50, &node_manager_),
        key_cmp_obj_{key_cmp_obj},
        key_eq_obj_{key_eq_obj},
        value_eq_obj_{val_eq_obj} {
    this->duplicate_support_ = true;
    this->GC_Interval_ = 50;
    this->max_level_ = SKIP_LIST_INITIAL_MAX_LEVEL_;
    this->skip_list_head_ = node_manager_.GetSkipListHead(0);
    this->epoch_manager_.StartEpochThread();
  }
  SkipList(bool duplicate, int GC_Interval,
           KeyComparator key_cmp_obj = KeyComparator{},
           KeyEqualityChecker key_eq_obj = KeyEqualityChecker{},
           ValueEqualityChecker val_eq_obj = ValueEqualityChecker{})
      :  // Construct Epoch Manager
        epoch_manager_(GC_Interval, &node_manager_),
        duplicate_support_(duplicate),
        GC_Interval_(GC_Interval),

        // Key comparator, equality checker
        key_cmp_obj_{key_cmp_obj},
        key_eq_obj_{key_eq_obj},
        value_eq_obj_{val_eq_obj}

  // Value equality checker and hasher
  {
    // LOG_INFO("SkipList constructed!");
    this->max_level_ = SKIP_LIST_INITIAL_MAX_LEVEL_;
    this->skip_list_head_ = node_manager_.GetSkipListHead(0);
    // Start epoch thread
    this->epoch_manager_.StartEpochThread();
  }

  /*
   * when calling this, need to make sure that no other threads exist
   */
  ~SkipList() {
    // TODO: deconstruct all nodes in the skip list
    // LOG_INFO("SkipList deconstructed!");
    // LOG_INFO("footprint: %lu", this->GetMemoryFootprint());
    node_manager_.PrintFootprint();
    auto head = this->skip_list_head_.load();
    std::vector<SkipListBaseNode *> tmp_vec;
    while (head) {
      auto cursor = head;
      while (cursor) {
        tmp_vec.push_back(cursor);
        cursor = GET_NEXT(cursor);
      }
      head = head->down_.load();
      for (auto tmp : tmp_vec) {
        node_manager_.ReturnSkipListNode(tmp);
      }
      tmp_vec.clear();
    }
    // LOG_INFO("footprint: %lu", this->GetMemoryFootprint());
    node_manager_.PrintFootprint();
    return;
  }
  void VerifyList() {
    if (DISABLE_VERIFY) {
      // verify list is more useful in single thread case
      return;
    }
    LOG_INFO("verifying!!!");
    auto head = this->skip_list_head_.load();

    size_t inner_node_count = 0;
    while (head) {
      LOG_INFO("verifying head: %p", head);
      auto cursor = reinterpret_cast<SkipListInnerNode *>(head->next_.load());
      while (cursor) {
        ++inner_node_count;
        LOG_INFO("verifying %p", cursor);
        auto next_cursor = reinterpret_cast<SkipListInnerNode *>(GET_NEXT(cursor));
        auto down_cursor = reinterpret_cast<SkipListInnerNode *>(cursor->down_.load());
        // verity flag bit and marker bit, should not see any delete node and flag node
        PL_ASSERT(next_cursor != cursor);
        PL_ASSERT(!CHECK_FLAG(cursor));
        PL_ASSERT(!CHECK_DELETE(cursor));
        if (next_cursor != nullptr) {
          PL_ASSERT(KeyCmpLessEqual(cursor->key_, next_cursor->key_));
        }
        if (down_cursor != nullptr) {
          PL_ASSERT(KeyCmpEqual(cursor->key_, down_cursor->key_));
        }
        cursor = next_cursor;
      }
      head = head->down_.load();
    }
    size_t count_actual = this->node_manager_.GetInnderNodeCount();
    PL_ASSERT(inner_node_count == count_actual);
  }
  /*
   * struct shouldn't exceed 64 bytes -- cache line
   * possible optimization: add a direct link to the root of the skiplist
   */
  class SkipListBaseNode {
   public:
    std::atomic<SkipListBaseNode *> next_, down_, back_link_;
    KeyType key_;
    bool isHead_;
    u_int32_t level_;

    SkipListBaseNode(bool isHead, u_int32_t level)
        : isHead_(isHead), level_(level) {
      this->next_ = this->down_ = this->back_link_ = nullptr;
    }

    SkipListBaseNode(SkipListBaseNode *next, SkipListBaseNode *down,
                     SkipListBaseNode *back_link, KeyType key, bool isHead)
        : next_(next),
          down_(down),
          back_link_(back_link),
          key_(key),
          isHead_(isHead) {}

    SkipListBaseNode(SkipListBaseNode *next, SkipListBaseNode *down,
                     SkipListBaseNode *back_link, KeyType key, bool isHead,
                     u_int32_t level)
        : next_(next),
          down_(down),
          back_link_(back_link),
          key_(key),
          isHead_(isHead),
          level_(level) {}

    virtual ~SkipListBaseNode(){};
  };

  class SkipListInnerNode : public SkipListBaseNode {
   public:
    SkipListInnerNode(SkipListBaseNode *next, SkipListBaseNode *down,
                      SkipListBaseNode *back_link, KeyType key, bool isHead)
        : SkipListBaseNode(next, down, back_link, key, isHead) {}
    SkipListInnerNode(SkipListBaseNode *next, SkipListBaseNode *down,
                      SkipListBaseNode *back_link, KeyType key, bool isHead,
                      u_int32_t level)
        : SkipListBaseNode(next, down, back_link, key, isHead, level) {}

    // set the union value
    void SetValue(ValueType value) {
      PL_ASSERT(this->down_ == nullptr);
      this->valueOrRoot_.value = value;
    }

    // set the union root
    void SetRoot(SkipListInnerNode *root) {
      PL_ASSERT(this->down_ != nullptr);
      this->valueOrRoot_.root = root;
    }

    // Get the value from the tower
    ValueType &GetRootValue() {
      if (this->down_ == nullptr) {
        return this->valueOrRoot_.value;
      } else {
        SkipListInnerNode *root = this->valueOrRoot_.root.load();
        return root->valueOrRoot_.value;
      }
    }

    ValueType &GetValue() {
      PL_ASSERT(this->down_ == nullptr);
      return this->valueOrRoot_.value;
    }

    std::atomic<SkipListInnerNode *> &GetRoot() {
      PL_ASSERT(this->down_ != nullptr);
      return this->valueOrRoot_.root;
    }

   private:
    // value when down is null
    // otherwise root
    union valueOrRoot {
      ValueType value;
      std::atomic<SkipListInnerNode *> root;
    } valueOrRoot_;
  };

  /*
   * class ForwardIterator - Iterator that supports forward iteration of the
   * skip list
   *
   */
  // TODO: modify the context
  class ForwardIterator {
   private:
    SkipListInnerNode *node_;
    KeyValuePair kv_p;
    SkipList *list_;
    OperationContext ctx_;

   public:
    /*
     * Default constructor - Create a default forward iterator
     */
    ForwardIterator() : node_{nullptr}, kv_p{}, list_{nullptr}, ctx_(nullptr) {}

    /*
     * Constructor - Create a forward iterator based on that skip list
     */
    ForwardIterator(SkipList *list) : list_{list}, ctx_(nullptr) {
      this->IterJoinEpoch();
      auto head = list->skip_list_head_.load();
      while (head->down_) head = head->down_;

      if (GET_NEXT(head) != nullptr) {
        node_ = reinterpret_cast<SkipListInnerNode *>(GET_NEXT(head));
      } else {
        node_ = nullptr;
      }

      kv_p = std::make_pair(&(node_->key_), &(node_->GetValue()));

      // list_->epoch_manager_.LeaveEpoch(epoch_node_p);
    }

    /*
     * Constructor - Create a forward iterator start from start_key in the list
     */

    ForwardIterator(SkipList *list, const KeyType &start_key)
        : list_(list), ctx_(nullptr) {
      this->IterJoinEpoch();

      // Get the pointer of the first node
      auto root_pair = list->Search(start_key, ctx_);
      if (root_pair.second != nullptr) {
        node_ =
            reinterpret_cast<SkipListInnerNode *>(GET_NEXT(root_pair.first));
        kv_p = std::make_pair(&(node_->key_), &(node_->GetValue()));
      } else {
        node_ = nullptr;
      }

      // list_->epoch_manager_.LeaveEpoch(epoch_node_p);
    }
    ~ForwardIterator() {}

    void IterJoinEpoch() {
      auto epoch_node = this->list_->epoch_manager_.JoinEpoch();
      this->ctx_.epoch_node_ = epoch_node;
    }

    void IterLeaveEpoch() {
      this->list_->epoch_manager_.LeaveEpoch(this->ctx_.epoch_node_);
    }
    /*
     * IsEnd() - If it's the end of the list
     */
    bool IsEnd() { return node_ == nullptr; }

    inline const KeyValuePair *operator->() { return &kv_p; }

    inline ForwardIterator &operator++() {
      if (IsEnd()) {
        return *this;
      }

      // auto epoch_node_p = list_->epoch_manager_.JoinEpoch();
      // OperationContext ctx{epoch_node_p};

      // TODO: Add logic for handling delete node
      node_ = reinterpret_cast<SkipListInnerNode *>(GET_NEXT(node_));
      if (node_) {
        if (CHECK_DELETE(node_)) {
          node_ = reinterpret_cast<SkipListInnerNode *>(GET_NEXT(node_));
        }
        kv_p.first = &(node_->key_);
        kv_p.second = &(node_->GetValue());
      }

      // list_->epoch_manager_.LeaveEpoch(epoch_node_p);

      return *this;
    }

    inline ForwardIterator operator++(int) {
      if (IsEnd() == true) {
        return *this;
      }

      // Make a copy of the current one before advancing
      // This will increase ref count temporarily, but it is always consistent
      ForwardIterator temp = *this;

      // auto epoch_node_p = list_->epoch_manager_.JoinEpoch();
      // OperationContext ctx{epoch_node_p};

      // TODO: Add logic for handling delete node
      node_ = reinterpret_cast<SkipListInnerNode *>(GET_NEXT(node_));
      if (node_) {
        if (CHECK_DELETE(node_)) {
          node_ = reinterpret_cast<SkipListInnerNode *>(GET_NEXT(node_));
        }
        kv_p.first = &(node_->key_);
        kv_p.second = &(node_->GetValue());
      }

      // list_->epoch_manager_.LeaveEpoch(epoch_node_p);

      return temp;
    }
    // TODO: more operation need to be override
  };

  /*
   * ReversedIterator - the iterator that supports reversed traverse of the
   * SkipList
   */
  class ReversedIterator {
   private:
    KeyValuePair kv_;
    SkipListInnerNode *node_;
    SkipList *list_;
    OperationContext ctx_;

   public:
    /*
     * Constructor() - create a default iterator
     */
    ReversedIterator() : kv_(), node_(nullptr), list_(nullptr), ctx_(nullptr) {}

    /*
     * Constructor() - create a iterator that starts from the last element in
     * the SkipList
     */
    ReversedIterator(SkipList *list) : list_(list), ctx_(nullptr) {
      // LOG_INFO("RI constructed");
      this->IterJoinEpoch();

      SkipListBaseNode *cursor = this->list_->skip_list_head_.load();
      while (cursor->down_.load()) {
        while (cursor->next_.load()) {
          cursor = GET_NEXT(cursor);
        }
        cursor = cursor->down_.load();
      }
      while (cursor->next_.load()) {
        cursor = cursor->next_.load();
      }
      this->node_ = reinterpret_cast<SkipListInnerNode *>(cursor);
      if (!IsEnd()) {
        kv_.first = &(this->node_->key_);
        kv_.second = &(this->node_->GetValue());
      }
    }

    /*
     * Constructor() - create a reversed iterator starts from specific key
     */
    ReversedIterator(SkipList *list, KeyType &key)
        : list_(list), ctx_(nullptr) {
      // LOG_INFO("RI constructed with %s", key.GetInfo().c_str());

      this->IterJoinEpoch();

      auto pair = this->list_->Search(key, ctx_);
      auto next = pair.second;
      node_ = reinterpret_cast<SkipListInnerNode *>(pair.first);
      while (next) {
        if (this->list_->KeyCmpGreater(next->key_, key)) {
          break;
        }
        node_ = reinterpret_cast<SkipListInnerNode *>(next);
        next = GET_NEXT(next);
      }
      if (!IsEnd()) {
        kv_.first = &(this->node_->key_);
        kv_.second = &(this->node_->GetValue());
      }
    }

    void IterJoinEpoch() {
      auto epoch_node = this->list_->epoch_manager_.JoinEpoch();
      this->ctx_.epoch_node_ = epoch_node;
    }

    void IterLeaveEpoch() {
      this->list_->epoch_manager_.LeaveEpoch(this->ctx_.epoch_node_);
    }

    ~ReversedIterator() {}

    /*
     * return whether the iteration is end or not
     */
    bool IsEnd() { return this->node_->isHead_; }

    inline KeyValuePair *operator->() { return &kv_; }
    /*
     * the ++ in reversed iterator means accessing a node's ancestor
     */
    inline ReversedIterator &operator++() {
      if (IsEnd()) {
        return *this;
      }
      do {
        this->node_ = this->list_->SearchPrevNode(node_, this->ctx_);
      } while (CHECK_DELETE(this->node_) && !IsEnd());
      if (!IsEnd()) {
        kv_.first = &(this->node_->key_);
        kv_.second = &(this->node_->GetValue());
      }
      return *this;
    }

    inline ReversedIterator operator++(int) {
      if (IsEnd()) {
        return *this;
      }
      ReversedIterator temp = *this;

      do {
        this->node_ = this->list_->SearchPrevNode(node_, this->ctx_);
      } while (CHECK_DELETE(this->node_) && !IsEnd());
      if (!IsEnd()) {
        kv_.first = &(this->node_->key_);
        kv_.second = &(this->node_->GetValue());
      }
      return temp;
    }
    /*
     * operator*() - Return the value reference currently pointed to by this
     *               iterator
     *
     * NOTE: We need to return a constant reference to both save a value copy
     * and also to prevent caller modifying value using the reference
     */
    inline const KeyValuePair &operator*() {
      // This itself is a ValueType reference
      return kv_;
    }
  };

  /*
   * Insert() - Insert a key-value pair
   *
   * This function returns false if value already exists
   * If CAS fails this function retries until it succeeds
   */
  bool Insert(const KeyType &key, const ValueType &value) {
    auto *epoch_node_p = epoch_manager_.JoinEpoch();
    OperationContext ctx{epoch_node_p};
    bool ret = InsertNode(key, value, ctx);
    epoch_manager_.LeaveEpoch(epoch_node_p);
    return ret;
  }

  /*
   * ConditionalInsert() - Insert a key-value pair only if a given
   *                       predicate fails for all values with a key
   *
   * If return true then the value has been inserted
   * If return false then the value is not inserted. The reason could be
   * predicates returning true for one of the values of a given key
   * or because the value is already in the index
   */
  bool ConditionalInsert(
      UNUSED_ATTRIBUTE const KeyType &key,
      UNUSED_ATTRIBUTE const ValueType &value,
      UNUSED_ATTRIBUTE std::function<bool(const void *)> predicate,
      UNUSED_ATTRIBUTE bool *predicate_satisfied) {
    // LOG_INFO("ConditionalInsert Called");
    auto *epoch_node_p = epoch_manager_.JoinEpoch();
    OperationContext ctx{epoch_node_p};
    bool ret = InsertNode(key, value, ctx, predicate, predicate_satisfied);
    epoch_manager_.LeaveEpoch(epoch_node_p);
    return ret;
  }

  /*
   * Delete() - Remove a key-value pair from the tree
   *
   * This function returns false if the key and value pair does not
   * exist. Return true if delete succeeds
   */
  bool Delete(const KeyType &key, const ValueType &value) {
    // LOG_INFO("Delete called!");
    auto *epoch_node_p = epoch_manager_.JoinEpoch();
    OperationContext ctx{epoch_node_p};
    bool ret = Delete(key, value, ctx);
    epoch_manager_.LeaveEpoch(epoch_node_p);
    return ret;
  }

  /*
   * GetValue() - Fill a value list with values stored
   *
   * This function accepts a value list as argument,
   * and will copy all values into the list
   *
   * The return value is used to indicate whether the value set
   * is empty or not
   */
  bool GetValue(const KeyType &search_key, std::vector<ValueType> &value_list) {
    // LOG_INFO("GetValue()");
    auto *epoch_node_p = epoch_manager_.JoinEpoch();
    OperationContext ctx{epoch_node_p};
    bool ret = Get(search_key, value_list, ctx);
    epoch_manager_.LeaveEpoch(epoch_node_p);
    return ret;
  }

  /*
   * GetValueLimit() - Fill a limited number of value to the list.
   *
   * This function accepts a value list as argument,
   * and will copy a limited number of values start from a offset into the list
   *
   * The return value is used to indicate whether the value set
   * is empty or not
   */
  bool GetValueLimit(const KeyType &search_key,
                     std::vector<ValueType> &value_list, uint64_t limit,
                     uint64_t offset) {
    LOG_INFO("GetValue()");
    auto *epoch_node_p = epoch_manager_.JoinEpoch();
    OperationContext ctx{epoch_node_p};
    bool ret = GetLimit(search_key, value_list, limit, offset, ctx);
    epoch_manager_.LeaveEpoch(epoch_node_p);
    return ret;
  }

  ForwardIterator ForwardBegin() { return ForwardIterator{this}; }

  // returns a forward iterator from the key
  ForwardIterator ForwardBegin(KeyType &starts_key) {
    return ForwardIterator{this, starts_key};
  }

  ReversedIterator ReverseBegin() { return ReversedIterator{this}; }

  ReversedIterator ReverseBegin(KeyType &startsKey) {
    return ReversedIterator{this, startsKey};
  };

  /*
   * PerformGC() - Interface function for external users to
   *                              force a garbage collection
   */
  void PerformGC() {
    // LOG_INFO("Perform garbage collection!");
    this->epoch_manager_.PerformGarbageCollection();
    this->epoch_manager_.need_gc = false;
  }

  /*
   * NeedGC() - Whether the skiplsit needs garbage collection
   */
  bool NeedGC() {
    // LOG_INFO("Need GC!");
    return this->epoch_manager_.need_gc;
  }

  size_t GetMemoryFootprint() {
    // LOG_INFO("Get Memory Footprint!");
    return node_manager_.GetFootprint();
  }

 public:
  // Key comparator
  const KeyComparator key_cmp_obj_;

  // Raw key eq checker
  const KeyEqualityChecker key_eq_obj_;

  // Check whether values are equivalent
  const ValueEqualityChecker value_eq_obj_;

  ///////////////////////////////////////////////////////////////////
  // Key Comparison Member Functions
  ///////////////////////////////////////////////////////////////////

  /*
   * KeyCmpLess() - Compare two keys for "less than" relation
   *
   * If key1 < key2 return true
   * If not return false
   *
   * NOTE: In older version of the implementation this might be defined
   * as the comparator to wrapped key type. However wrapped key has
   * been removed from the newest implementation, and this function
   * compares KeyType specified in template argument.
   */
  inline bool KeyCmpLess(const KeyType &key1, const KeyType &key2) const {
    return key_cmp_obj_(key1, key2);
  }

  /*
   * KeyCmpEqual() - Compare a pair of keys for equality
   *
   * This functions compares keys for equality relation
   */
  inline bool KeyCmpEqual(const KeyType &key1, const KeyType &key2) const {
    return key_eq_obj_(key1, key2);
  }

  /*
   * KeyCmpGreaterEqual() - Compare a pair of keys for >= relation
   *
   * It negates result of keyCmpLess()
   */
  inline bool KeyCmpGreaterEqual(const KeyType &key1,
                                 const KeyType &key2) const {
    return !KeyCmpLess(key1, key2);
  }

  /*
   * KeyCmpGreater() - Compare a pair of keys for > relation
   *
   * It flips input for keyCmpLess()
   */
  inline bool KeyCmpGreater(const KeyType &key1, const KeyType &key2) const {
    return KeyCmpLess(key2, key1);
  }

  /*
   * KeyCmpLessEqual() - Compare a pair of keys for <= relation
   */
  inline bool KeyCmpLessEqual(const KeyType &key1, const KeyType &key2) const {
    return !KeyCmpGreater(key1, key2);
  }

  /*
   * ValueCmpEqual() - Compares whether two values are equal
   */
  inline bool ValueCmpEqual(const ValueType &v1, const ValueType &v2) {
    return value_eq_obj_(v1, v2);
  }

  /**
    * class GarbageNode - A linked list of garbage nodes
    */
  class GarbageNode {
   public:
    SkipListBaseNode *node_p;
    // Only insert at head of garbage list, no need to be atomic
    GarbageNode *next_p;
  };

  /*
   * class EpochNode - A linked list of epoch node that records thread count
   * and start garbage node
   */
  class EpochNode {
   public:
    // Track number of active threads to determine GC or not
    std::atomic<int> active_txn_count;

    // Head of garbage node list, GC nodes are CASed onto this pointer
    std::atomic<GarbageNode *> garbage_list_p;

    EpochNode *next_p;
  };
  // maintains Epoch
  // has a inside linked list in which every node represents an epoch
  class EpochManager {
   public:
    // Indicates whether destructor is running
    std::atomic<bool> destruct_flag;

    // Threshold for gc
    const static int gc_threshold = 5000;

    // Only accessed by epoch manager
    EpochNode *head_epoch_p;

    // Worker thread only read
    EpochNode *current_epoch_p;

    // Pointer to thread created by EpochManager internally
    std::thread *epoch_thread_p;

    // Need GC boolean
    bool need_gc;

    // GC Interval equals to Epoch Interval
    int GC_Interval_;

    NodeManager *node_manager_epoch_;

    EpochManager(int GC_interval, NodeManager *node_manager_)
        : GC_Interval_(GC_interval), node_manager_epoch_(node_manager_) {
      current_epoch_p = node_manager_epoch_->GetEpochNode();

      current_epoch_p->active_txn_count = 0;
      current_epoch_p->garbage_list_p = nullptr;

      head_epoch_p = current_epoch_p;

      destruct_flag.store(false);
      need_gc = false;

      GC_Interval_ = GC_interval;
    }

    ~EpochManager() {
      if (DISABLE_CLEAR_EPOCH) return;
      // LOG_INFO("footprint before deconstruct epoch:");
      this->node_manager_epoch_->PrintFootprint();
      destruct_flag.store(true);

      if (epoch_thread_p != nullptr) {
        epoch_thread_p->join();

        delete epoch_thread_p;
      }

      current_epoch_p = nullptr;

      ClearEpoch();

      if (head_epoch_p != nullptr) {
        for (EpochNode *epoch_node_p = head_epoch_p; epoch_node_p != nullptr;
             epoch_node_p = epoch_node_p->next_p) {
          epoch_node_p->active_txn_count = 0;
        }

        ClearEpoch();
      }

      PL_ASSERT(head_epoch_p == nullptr);
      // LOG_INFO("footprint after deconstruct epoch:");
      this->node_manager_epoch_->PrintFootprint();
      return;
    }

    // This function is called by worker threads, so has to consider race
    // conditions.
    void AddGarbageNode(SkipListBaseNode *node) {
      EpochNode *epoch_p = current_epoch_p;

      GarbageNode *garbage_node_p = node_manager_epoch_->GetGarbageNode();
      garbage_node_p->node_p = node;
      garbage_node_p->next_p = epoch_p->garbage_list_p.load();

      while (1) {
        bool ret = epoch_p->garbage_list_p.compare_exchange_strong(
            garbage_node_p->next_p, garbage_node_p);

        if (ret == true) {
          break;
        } else {
          LOG_TRACE("Add garbage node CAS failed. Retry");
        }
      }  // while 1
      int cur_counter = node_manager_epoch_->GetGarbageNodeCount();
      if (cur_counter > gc_threshold) {
        need_gc = true;
        // LOG_INFO("Need gc set to true!!!!");
      }
      return;
    }

    /*
     * return the current EpochNode
     * need to add the reference count of current EpochNode
     */

    EpochNode *JoinEpoch() {
    try_join_again:
      EpochNode *epoch_p = current_epoch_p;
      int64_t prev_count = epoch_p->active_txn_count.fetch_add(1);
      if (prev_count < 0) {
        epoch_p->active_txn_count.fetch_sub(1);
        goto try_join_again;
      }
      return epoch_p;
    };

    /*
     * leaves current EpochNode
     * should maintain atomicity when counting the reference
     */

    void LeaveEpoch(EpochNode *node) {
      node->active_txn_count.fetch_sub(1);
      return;
    };

    /*
     * NewEpoch() - start new epoch after the call
     *
     * begins new Epoch that caused by the
     * Need to atomically maintain the epoch list
     */
    void NewEpoch() {
      EpochNode *epoch_node_p = node_manager_epoch_->GetEpochNode();

      epoch_node_p->active_txn_count = 0;
      epoch_node_p->garbage_list_p = nullptr;

      epoch_node_p->next_p = nullptr;

      current_epoch_p->next_p = epoch_node_p;
      current_epoch_p = epoch_node_p;

      return;
    };

    /*
     * ClearEpoch() - Sweep the chain of epoch and free memory
     *
     * The minimum number of epoch we must maintain is 1 which means
     * when current epoch is the head epoch we should stop scanning
     *
     * NOTE: There is no race condition in this function since it is
     * only called by the cleaner thread
     */

    void ClearEpoch() {
      if (DISABLE_CLEAR_EPOCH) return;
      // LOG_INFO("Start to clear epoch");

      while (1) {
        if (head_epoch_p == current_epoch_p) {
          // LOG_TRACE("Current epoch is head epoch. Do not clean");
          break;
        }

        int active_txn_count = head_epoch_p->active_txn_count.load();
        PL_ASSERT(active_txn_count >= 0);

        if (active_txn_count != 0) {
          // LOG_TRACE("Head epoch is not empty. Return");
          break;
        }

        if (head_epoch_p->active_txn_count.fetch_sub(MAX_THREAD_COUNT) > 0) {
          // LOG_TRACE(
          //    "Some thread sneaks in after we have decided"
          //    " to clean. Return");

          head_epoch_p->active_txn_count.fetch_add(MAX_THREAD_COUNT);

          break;
        }

        GarbageNode *next_garbage_node_p = nullptr;

        for (GarbageNode *garbage_node_p = head_epoch_p->garbage_list_p.load();
             garbage_node_p != nullptr; garbage_node_p = next_garbage_node_p) {
          next_garbage_node_p = garbage_node_p->next_p;

          // LOG_INFO("Delete garbage node, key:%s",
          // garbage_node_p->node_p->key_.GetInfo().c_str());
          node_manager_epoch_->ReturnSkipListNode(garbage_node_p->node_p);
          node_manager_epoch_->ReturnGarbageNode(garbage_node_p);
        }  // for

        EpochNode *next_epoch_node_p = head_epoch_p->next_p;
        node_manager_epoch_->ReturnEpochNode(head_epoch_p);

        head_epoch_p = next_epoch_node_p;
      }  // while(1) through epoch nodes
      return;
    }

    void PerformGarbageCollection() {
      // LOG_INFO("Call Perform Garbage Collection!!!!!!");
      // LOG_INFO("footprint before gc: %lu",
      //         this->node_manager_epoch_->GetFootprint());
      ClearEpoch();
      // LOG_INFO("footprint after gc: %lu",
      //         this->node_manager_epoch_->GetFootprint());
      need_gc = false;
      return;
    }

    /**
    * EpochThreadFunc() - The Epoch thread executes this every GC_INTERVAL ms
    */
    void EpochThreadFunc() {
      while (!destruct_flag.load()) {
        NewEpoch();
        // if(need_gc) PerformGarbageCollection();
        // else LOG_INFO("Don't need gc, Garbage node count:%lu\n",
        // node_manager_epoch_->GetGarbageNodeCount());

        // Sleep for 50 ms
        std::chrono::microseconds duration(GC_Interval_);
        std::this_thread::sleep_for(duration);
      }

      // LOG_TRACE("Epoch manager exits; thread return");

      return;
    }

    /**
     * StartEpochThread() - Start epoch thread for garbage collection
     */
    void StartEpochThread() {
      epoch_thread_p = new std::thread{[this]() { this->EpochThreadFunc(); }};
    }
  };

  /*
   * NodeManager - maintains the SkipList Node pool
   *
   */
  class NodeManager {
   private:
    std::atomic<size_t> inner_node_count_;
    std::atomic<size_t> head_node_count_;
    std::atomic<size_t> epoch_node_count_;
    std::atomic<size_t> garbage_node_count_;

   public:
    NodeManager()
        : inner_node_count_(0),
          head_node_count_(0),
          epoch_node_count_(0),
          garbage_node_count_(0) {}
    void PrintFootprint() {
      // LOG_INFO("FootPrint: inner: %lu head %lu epoch: %lu gn: %lu",
      //         inner_node_count_.load(), head_node_count_.load(),
      //         epoch_node_count_.load(), garbage_node_count_.load());
    }
    size_t GetFootprint() {
      return sizeof(SkipListBaseNode) * head_node_count_.load() +
             sizeof(SkipListInnerNode) * inner_node_count_.load() +
             sizeof(EpochNode) * epoch_node_count_.load() +
             sizeof(GarbageNode) * garbage_node_count_.load();
    }
    size_t GetGarbageNodeCount() { return garbage_node_count_.load(); }
    size_t GetInnderNodeCount() { return inner_node_count_.load(); }
    /*
     *
     */
    SkipListBaseNode *GetSkipListHead(u_int32_t level) {
      head_node_count_.fetch_add(1);
      return new SkipListBaseNode(true, level);
    }
    /*
     * GetSkipListNode() - getSkipListNode with only key and isHead settled
     */
    SkipListBaseNode *GetSkipListNode(KeyType key, bool isHead,
                                      u_int32_t level) {
      return new SkipListBaseNode(nullptr, nullptr, nullptr, key, isHead,
                                  level);
    }
    /*
     * GetSkipListNode() - get SkipListNode full equiped
     */
    SkipListBaseNode *GetSkipListNode(SkipListBaseNode *next,
                                      SkipListBaseNode *down,
                                      SkipListBaseNode *back_link, KeyType key,
                                      bool isHead, u_int32_t level) {
      return new SkipListBaseNode(next, down, back_link, key, isHead, level);
    }
    /*
     * GetSkipListInnerNode() - get a SkipListInnerNode using key and value
     */
    SkipListInnerNode *GetSkipListInnerNode(KeyType key, ValueType value,
                                            u_int32_t level) {
      inner_node_count_.fetch_add(1);
      auto tmp =
          new SkipListInnerNode(nullptr, nullptr, nullptr, key, false, level);
      tmp->SetValue(value);
      return tmp;
    }
    /*
     * GetSkipListInnerNode() - get a SkipListInnerNode using key and root
     * pointer
     */
    SkipListInnerNode *GetSkipListInnerNode(KeyType key,
                                            SkipListInnerNode *root,
                                            SkipListInnerNode *down,
                                            u_int32_t level) {
      inner_node_count_.fetch_add(1);
      auto tmp =
          new SkipListInnerNode(nullptr, down, nullptr, key, false, level);
      tmp->SetRoot(root);
      return tmp;
    }
    /*
     * GetEpochNode() - get an EpochNode
     */
    EpochNode *GetEpochNode() {
      epoch_node_count_.fetch_add(1);
      auto tmp = new EpochNode();
      return tmp;
    }
    /*
     * GetGarbageNode() - get a garbage node
     */
    GarbageNode *GetGarbageNode() {
      garbage_node_count_.fetch_add(1);
      auto tmp = new GarbageNode();
      return tmp;
    }
    void ReturnSkipListNode(SkipListBaseNode *node) {
      if (node->isHead_) {
        head_node_count_.fetch_sub(1);
        // LOG_INFO("Deleted Head Node, is it right time??");
      } else
        inner_node_count_.fetch_sub(1);
      delete node;
    }
    void ReturnEpochNode(UNUSED_ATTRIBUTE EpochNode *node) {
      epoch_node_count_.fetch_sub(1);
      delete node;
    }
    void ReturnGarbageNode(UNUSED_ATTRIBUTE GarbageNode *node) {
      garbage_node_count_.fetch_sub(1);
      delete node;
    }
  };

  /*
   * OperationContext - maintains info and context of each thread
   *
   * EpochNode: the epoch node that the thread is in
   */
  class OperationContext {
   public:
    EpochNode *epoch_node_;
    OperationContext(EpochNode *epoch_node) : epoch_node_(epoch_node) {}
  };
};
}  // namespace index
}  // namespace peloton

#endif
