//===----------------------------------------------------------------------===//
//
//                         Peloton
//
// plan_comparator.h
//
// Identification: src/include/codegen/plan_comparator.h
//
// Copyright (c) 2015-17, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#pragma once

#include "planner/abstract_plan.h"
#include "planner/abstract_scan_plan.h"
#include "planner/aggregate_plan.h"
#include "planner/delete_plan.h"
#include "planner/hash_join_plan.h"
#include "planner/hash_plan.h"
#include "planner/insert_plan.h"
#include "planner/projection_plan.h"
#include "planner/order_by_plan.h"
#include "planner/seq_scan_plan.h"
#include "planner/update_plan.h"

namespace peloton {
namespace codegen {
//===----------------------------------------------------------------------===//
// PlanComparator provides compare function for all plans that are
// supported by codegen including SeqScanPlan, OrderByPlan, etc.
//===----------------------------------------------------------------------===//
class PlanComparator {
 public:
  // No constructor
  PlanComparator() = delete;

  static int Compare(const planner::AbstractPlan &,
                     const planner::AbstractPlan &);
 private:
  // Compare functions
  static int CompareAggregate(const planner::AggregatePlan &,
                              const planner::AggregatePlan &);
  static int CompareDelete(const planner::DeletePlan &,
                           const planner::DeletePlan &);
  static int CompareHash(const planner::HashPlan &, const planner::HashPlan &);
  static int CompareHashJoin(const planner::HashJoinPlan &,
                             const planner::HashJoinPlan &);
  static int CompareInsert(const planner::InsertPlan &,
                           const planner::InsertPlan &);
  static int CompareOrderBy(const planner::OrderByPlan &,
                            const planner::OrderByPlan &);
  static int CompareProjection(const planner::ProjectionPlan &,
                               const planner::ProjectionPlan &);
  static int CompareSeqScan(const planner::SeqScanPlan &,
                            const planner::SeqScanPlan &);
  static int CompareUpdate(const planner::UpdatePlan &,
                           const planner::UpdatePlan &);

  // Helper functions for comparison
  static int CompareAggType(
      const std::vector<planner::AggregatePlan::AggTerm> &,
      const std::vector<planner::AggregatePlan::AggTerm> &);
  static int CompareChildren(const planner::AbstractPlan &,
                             const planner::AbstractPlan &);
  static int CompareDerivedAttr(const planner::DerivedAttribute &,
                                const planner::DerivedAttribute &);
  static int CompareProjectInfo(const planner::ProjectInfo *,
                                const planner::ProjectInfo *);
  static int CompareSchema(const catalog::Schema &, const catalog::Schema &);
};

//===----------------------------------------------------------------------===//
// ExpressionComparator provides compare function for all expressions that are
// supported by codegen.
//===----------------------------------------------------------------------===//
class ExpressionComparator {
 public:
  // No constructor
  ExpressionComparator() = delete;

  static int Compare(const expression::AbstractExpression *,
                     const expression::AbstractExpression *);
 private:
  static int CompareChildren(const expression::AbstractExpression *,
                             const expression::AbstractExpression *);
};

}  // namespace codegen
}  // namespace peloton