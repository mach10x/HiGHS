/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2023 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"

#include <algorithm>
/*
#include <atomic>
#include <cmath>
#include <limits>

#include "Highs.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsSolution.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsImplications.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsObjectiveFunction.h"
#include "pdqsort/pdqsort.h"
#include "presolve/HighsPostsolveStack.h"
#include "test/DevKkt.h"
#include "util/HFactor.h"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"
#include "util/HighsLinearSumBounds.h"
#include "util/HighsSplay.h"
#include "util/HighsUtils.h"

#define ENABLE_SPARSIFY_FOR_LP 0

#define HPRESOLVE_CHECKED_CALL(presolveCall)                           \
  do {                                                                 \
    HPresolve::Result __result = presolveCall;                         \
    if (__result != presolve::HPresolve::Result::kOk) return __result; \
  } while (0)
*/
namespace presolve {

#ifndef NDEBUG
void HPresolve::debugPrintRow(HighsPostsolveStack& postsolve_stack,
                              HighsInt row) {
  printf("(row %" HIGHSINT_FORMAT ") %.15g (impl: %.15g) <= ",
         postsolve_stack.getOrigRowIndex(row), model->row_lower_[row],
         impliedRowBounds.getSumLower(row));

  for (const HighsSliceNonzero& nonzero : getSortedRowVector(row)) {
    // for (HighsInt rowiter = rowhead[row]; rowiter != -1; rowiter =
    // ARnext[rowiter]) {
    char colchar =
        model->integrality_[nonzero.index()] == HighsVarType::kInteger ? 'y'
                                                                       : 'x';
    char signchar = nonzero.value() < 0 ? '-' : '+';
    printf("%c%g %c%" HIGHSINT_FORMAT " ", signchar, std::abs(nonzero.value()),
           colchar, postsolve_stack.getOrigColIndex(nonzero.index()));
  }

  printf("<= %.15g (impl: %.15g)\n", model->row_upper_[row],
         impliedRowBounds.getSumUpper(row));
}

void HPresolve::debugReportRowWrtCol(const HighsPostsolveStack& postsolve_stack,
				     const HighsInt row, const HighsInt col,
				     const std::string message) {
  const HighsInt original_col_index = postsolve_stack.getOrigColIndex(col);
  const HighsInt original_row_index = postsolve_stack.getOrigRowIndex(row);
  if (original_col_index < 0 || original_row_index < 0) return;
  printf("\nFor Row %d (originally %d) WRT col %d (originally %d) %s\n",
	 int(row), int(original_row_index),
	 int(col), int(original_col_index),
	 message.c_str());
  printf("(row %" HIGHSINT_FORMAT ") %11.4g (impl: %11.4g) <= ",
         postsolve_stack.getOrigRowIndex(row), model->row_lower_[row],
         impliedRowBounds.getSumLower(row));

  for (const HighsSliceNonzero& nonzero : getSortedRowVector(row)) {
    char colchar =
        model->integrality_[nonzero.index()] == HighsVarType::kInteger ? 'y'
                                                                       : 'x';
    char signchar = nonzero.value() < 0 ? '-' : '+';
    printf("%c%g %c%d ", signchar, std::abs(nonzero.value()),
           colchar, int(postsolve_stack.getOrigColIndex(nonzero.index())));
    
  }
  printf("<= %11.4g (impl: %11.4g)\n", model->row_upper_[row],
         impliedRowBounds.getSumUpper(row));
  double row_min_activity_wo_col = 0;
  double row_max_activity_wo_col = 0;
  double col_val = 0;
  for (const HighsSliceNonzero& nonzero : getSortedRowVector(row)) {
    HighsInt iCol = nonzero.index();
    printf("Col %2d bounds [%11.4g, %11.4g] implied [%11.4g, %11.4g] source [%2d, %2d] activity [%11.4g, %11.4g]\n",
	   int(iCol),
	   model->col_lower_[iCol], model->col_upper_[iCol],
	   implColLower[iCol], implColUpper[iCol],
	   int(colLowerSource[iCol]), int(colUpperSource[iCol]),
	   row_min_activity_wo_col, row_max_activity_wo_col);
    if (col != iCol) {
      if (nonzero.value() > 0) {
	row_min_activity_wo_col += nonzero.value() * model->col_lower_[iCol];
	row_max_activity_wo_col += nonzero.value() * model->col_upper_[iCol];
      } else {
	row_min_activity_wo_col += nonzero.value() * model->col_upper_[iCol];
	row_max_activity_wo_col += nonzero.value() * model->col_lower_[iCol];
      }
    } else {
      col_val = nonzero.value();
    }
  }
  if (col_val) {
    double local_implied_col_lb = col_val>0 ?
      (model->row_lower_[row]-row_max_activity_wo_col)/col_val :
      (model->row_upper_[row]-row_min_activity_wo_col)/col_val;
    double local_implied_col_ub = col_val>0 ?
      (model->row_upper_[row]-row_min_activity_wo_col)/col_val :
      (model->row_lower_[row]-row_max_activity_wo_col)/col_val;
    printf("Row actvity without col is in [%11.4g, %11.4g] Local implied col bounds [%11.4g, %11.4g]\n",
	   row_min_activity_wo_col, row_max_activity_wo_col,
	   local_implied_col_lb, local_implied_col_ub);
  }
  printf("\n");
}

bool HPresolve::debugColImpliedBoundsNotUpToDate(HighsInt row, HighsInt col,
                                                 double val) {
  // for debugging!  implied upper bound on variable x_6 (computed
  // using row 1) is incorrect when entering this method.  as a
  // consequence, this subroutine incorrectly decides that the bounds
  // on x_6 (singleton!) are redundant and computes incorrect row
  // duals.  Why are 'implColLower' and 'implColUpper' not re-computed
  // when the corresponding rows stored in 'colLowerSource' or
  // 'colUpperSource' are modified (e.g., non-zeros added by
  // substitution or constraint bounds modified)?
  const HighsInt check_col = 6; 
  const HighsInt check_row = 1;
  bool notUpToDate = false;
  if (col != check_col || row != check_row) return notUpToDate;
  printf("debugColImpliedBoundsNotUpToDate(%d, %d, %g)\n",
	 int(row), int(col), val);
  assert(colsize[col] == 1);
  bool oldColLowerUpToDate = true;
  bool oldColUpperUpToDate = true;
  // save old implied bounds
  double oldImplLower = implColLower[col];
  double oldImplUpper = implColUpper[col];
  HighsInt oldImplLowerSource = colLowerSource[col];
  HighsInt oldImplUpperSource = colUpperSource[col];

  if ((colLowerSource[col] == row) || (colUpperSource[col] == row)) {
    // set implied bounds to infinity
    if (colLowerSource[col] == row)
      changeImplColLower(col, -kHighsInf, -1);
    if (colUpperSource[col] == row)
      changeImplColUpper(col, kHighsInf, -1);
    // recompute implied bounds
    updateColImpliedBounds(row, col, val);
    // check if bounds were correct / up-to-date when entering this method
    const double dl_lower = std::fabs(oldImplLower - implColLower[col]);
    const double dl_upper = std::fabs(oldImplUpper - implColUpper[col]);
    oldColLowerUpToDate = 
      ((oldImplLower == -kHighsInf) && (implColLower[col] == -kHighsInf)) ||
      (dl_lower <= primal_feastol);
    oldColUpperUpToDate =
      ((oldImplUpper == kHighsInf) && (implColUpper[col] == kHighsInf)) ||
      (dl_upper <= primal_feastol);
    notUpToDate = (!oldColLowerUpToDate) || (!oldColUpperUpToDate);
    if (notUpToDate) {
      printf("ColImpliedBoundsNotUpToDate for row = %d; col = %d: lower(%g, %g, %g); upper(%g, %g, %g)\n",
	     int(row), int(col),
	     oldImplLower, implColLower[col], dl_lower,
	     oldImplUpper, implColUpper[col], dl_upper);
    }
  }

  return notUpToDate;
}

bool HPresolve::debugCheckColImpliedBoundsOk(HighsPostsolveStack& postsolve_stack) {
  const double implied_bound_error_tolerance = 1e-6;
  std::vector<HighsInt> check_colLowerSource(model->num_col_, -1);
  std::vector<HighsInt> check_colUpperSource(model->num_col_, -1);
  std::vector<double> check_implColLower(model->num_col_, -kHighsInf);
  std::vector<double> check_implColUpper(model->num_col_, kHighsInf);
  std::vector<double> row_min_activity(model->num_row_, 0);
  std::vector<double> row_max_activity(model->num_row_, 0);
  for (HighsInt iRow = 0; iRow < model->num_row_; iRow++) {
    for (const HighsSliceNonzero& nonzero : getSortedRowVector(iRow)) {
      HighsInt iCol = nonzero.index();
      assert(nonzero.value());
      if (nonzero.value() > 0) {
	row_min_activity[iRow] += nonzero.value() * model->col_lower_[iCol];
	row_max_activity[iRow] += nonzero.value() * model->col_upper_[iCol];
      } else {
	row_min_activity[iRow] += nonzero.value() * model->col_upper_[iCol];
	row_max_activity[iRow] += nonzero.value() * model->col_lower_[iCol];
      }
    }
  }
  for (HighsInt iCol = 0; iCol < model->num_col_; iCol++) {
    for (HighsInt iEl = model->a_matrix_.start_[iCol]; iEl < model->a_matrix_.start_[iCol+1]; iEl++) {
      HighsInt iRow = model->a_matrix_.index_[iEl];
      double col_val = model->a_matrix_.value_[iEl];  
      double row_min_activity_wo_col = row_min_activity[iRow];
      double row_max_activity_wo_col = row_max_activity[iRow];
      if (col_val > 0) {
	row_min_activity_wo_col -= col_val * model->col_lower_[iCol];
	row_max_activity_wo_col -= col_val * model->col_upper_[iCol];
      } else {
	row_min_activity_wo_col -= col_val * model->col_upper_[iCol];
	row_max_activity_wo_col -= col_val * model->col_lower_[iCol];
      }
      double local_implied_col_lower = col_val>0 ?
	(model->row_lower_[iRow]-row_max_activity_wo_col)/col_val :
	(model->row_upper_[iRow]-row_min_activity_wo_col)/col_val;
      if (check_implColLower[iCol] < local_implied_col_lower) {
	check_implColLower[iCol] = local_implied_col_lower;
	check_colLowerSource[iCol] = iRow;
      }
      double local_implied_col_upper = col_val>0 ?
	(model->row_upper_[iRow]-row_min_activity_wo_col)/col_val :
	(model->row_lower_[iRow]-row_max_activity_wo_col)/col_val;
      if (check_implColUpper[iCol] > local_implied_col_upper) {
	check_implColUpper[iCol] = local_implied_col_upper;
	check_colUpperSource[iCol] = iRow;
      }
    }
  }
  double norm_implied_col_lower_er = 0;
  double norm_implied_col_upper_er = 0;
  for (HighsInt iCol = 0; iCol < model->num_col_; iCol++) {
    double implied_col_lower_er = (check_implColLower[iCol] == -kHighsInf &&
				   implColLower[iCol] == -kHighsInf) ? 0 :
      std::fabs(check_implColLower[iCol]-implColLower[iCol]);
    double implied_col_upper_er = (check_implColUpper[iCol] == kHighsInf &&
				   implColUpper[iCol] == kHighsInf)  ? 0 :
      std::fabs(check_implColUpper[iCol]-implColUpper[iCol]);
    norm_implied_col_lower_er = std::max(implied_col_lower_er, norm_implied_col_lower_er);
    norm_implied_col_upper_er = std::max(implied_col_upper_er, norm_implied_col_upper_er);
    if (implied_col_lower_er > implied_bound_error_tolerance ||
	implied_col_upper_er > implied_bound_error_tolerance) {
      printf("Column %4d has implied bounds [%11.4g, %11.4g] check implied bounds [%11.4g, %11.4g] errors [%11.4g, %11.4g]\n",
	     int(iCol),
	     implColLower[iCol], implColUpper[iCol], 
	     check_implColLower[iCol], check_implColUpper[iCol], 
	     implied_col_lower_er, implied_col_upper_er);
    }
  }
  if (norm_implied_col_lower_er || norm_implied_col_upper_er) {
    printf("Columns have ||implied_lower_bound_error|| = %11.4g; ||implied_upper_bound_error|| = %11.4g\n",
	   norm_implied_col_lower_er, norm_implied_col_upper_er);
   return false;
  }
  return true;
}

#endif
}  // namespace presolve
