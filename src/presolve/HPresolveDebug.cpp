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
#include <algorithm>

#include "presolve/HPresolve.h"
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
         int(row), int(original_row_index), int(col), int(original_col_index),
         message.c_str());
  printf("(row %" HIGHSINT_FORMAT ") %11.4g (impl: %11.4g) <= ",
         postsolve_stack.getOrigRowIndex(row), model->row_lower_[row],
         impliedRowBounds.getSumLower(row));

  for (const HighsSliceNonzero& nonzero : getSortedRowVector(row)) {
    char colchar =
        model->integrality_[nonzero.index()] == HighsVarType::kInteger ? 'y'
                                                                       : 'x';
    char signchar = nonzero.value() < 0 ? '-' : '+';
    printf("%c%g %c%d ", signchar, std::abs(nonzero.value()), colchar,
           int(postsolve_stack.getOrigColIndex(nonzero.index())));
  }
  printf("<= %11.4g (impl: %11.4g)\n", model->row_upper_[row],
         impliedRowBounds.getSumUpper(row));
  double row_min_activity_wo_col = 0;
  double row_max_activity_wo_col = 0;
  double col_val = 0;
  for (const HighsSliceNonzero& nonzero : getSortedRowVector(row)) {
    HighsInt iCol = nonzero.index();
    printf(
        "Col %2d bounds [%11.4g, %11.4g] implied [%11.4g, %11.4g] source [%2d, "
        "%2d] activity [%11.4g, %11.4g]\n",
        int(iCol), model->col_lower_[iCol], model->col_upper_[iCol],
        implColLower[iCol], implColUpper[iCol], int(colLowerSource[iCol]),
        int(colUpperSource[iCol]), row_min_activity_wo_col,
        row_max_activity_wo_col);
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
    double local_implied_col_lb =
        col_val > 0
            ? (model->row_lower_[row] - row_max_activity_wo_col) / col_val
            : (model->row_upper_[row] - row_min_activity_wo_col) / col_val;
    double local_implied_col_ub =
        col_val > 0
            ? (model->row_upper_[row] - row_min_activity_wo_col) / col_val
            : (model->row_lower_[row] - row_max_activity_wo_col) / col_val;
    printf(
        "Row actvity without col is in [%11.4g, %11.4g] Local implied col "
        "bounds [%11.4g, %11.4g]\n",
        row_min_activity_wo_col, row_max_activity_wo_col, local_implied_col_lb,
        local_implied_col_ub);
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
  printf("debugColImpliedBoundsNotUpToDate(%d, %d, %g)\n", int(row), int(col),
         val);
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
    if (colLowerSource[col] == row) changeImplColLower(col, -kHighsInf, -1);
    if (colUpperSource[col] == row) changeImplColUpper(col, kHighsInf, -1);
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
      printf(
          "ColImpliedBoundsNotUpToDate for row = %d; col = %d: lower(%g, %g, "
          "%g); upper(%g, %g, %g)\n",
          int(row), int(col), oldImplLower, implColLower[col], dl_lower,
          oldImplUpper, implColUpper[col], dl_upper);
    }
  }

  return notUpToDate;
}

#endif
}  // namespace presolve
