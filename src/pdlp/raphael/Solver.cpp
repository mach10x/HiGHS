/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/raphael/Solver.cpp
 */
#include "pdlp/raphael/Solver.h"
#include "lp_data/HighsLpUtils.h"

HighsStatus solveLpRaphael(HighsLpSolverObject& solver_object) {
  return solveLpRaphael(solver_object.options_, solver_object.timer_,
			solver_object.lp_, solver_object.basis_,
			solver_object.solution_, solver_object.model_status_,
			solver_object.highs_info_, solver_object.callback_);
}

HighsStatus solveLpRaphael(const HighsOptions& options, HighsTimer& timer,
			   const HighsLp& lp, HighsBasis& highs_basis,
			   HighsSolution& highs_solution,
			   HighsModelStatus& model_status, HighsInfo& highs_info,
			   HighsCallback& callback) {
  // Indicate that there is no valid primal solution, dual solution or basis
  highs_basis.valid = false;
  highs_solution.value_valid = false;
  highs_solution.dual_valid = false;
  // Indicate that no imprecise solution has (yet) been found
  resetModelStatusAndHighsInfo(model_status, highs_info);

  double standard_form_offset;
  std::vector<double> standard_form_cost;
  std::vector<double> standard_form_rhs;
  HighsSparseMatrix standard_form_matrix;

  formStandardFormLp(lp,
		     options.log_options,
		     standard_form_offset,
		     standard_form_cost,
		     standard_form_rhs,
		     standard_form_matrix);
  const HighsInt num_col = standard_form_cost.size();
  const HighsInt num_row = standard_form_rhs.size();
  const HighsInt num_nz = standard_form_matrix.numNz();
  printf("Standard form LP has %d columns, %d rows and %d nonzeros\n",
	 int(num_col), int(num_row), int(num_nz)); 

  // Now solve the LP in standard form using PDLP

  // Once solved, the solution for the LP in standard form obtained
  // with PDLP needs to be converted to a solution to the original
  // LP. Do this with a call in the following line to be written by
  // Julian

  // For the moment, return the model status as kSolveError, and HiGHS
  // status as error, so HiGHS doesn't expect anything in terms of a
  // primal or dual solution
  model_status = HighsModelStatus::kSolveError;
  return HighsStatus::kError;
}
