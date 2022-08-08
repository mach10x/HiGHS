
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkSifting.cpp
 * @brief
 */
#include "io/HMPSIO.h"
#include "lp_data/HighsLpSolverObject.h"
#include "simplex/HEkk.h"
#include "util/HighsRandom.h"
#include "util/HighsSort.h"

const bool check_duals = true;
const double purge_threshold = 1.5;
const double purge_multiplier = 0.8;

HighsStatus HEkk::sifting() {
  HighsStatus return_status = HighsStatus::kOk;
  // Need to start from a primal feasible solution
  assert(info_.num_primal_infeasibilities >= 0);

  HighsInt sifted_list_max_count = lp_.num_row_;

  std::vector<HighsInt> sifted_list;
  std::vector<bool> in_sifted_list;
  in_sifted_list.assign(lp_.num_col_, false);
  HEkk sifted_ekk_instance;
  HighsLp sifted_lp;
  HighsBasis sifted_basis;
  HighsSolution sifted_solution;
  HighsInfo sifted_highs_info;
  HighsOptions sifted_options = *options_;
  HighsLpSolverObject sifted_solver_object(
      sifted_lp, sifted_basis, sifted_solution, sifted_highs_info,
      sifted_ekk_instance, sifted_options, *timer_);

  // Prevent recursive sifting!
  sifted_options.sifting_strategy = kSiftingStrategyOff;
  //  sifted_options.log_dev_level = 3;
  //  sifted_options.highs_analysis_level = 4;
  SimplexAlgorithm last_algorithm = SimplexAlgorithm::kNone;

  getInitialRowStatusAndDual(sifted_solver_object);

  HighsInt sifting_iter = 0;
  bool first_sifted_lp = true;
  this->called_return_from_solve_ = false;
  for (;;) {
    HighsInt num_add_to_sifted_list =
        addToSiftedList(lp_.num_row_, sifted_solver_object, sifted_list,
                        in_sifted_list, first_sifted_lp);
    if (num_add_to_sifted_list == 0) {
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   "Optimal after %d sifting iterations\n", (int)sifting_iter);
      model_status_ = HighsModelStatus::kOptimal;
      exit_algorithm_ = last_algorithm;

      getSiftedBasisSolutionInvertWeightsInfo(sifted_solver_object,
                                              sifted_list);
      return returnFromSolve(HighsStatus::kOk);
    }
    assert(debugOkSiftedList(sifted_list, in_sifted_list));
    first_sifted_lp = false;
    sifting_iter++;
    const bool write_lp = false;
    if (write_lp) {
      HighsModel model;
      model.lp_ = sifted_solver_object.lp_;
      writeModelAsMps(*options_, "sifted.mps", model);
    }
    return_status = sifted_ekk_instance.solve();
    assert(return_status == HighsStatus::kOk);

    last_algorithm = sifted_ekk_instance.exit_algorithm_;
    getSiftedLpInfo(sifted_solver_object);

    assert(sifted_ekk_instance.debugRetainedDataOk(sifted_ekk_instance.lp_) !=
           HighsDebugStatus::kLogicalError);

    highsLogUser(options_->log_options, HighsLogType::kInfo,
                 "Sifting iter %3d: LP has %6d rows and %6d columns. "
                 "After %6d simplex iters: %6d primal and %6d dual infeas; obj "
                 "= %15.8g\n",
                 int(sifting_iter), int(lp_.num_row_), int(sifted_list.size()),
                 int(iteration_count_),
                 sifted_ekk_instance.info_.primal_objective_value,
                 int(info_.num_primal_infeasibilities),
                 int(info_.num_dual_infeasibilities));
    if (sifting_iter > 1000) break;

    const HighsInt sifted_lp_num_col = sifted_list.size();
    if (sifted_lp_num_col < purge_threshold * lp_.num_row_) continue;
    // Purge some entries
    const HighsInt purge_num_col = purge_multiplier * lp_.num_row_;
    
    highsLogUser(options_->log_options, HighsLogType::kInfo,
                 "Sifting iter %3d: Purge %6d of %6d columns from the LP\n",
                 int(sifting_iter), int(purge_num_col), int(sifted_lp_num_col));
    purgeSiftedList(purge_num_col, sifted_solver_object,
		    sifted_list, in_sifted_list);
    
  }

  assert(1 == 0);
  return HighsStatus::kError;
}

void HEkk::purgeSiftedList(const HighsInt num_purge_from_sifted_list,
			   HighsLpSolverObject& sifted_solver_object,
			   std::vector<HighsInt>& sifted_list,
			   std::vector<bool>& in_sifted_list) {
  const HighsInt sifted_lp_num_col = sifted_list.size();
  const HighsInt sifted_lp_num_row = lp_.num_row_;
  HEkk& sifted_ekk_instance = sifted_solver_object.ekk_instance_;
  std::vector<double>& workDual =
      sifted_solver_object.ekk_instance_.info_.workDual_;
  std::vector<HighsInt> heap_index;
  std::vector<double> heap_value;
  heap_index.push_back(0);
  heap_value.push_back(0);
  for (HighsInt iX = 0; iX < sifted_lp_num_col; iX++) {
    HighsInt iCol = sifted_list[iX];
    double dual = workDual[iX];
    // Determine the dual feasibility for this column
    const double lower = info_.workLower_[iCol];
    const double upper = info_.workUpper_[iCol];
    assert(lower<upper);
    double dual_feasibility = 0;
    if (lower > -kHighsInf || upper < kHighsInf) {
      // Not free: dual feasibility is given by the dual value signed
      // by nonbasicMove since there must be no equalities
      dual_feasibility = basis_.nonbasicMove_[iCol] * dual;
    }
    assert(dual_feasibility >= -options_->dual_feasibility_tolerance);
    heap_index.push_back(iX);
    heap_value.push_back(-dual_feasibility);
  }  
  // Sort candidates to leave the sifted list
  HighsInt heap_num_en = heap_index.size() - 1;
  assert(heap_num_en == sifted_lp_num_col);
  // Sort by increasing dual merit - since decreasing sort isn't
  // available.
  maxheapsort(&heap_value[0], &heap_index[0], heap_num_en);
  std::vector<HighsInt> mask;
  mask.assign(sifted_lp_num_col, 0);
  for (HighsInt k = 1; k <= num_purge_from_sifted_list; k++)
    mask[heap_index[k]] = 1;
  HighsInt sifted_lp_new_num_col = 0;
  for (HighsInt iX = 0; iX < sifted_lp_num_col; iX++) {
    HighsInt iCol = sifted_list[iX];
    if (mask[iX]) {
      // Purge
      in_sifted_list[iCol] = false;
      // Need to ensure that the main dual value is up-to-date, since
      // they are used to consider adding the variable to the sifted
      // LP
      this->info_.workDual_[iCol] = workDual[iX];
    } else {
      // Keep
      sifted_list[sifted_lp_new_num_col] = iCol;
      // Shift workDual, as it's needed for to check duals in
      // HEkk::addToSiftedList
      if (check_duals)
	workDual[sifted_lp_new_num_col] = workDual[iX];
      sifted_lp_new_num_col++;
    }
  }
  // Shift the row duals as they're needed for HEkk::addToSiftedList
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) 
    workDual[sifted_lp_new_num_col+iRow] = workDual[sifted_lp_num_col+iRow];
  workDual.resize(sifted_lp_new_num_col + lp_.num_row_);
  HighsIndexCollection delete_sifted_cols;
  create(delete_sifted_cols, &mask[0], sifted_lp_num_col);
  sifted_ekk_instance.deleteNonbasicColsFromLp(delete_sifted_cols);
  sifted_list.resize(sifted_lp_new_num_col);
}

void HEkk::getInitialRowStatusAndDual(
    HighsLpSolverObject& sifted_solver_object) {
  // Set up the initial basis status for rows
  assert(status_.has_basis);
  HighsBasis& sifted_basis = sifted_solver_object.basis_;
  sifted_basis.row_status.resize(lp_.num_row_);
  std::vector<double>& workDual =
      sifted_solver_object.ekk_instance_.info_.workDual_;
  workDual.resize(lp_.num_row_);

  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iVar = lp_.num_col_ + iRow;
    if (basis_.nonbasicFlag_[iVar]) {
      if (basis_.nonbasicMove_[iVar] > 0) {
        sifted_basis.row_status[iRow] = HighsBasisStatus::kUpper;
      } else if (basis_.nonbasicMove_[iVar] < 0) {
        sifted_basis.row_status[iRow] = HighsBasisStatus::kLower;
      } else if (info_.workLower_[iVar] == info_.workUpper_[iVar]) {
        sifted_basis.row_status[iRow] = HighsBasisStatus::kUpper;
      } else {
        sifted_basis.row_status[iRow] = HighsBasisStatus::kZero;
      }
    } else {
      sifted_basis.row_status[iRow] = HighsBasisStatus::kBasic;
    }
    workDual[iRow] = info_.workDual_[iVar];
  }
}

HighsInt HEkk::addToSiftedList(const HighsInt max_add_to_sifted_list,
                               HighsLpSolverObject& sifted_solver_object,
                               std::vector<HighsInt>& sifted_list,
                               std::vector<bool>& in_sifted_list,
                               const bool first_sifted_lp) {
  HighsLp& sifted_lp = sifted_solver_object.lp_;
  HighsBasis& sifted_basis = sifted_solver_object.basis_;
  HEkk& sifted_ekk_instance = sifted_solver_object.ekk_instance_;
  HighsInt sifted_lp_num_col = sifted_list.size();
  const HighsInt sifted_lp_num_row = lp_.num_row_;
  const bool primal_feasible = info_.num_primal_infeasibilities == 0;
  std::vector<double>& workDual = sifted_ekk_instance.info_.workDual_;
  HighsInt& num_dual_infeasibilities = info_.num_dual_infeasibilities;
  if (!primal_feasible) assert(first_sifted_lp);
  if (first_sifted_lp) {
    assert(sifted_basis.col_status.size() == 0);
    assert(sifted_lp.num_col_ == 0);
    assert(sifted_lp.num_row_ == 0);
    assert(sifted_lp.col_cost_.size() == 0);
    assert(sifted_lp.col_lower_.size() == 0);
    assert(sifted_lp.col_upper_.size() == 0);
    assert(sifted_lp.row_lower_.size() == 0);
    assert(sifted_lp.row_upper_.size() == 0);
  }
  if (check_duals && primal_feasible) {
    // Check the computation of column duals by re-computing the duals
    // for sifted columns
    for (HighsInt iX = 0; iX < sifted_list.size(); iX++) {
      HighsInt iCol = sifted_list[iX];
      double dual = info_.workCost_[iCol];
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol];
           iEl < lp_.a_matrix_.start_[iCol + 1]; iEl++) {
        HighsInt iRow = lp_.a_matrix_.index_[iEl];
        dual += lp_.a_matrix_.value_[iEl] * workDual[sifted_lp_num_col + iRow];
      }
      const bool dual_ok = std::fabs(dual - workDual[iX]) < 1e-4;
      if (!dual_ok) {
	printf("iX %2d; iCol %4d: dual = %11.4g; workDual[iX] =  %11.4g\n",
	       int(iX), int(iCol), dual, workDual[iX]);
	fflush(stdout);
      }
      assert(dual_ok);
    }
  }

  HighsInt num_add_to_sifted_list = 0;
  std::vector<HighsInt> heap_index;
  std::vector<double> heap_value;
  heap_index.push_back(0);
  heap_value.push_back(0);
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    if (in_sifted_list[iCol]) continue;
    // Not in the sifted list, so consider adding this column
    //
    // A basic column must be added - but this should only happen when
    // setting up the first sifted LP
    //
    // For a nonbasic column
    //
    // ... if the basis is primal feasible and the column is dual
    // infeasible, add the index and dual merit to the heap
    //
    // ... if the basis is not primal feasible, the sifted list is a
    // random selection of columns constructed in a later loop
    if (basis_.nonbasicFlag_[iCol] == 0) {
      // Basic, so this must be the first LP and column must be added
      // to the sifted list.
      assert(1 == 0);
      assert(first_sifted_lp);
      num_add_to_sifted_list++;
      sifted_list.push_back(iCol);
      in_sifted_list[iCol] = true;
      sifted_lp.col_cost_.push_back(info_.workCost_[iCol]);
      sifted_lp.col_lower_.push_back(info_.workLower_[iCol]);
      sifted_lp.col_upper_.push_back(info_.workUpper_[iCol]);
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol];
           iEl < lp_.a_matrix_.start_[iCol + 1]; iEl++) {
        sifted_lp.a_matrix_.index_.push_back(lp_.a_matrix_.index_[iEl]);
        sifted_lp.a_matrix_.value_.push_back(lp_.a_matrix_.value_[iEl]);
      }
      sifted_lp.a_matrix_.start_.push_back(sifted_lp.a_matrix_.index_.size());
      sifted_basis.col_status.push_back(HighsBasisStatus::kBasic);
      continue;
    }
    // If not primal feasible, the sifted list is a random selection
    // of columns constructed in a later loop
    if (!primal_feasible) continue;
    // Nonbasic, and not in sifted list, so possible new entry
    //
    // Dual is known for first sifted LP
    double dual = info_.workDual_[iCol];
    if (check_duals || !first_sifted_lp) {
      // Compute the dual for this column
      dual = info_.workCost_[iCol];
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol];
           iEl < lp_.a_matrix_.start_[iCol + 1]; iEl++) {
        HighsInt iRow = lp_.a_matrix_.index_[iEl];
        dual += lp_.a_matrix_.value_[iEl] * workDual[sifted_lp_num_col + iRow];
      }
      if (first_sifted_lp)
        assert(std::fabs(dual - info_.workDual_[iCol]) < 1e-4);
      info_.workDual_[iCol] = dual;
    }
    // Determine the dual infeasibility for this column
    const double lower = info_.workLower_[iCol];
    const double upper = info_.workUpper_[iCol];
    double dual_infeasibility = 0;
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the negation of
      // the dual value signed by nonbasicMove
      dual_infeasibility = -basis_.nonbasicMove_[iCol] * dual;
    }
    if (dual_infeasibility > options_->dual_feasibility_tolerance) {
      // Mustn't be an equality
      assert(lower<upper);
      num_dual_infeasibilities++;
      // Dual infeasible, so store the index and dual merit
      heap_index.push_back(iCol);
      heap_value.push_back(dual_infeasibility);
    }
  }
  // sifted list contains basic columns and, if primal feasible, a
  // selection of candidates for the sifted list is ready for sorting
  std::vector<double> new_col_cost;
  std::vector<double> new_col_lower;
  std::vector<double> new_col_upper;
  HighsSparseMatrix new_a_matrix;
  std::vector<HighsBasisStatus> new_col_status;
  if (!primal_feasible) {
    // Add a random collection of max_add_to_sifted_list nonbasic
    // columns to the sifted list, unless the sifted list is approaching
    // the number of columns in the LP
    HighsRandom random;
    const bool use_random =
        sifted_list.size() + 2 * max_add_to_sifted_list < lp_.num_col_;
    assert(use_random);
    HighsInt iCol = -1;
    // ToDo: Compile initial LP and basis in new data structures
    for (;;) {
      iCol = use_random ? random.integer(lp_.num_col_) : iCol++;
      if (in_sifted_list[iCol]) continue;
      // Add to sifted list
      num_add_to_sifted_list++;
      sifted_list.push_back(iCol);
      in_sifted_list[iCol] = true;
      new_col_cost.push_back(info_.workCost_[iCol]);
      new_col_lower.push_back(info_.workLower_[iCol]);
      new_col_upper.push_back(info_.workUpper_[iCol]);
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol];
           iEl < lp_.a_matrix_.start_[iCol + 1]; iEl++) {
        new_a_matrix.index_.push_back(lp_.a_matrix_.index_[iEl]);
        new_a_matrix.value_.push_back(lp_.a_matrix_.value_[iEl]);
      }
      new_a_matrix.start_.push_back(new_a_matrix.index_.size());
      assert(basis_.nonbasicFlag_[iCol]);
      HighsBasisStatus status = HighsBasisStatus::kBasic;
      if (basis_.nonbasicMove_[iCol] > 0) {
        status = HighsBasisStatus::kLower;
      } else if (basis_.nonbasicMove_[iCol] < 0) {
        status = HighsBasisStatus::kUpper;
      } else if (info_.workLower_[iCol] == info_.workUpper_[iCol]) {
        status = HighsBasisStatus::kLower;
      } else {
        status = HighsBasisStatus::kZero;
      }
      new_col_status.push_back(status);
      // Break out if the number of columns added to the sifted list
      // has reached max_add_to_sifted_list, or if the sifted list
      // contains all columns
      if (num_add_to_sifted_list == max_add_to_sifted_list ||
          sifted_list.size() >= lp_.num_col_)
        break;
    }
  } else {
    // Sort any candidates for the sifted list
    HighsInt heap_num_en = heap_index.size() - 1;
    assert(heap_num_en >= 0);
    HighsInt sifted_list_count = sifted_list.size();
    if (heap_num_en == 0) {
      // No dual infeasibilities, so should be optimal!
      assert(primal_feasible);
      assert(num_add_to_sifted_list == 0);
      return num_add_to_sifted_list;
    }
    // There are dual infeasibilities, so sort by increasing dual
    // merit - since decreasing sort isn't available.
    maxheapsort(&heap_value[0], &heap_index[0], heap_num_en);
    // Take from the end of the heap until number of columns added to
    // the sifted list has reached max_add_to_sifted_list, or all heap
    // entries have been added
    for (HighsInt iEl = heap_num_en; iEl > 0; iEl--) {
      HighsInt iCol = heap_index[iEl];
      num_add_to_sifted_list++;
      sifted_list.push_back(iCol);
      in_sifted_list[iCol] = true;
      new_col_cost.push_back(info_.workCost_[iCol]);
      new_col_lower.push_back(info_.workLower_[iCol]);
      new_col_upper.push_back(info_.workUpper_[iCol]);
      for (HighsInt iEl = lp_.a_matrix_.start_[iCol];
           iEl < lp_.a_matrix_.start_[iCol + 1]; iEl++) {
        new_a_matrix.index_.push_back(lp_.a_matrix_.index_[iEl]);
        new_a_matrix.value_.push_back(lp_.a_matrix_.value_[iEl]);
      }
      new_a_matrix.start_.push_back(new_a_matrix.index_.size());
      assert(basis_.nonbasicFlag_[iCol]);
      HighsBasisStatus status = HighsBasisStatus::kBasic;
      if (basis_.nonbasicMove_[iCol] > 0) {
        status = HighsBasisStatus::kLower;
      } else if (basis_.nonbasicMove_[iCol] < 0) {
        status = HighsBasisStatus::kUpper;
      } else if (info_.workLower_[iCol] == info_.workUpper_[iCol]) {
        status = HighsBasisStatus::kLower;
      } else {
        status = HighsBasisStatus::kZero;
      }
      new_col_status.push_back(status);
      if (num_add_to_sifted_list == max_add_to_sifted_list) break;
    }
  }
  sifted_lp.col_cost_ = new_col_cost;
  sifted_lp.col_lower_ = new_col_lower;
  sifted_lp.col_upper_ = new_col_upper;
  sifted_lp.a_matrix_ = new_a_matrix;
  sifted_basis.col_status = new_col_status;
  // Determine the row activity corresponding to columns not in the
  // sifted list so that the row bounds can be modified appropriately
  std::vector<double> unsifted_row_activity;
  unsifted_row_activity.assign(sifted_lp_num_row, 0);
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    if (in_sifted_list[iCol]) continue;
    double value = info_.workValue_[iCol];
    for (HighsInt iEl = lp_.a_matrix_.start_[iCol];
         iEl < lp_.a_matrix_.start_[iCol + 1]; iEl++) {
      HighsInt iRow = lp_.a_matrix_.index_[iEl];
      unsifted_row_activity[iRow] += value * lp_.a_matrix_.value_[iEl];
    }
  }
  if (first_sifted_lp) {
    // Complete the first sifted LP by adding its rows
    sifted_lp.row_lower_.resize(sifted_lp_num_row);
    sifted_lp.row_upper_.resize(sifted_lp_num_row);
    for (HighsInt iRow = 0; iRow < sifted_lp_num_row; iRow++) {
      HighsInt iVar = lp_.num_col_ + iRow;
      sifted_lp.row_lower_[iRow] =
          -info_.workUpper_[iVar] - unsifted_row_activity[iRow];
      sifted_lp.row_upper_[iRow] =
          -info_.workLower_[iVar] - unsifted_row_activity[iRow];
    }
    // Set the dimensions of the first sifted LP
    sifted_lp.num_col_ = sifted_lp.col_lower_.size();
    sifted_lp.num_row_ = sifted_lp.row_lower_.size();
    assert(sifted_lp.num_col_ == sifted_list.size());
    assert(sifted_basis.col_status.size() == sifted_lp.num_col_);
    assert(sifted_lp.num_row_ == sifted_lp_num_row);
    sifted_lp.setMatrixDimensions();
    // Move the first sifted LP to HEkk
    sifted_solver_object.ekk_instance_.moveLp(sifted_solver_object);
    // Local sifted LP is not updated - so clear it
    sifted_lp.clear();
  } else {
    // Update the row bounds of the sifted LP
    sifted_lp_num_col = sifted_list.size();
    const HighsInt sifted_lp_num_tot = sifted_lp_num_col + sifted_lp_num_row;
    sifted_solver_object.ekk_instance_.info_.workLower_.resize(
        sifted_lp_num_tot);
    sifted_solver_object.ekk_instance_.info_.workUpper_.resize(
        sifted_lp_num_tot);
    for (HighsInt iRow = 0; iRow < sifted_lp_num_row; iRow++) {
      HighsInt iVar = lp_.num_col_ + iRow;
      HighsInt sifted_iVar = sifted_lp_num_col + iRow;
      sifted_solver_object.ekk_instance_.info_.workLower_[sifted_iVar] =
          info_.workLower_[iVar] - unsifted_row_activity[iRow];
      sifted_solver_object.ekk_instance_.info_.workUpper_[sifted_iVar] =
          info_.workUpper_[iVar] - unsifted_row_activity[iRow];
    }
    // Add new columns to the LP in HEkk
    sifted_solver_object.ekk_instance_.addColsToLp(
        num_add_to_sifted_list, new_col_cost, new_col_lower, new_col_upper,
        new_a_matrix, new_col_status);
  }
  return num_add_to_sifted_list;
}

void HEkk::getSiftedLpInfo(HighsLpSolverObject& sifted_solver_object) {
  const HEkk& sifted_ekk_instance = sifted_solver_object.ekk_instance_;
  // Save data on infeasibilities
  info_.num_primal_infeasibilities =
      sifted_ekk_instance.info_.num_primal_infeasibilities;
  info_.num_dual_infeasibilities =
      sifted_ekk_instance.info_.num_dual_infeasibilities;

  // Copy over iteration count data
  iteration_count_ = sifted_ekk_instance.iteration_count_;
  info_.dual_phase1_iteration_count =
      sifted_ekk_instance.info_.dual_phase1_iteration_count;
  info_.dual_phase2_iteration_count =
      sifted_ekk_instance.info_.dual_phase2_iteration_count;
  info_.primal_phase1_iteration_count =
      sifted_ekk_instance.info_.primal_phase1_iteration_count;
  info_.primal_phase2_iteration_count =
      sifted_ekk_instance.info_.primal_phase2_iteration_count;
  info_.primal_bound_swap = sifted_ekk_instance.info_.primal_bound_swap;

  // Set analysis data so that reporting is correct
  analysis_.simplex_iteration_count = iteration_count_;
  analysis_.objective_value = sifted_ekk_instance.info_.primal_objective_value;
}

void HEkk::getSiftedBasisSolutionInvertWeightsInfo(
    HighsLpSolverObject& sifted_solver_object,
    const std::vector<HighsInt>& sifted_list) {
  const HEkk& sifted_ekk_instance = sifted_solver_object.ekk_instance_;
  const SimplexBasis& sifted_basis = sifted_ekk_instance.basis_;
  const HighsInt sifted_lp_num_col = sifted_list.size();
  const HighsInt sifted_lp_num_row = lp_.num_row_;

  // Get the nonbasic info for columns in the final sifted LP
  for (HighsInt iX = 0; iX < sifted_lp_num_col; iX++) {
    HighsInt iCol = sifted_list[iX];
    basis_.nonbasicFlag_[iCol] = sifted_basis.nonbasicFlag_[iX];
    basis_.nonbasicMove_[iCol] = sifted_basis.nonbasicMove_[iX];
  }
  const HighsInt var_offset = lp_.num_col_ - sifted_lp_num_col;
  // Get the nonbasic info for rows in the final sifted LP, and basic
  // info for variables in the final sifted LP
  for (HighsInt iRow = 0; iRow < sifted_lp_num_row; iRow++) {
    basis_.nonbasicFlag_[lp_.num_col_ + iRow] =
        sifted_basis.nonbasicFlag_[sifted_lp_num_col + iRow];
    basis_.nonbasicMove_[lp_.num_col_ + iRow] =
        sifted_basis.nonbasicMove_[sifted_lp_num_col + iRow];
    HighsInt iX = sifted_basis.basicIndex_[iRow];
    if (iX < sifted_lp_num_col) {
      basis_.basicIndex_[iRow] = sifted_list[iX];
    } else {
      basis_.basicIndex_[iRow] = var_offset + iX;
    }
  }
  // Clear statuses for a new basis
  this->updateStatus(LpAction::kNewBasis);
  this->status_.has_basis = true;
  // Copy any dual edge weights
  this->getDualEdgeWeights(sifted_solver_object.ekk_instance_);
  // Copy the INVERT
  this->getInvert(sifted_solver_object.ekk_instance_);

  getSiftedLpInfo(sifted_solver_object);
}

bool HEkk::debugOkSiftedList(const std::vector<HighsInt>& sifted_list,
                             const std::vector<bool>& in_sifted_list) {
  std::vector<bool> local_in_sifted_list = in_sifted_list;
  for (HighsInt iX = 0; iX < sifted_list.size(); iX++) {
    if (!local_in_sifted_list[sifted_list[iX]]) {
      printf("local_in_sifted_list[sifted_list[%d]] is false\n", (int)iX);
      return false;
    }
    local_in_sifted_list[sifted_list[iX]] = false;
  }
  for (HighsInt iCol = 0; iCol < sifted_list.size(); iCol++) {
    if (local_in_sifted_list[iCol]) {
      printf("local_in_sifted_list[%d] is true\n", (int)iCol);
      return false;
    }
  }
  return true;
}
