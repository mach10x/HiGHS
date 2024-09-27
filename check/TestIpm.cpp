#include <cmath>

#include "Highs.h"
#include "catch.hpp"

// I use dev_run to switch on/off printing and logging used for
// development of the unit test
const bool dev_run = false;
const double inf = kHighsInf;

TEST_CASE("test-analytic-centre", "[highs_ipm]") {
  //  std::string model = "greenbea.mps";
  //  std::string model = "adlittle.mps";
  std::string model = "afiro.mps";
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/" + model;
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(filename);
  HighsLp lp = highs.getLp();
  lp.col_cost_.assign(lp.num_col_, 0);
  highs.passModel(lp);
  highs.setOptionValue("run_centring", true);
  highs.setOptionValue("ipm_optimality_tolerance", 1e-2);
  HighsStatus run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);
}

TEST_CASE("test-analytic-centre-infeasible", "[highs_ipm]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_.assign(lp.num_col_, 0);
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, inf);
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {-1};
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};
  highs.passModel(lp);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("run_centring", true);
  highs.setOptionValue("ipm_optimality_tolerance", 1e-2);
  HighsStatus run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
}

TEST_CASE("test-analytic-centre-box", "[highs_ipm]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsInt dim = 4;
  HighsLp lp;
  lp.num_col_ = dim;
  lp.col_cost_.assign(dim, 0);
  lp.col_lower_.assign(dim, -1);
  lp.col_upper_.assign(dim, 1);
  highs.passModel(lp);

  std::vector<HighsInt> index = {0, 1};
  std::vector<double> value = {1, 1};

  const double root2 = std::sqrt(2.0);
  highs.addRow(-root2, root2, 2, index.data(), value.data());
  value[1] = -1;
  highs.addRow(-root2, root2, 2, index.data(), value.data());
  highs.setOptionValue("run_centring", true);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("ipm_optimality_tolerance", 1e-2);
  HighsStatus run_status = highs.run();
  const HighsSolution& solution = highs.getSolution();
  double solution_norm = 0;
  for (HighsInt ix = 0; ix < dim; ix++) {
    if (dev_run)
      printf("Analytic centre solution %d is %g\n", int(ix),
             solution.col_value[ix]);
    solution_norm += std::fabs(solution.col_value[ix]);
  }
  REQUIRE(solution_norm < 1e-6);
  if (dev_run) printf("Analytic centre solution norm is %g\n", solution_norm);
}


TEST_CASE("ipm-primal-dual-infeasible-instance", "[highs_ipm]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsInt dim = 2;
  HighsLp lp;

  lp.num_col_ = dim;
  lp.num_row_ = dim;

  lp.col_cost_ = {2.0, -1.0};
  lp.sense_ = ObjSense::kMaximize;

  lp.col_lower_ = {dim, 0.0};
  lp.col_upper_= {dim, kHighsInf};

  lp.row_lower_= {-kHighsInf, 2.0};
  lp.row_upper_= {1.0, kHighsInf};

  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, -1, -1};

  HighsStatus status = highs.passModel(lp);
  REQUIRE(status == HighsStatus::kOk);

  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("solver", "ipm");

  const HighsOptions& options = highs.getOptions();

  HighsStatus run_status = highs.run();
  const HighsSolution& solution = highs.getSolution();
  const HighsInfo& info = highs.getInfo();

  std::cout << solution.value_valid << std::endl;
  std::cout << solution.dual_valid << std::endl;

    
  std::cout << info.primal_solution_status << std::endl;
  std::cout << info.dual_solution_status << std::endl;
  std::cout << solution.col_value[0] << " " << solution.col_value[1] << std::endl;
  std::cout << solution.col_dual[0] << " " << solution.col_dual[1] << std::endl;
  REQUIRE(info.primal_solution_status == kSolutionStatusInfeasible);
}
