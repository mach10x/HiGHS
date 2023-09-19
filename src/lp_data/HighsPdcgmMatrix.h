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
/**@file lp_data/HighsPdcgmMatrix.h
 * @brief
 */
#ifndef LP_DATA_HIGHS_PDCGMMATRIX_H_
#define LP_DATA_HIGHS_PDCGMMATRIX_H_

#include "util/HighsInt.h"
#include <vector>

class HighsPdcgm_SMatrix_CW {
 public:
  HighsPdcgm_SMatrix_CW() { clear(); }
  HighsInt m_m{};
  HighsInt m_n{};
  HighsInt m_nz{};

  HighsInt m_max_m{};
  HighsInt m_max_n{};
  HighsInt m_max_nz{};

  std::vector<double> m_obj{};

  std::vector<double> m_primal_lb{};
  std::vector<double> m_primal_ub{};
  std::vector<double> m_primal_lb_fxd{};
  std::vector<double> m_primal_ub_fxd{};

  std::vector<double> m_dual_lb{};
  std::vector<double> m_dual_ub{};
  std::vector<double> m_dual_lb_fxd{};
  std::vector<double> m_dual_ub_fxd{};

  std::vector<double> m_coeff{};
  std::vector<HighsInt> m_rwnmbs{};
  std::vector<HighsInt> m_clpnts{};

  std::vector<double> m_rhs{};
  std::vector<char> m_constr_type{};

  void clear();
};

#endif
