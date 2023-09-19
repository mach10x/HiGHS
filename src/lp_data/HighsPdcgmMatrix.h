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

class HighsPdcgmMatrix {
 public:
  HighsPdcgmMatrix() { clear(); }
  void clear();
};

#endif
