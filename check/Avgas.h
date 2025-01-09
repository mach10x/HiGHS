/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/Avgas.h
 * @brief Utilities for tests with AVGAS
 */
#ifndef SIMPLEX_AVGAS_H_
#define SIMPLEX_AVGAS_H_

#include <vector>

#include "util/HighsInt.h"

/**
 * @brief Utilities for tests with AVGAS
 */
class Avgas {
 public:
  void addRow(HighsInt row, HighsInt& num_row, HighsInt& num_row_nz,
	      std::vector<double>& rowLower, std::vector<double>& rowUpper,
	      std::vector<HighsInt>& ARstart, std::vector<HighsInt>& ARindex,
	      std::vector<double>& ARvalue);

  void getRow(HighsInt row, 
	      std::vector<double>& lower, std::vector<double>& upper,
	      HighsInt& num_nz,
	      std::vector<HighsInt>& index,
	      std::vector<double>& value);

  void addCol(HighsInt col, HighsInt& num_col, HighsInt& num_col_nz,
	      std::vector<double>& colCost, std::vector<double>& colLower,
	      std::vector<double>& colUpper, std::vector<HighsInt>& Astart,
	      std::vector<HighsInt>& Aindex, std::vector<double>& Avalue);

  void getCol(HighsInt col, 
	      std::vector<double>& colCost, std::vector<double>& colLower,
	      std::vector<double>& colUpper,
	      HighsInt& num_nz, 
	      std::vector<HighsInt>& index, std::vector<double>& value);
};
#endif /* SIMPLEX_AVGAS_H_ */
