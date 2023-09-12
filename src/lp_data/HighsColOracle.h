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
/**@file lp_data/HighsColOracle.h
 * @brief
 */
#ifndef LP_DATA_HIGHSCOLORACLE_H_
#define LP_DATA_HIGHSCOLORACLE_H_

#include "lp_data/HStruct.h"

struct HighsColOracleDataOut {
  void clear();
};

struct HighsColOracleDataIn {
  bool user_interrupt;
  void clear();
};

struct HighsColOracle {
  void (*col_oracle)(const int, const char*, const HighsColOracleDataOut*,
                        HighsColOracleDataIn*, void*) = nullptr;
  void* col_oracle_data = nullptr;
  std::vector<bool> active;
  HighsColOracleDataOut data_out;
  HighsColOracleDataIn data_in;
  bool colOracleActive(const int col_oracle_type);
  bool colOracleAction(const int col_oracle_type, std::string message = "");
  void clear();
};
#endif /* LP_DATA_HIGHSCOLORACLE_H_ */
