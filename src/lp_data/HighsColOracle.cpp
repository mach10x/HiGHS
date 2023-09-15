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
/**@file lp_data/HighsColOracle.cpp
 * @brief
 */
#include "HighsColOracle.h"

#include <cassert>

void HighsColOracleDataOut::clear() {}

void HighsColOracleDataIn::clear() { this->user_interrupt = false; }

void HighsColOracle::clear() {
  this->col_oracle = nullptr;
  this->col_oracle_data = nullptr;
  this->active.assign(kNumHighsColOracleType, false);
  this->data_out.clear();
  this->data_in.clear();
}

bool HighsColOracle::colOracleActive(const int col_oracle_type) {
  // Check that col_oracle function has been defined
  if (!this->col_oracle) return false;
  // Check that col_oracle_type is within range
  const bool col_oracle_type_ok = col_oracle_type >= kHighsColOracleMin &&
                                  col_oracle_type <= kHighsColOracleMax;
  assert(col_oracle_type_ok);
  if (!col_oracle_type_ok) return false;
  // Don't call col_oracle if it is not active
  assert(this->active.size() > 0);
  if (!this->active[col_oracle_type]) return false;
  return true;
}

bool HighsColOracle::colOracleAction(const int col_oracle_type,
                                     std::string message) {
  if (!colOracleActive(col_oracle_type)) return false;
  this->col_oracle(col_oracle_type, message.c_str(), &this->data_out,
                   &this->data_in, this->col_oracle_data);
  // Assess any action
  bool action = this->data_in.user_interrupt;

  return action;
}
