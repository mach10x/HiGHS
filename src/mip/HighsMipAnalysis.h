/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsMipAnalysis.h
 * @brief Analyse MIP iterations, both for run-time control and data
 * gathering
 */
#ifndef MIP_HIGHSMIPANALYSIS_H_
#define MIP_HIGHSMIPANALYSIS_H_

#include "lp_data/HighsAnalysis.h"
#include "lp_data/HighsLp.h"
#include "util/HighsTimer.h"

class HighsMipAnalysis {
 public:
  HighsMipAnalysis() : timer_(nullptr), analyse_mip_time(false) {}

  HighsTimer* timer_;
  void setup(const HighsLp& lp, const HighsOptions& options);

  void setupMipTime(const HighsOptions& options);
  void mipTimerStart(const HighsInt mip_clock
                     //		     , const HighsInt thread_id = 0
  );
  void mipTimerStop(const HighsInt mip_clock
                    //		    , const HighsInt thread_id = 0
  );
  void reportMipTimer();

  HighsTimerClock mip_clocks;
  bool analyse_mip_time;
};

#endif /* MIP_HIGHSMIPANALYSIS_H_ */