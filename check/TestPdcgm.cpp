#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;
//const double inf = kHighsInf;



// No commas in test case name.
TEST_CASE("pdcgm", "[highs_pdcgm]") {

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  

}
