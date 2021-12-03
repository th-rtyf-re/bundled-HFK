#ifndef TEST_FOREST_OPTIONS_H_
#define TEST_FOREST_OPTIONS_H_

#include <vector>
#include <utility>  // pair

#include "Bordered_algebra/Bordered_algebra.h"
#include "Bordered_algebra/Idempotent.h"

struct Forest_options_default_short {
  using Idem = Idempotent_short;
  using Bordered_algebra = Bordered_algebra< Idem >;
  using Alg_el = typename Bordered_algebra::Element;
  using Gen_type = unsigned char;  // no need to pass by reference excessively
  using Weights = std::pair< int, int >;
};

struct Forest_options_default_long {
  using Idem = Idempotent_long< std::vector< bool > >;
  using Bordered_algebra = Bordered_algebra< Idem >;
  using Alg_el = typename Bordered_algebra::Element;
  using Gen_type = unsigned char;  // no need to pass by reference excessively
  using Weights = std::pair< int, int >;
};

#endif  // TEST_FOREST_OPTIONS_H_