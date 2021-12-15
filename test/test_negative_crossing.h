#ifndef TEST_NEGATIVE_CROSSING_H_
#define TEST_NEGATIVE_CROSSING_H_

#include <string>
#include <vector>

#ifdef VERBOSE
#include <iostream>
#endif  // VERBOSE

#include <boost/any.hpp>

#include "test_positive_crossing.h"
#include "test_reverse_d_module.h"

template< class D_module >
class Negative_crossing {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Algebra = typename D_module::Bordered_algebra;
  using Coef_bundle = typename D_module::Coef_bundle;
  using Weights = typename D_module::Weights;  // currently, pair of int
  
  using Positive_crossing = Positive_crossing< Reverse_D_module< D_module > >;
  
  Negative_crossing(const std::vector< boost::any >& args) :
    positive_crossing_(args)
  { }
  
    /* Topological methods */
  
  std::vector< int > lower_matchings(std::vector< int > matchings) const {
    return positive_crossing_.lower_matchings(matchings);
  }
  
  std::vector< bool > upper_orientations(
    std::vector< bool > orientations,
    const std::vector< int >& upper_matchings
  ) const {
    return positive_crossing_.upper_orientations(
      orientations,
      upper_matchings
    );
  }
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const {
    return positive_crossing_.update_margins(margins);
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(
    const std::pair< int, int >& margins,
    const std::pair< int, int >& n_strands,
    std::string symbol = "-"
  ) const {
    return positive_crossing_.to_string(margins, n_strands, symbol);
  }
  
  /* Algebraic methods */
  
  std::vector< Weights > get_weights(
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    std::vector< Weights > weights = positive_crossing_.get_weights(
      upper_algebra,
      lower_algebra
    );
    for (auto& value_pair : weights) {
      value_pair.first = -value_pair.first;
      value_pair.second = -value_pair.second;
    }
    return weights;
  }
  
  std::vector< std::string > get_labels(
    const Algebra& upper_algebra,
    const Algebra& lower_algebra,
    std::string symbols = "nesw"
  ) const {
    return positive_crossing_.get_labels(
      upper_algebra,
      lower_algebra,
      symbols
    );
  }
  
  D_module tensor_generators(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    auto reverse_new_d_module = reverse_view(new_d_module);
    const auto reverse_old_d_module = reverse_view(old_d_module);
    positive_crossing_.tensor_generators(
      reverse_new_d_module,
      reverse_old_d_module,
      upper_algebra,
      lower_algebra
    );
    return new_d_module;
  }
  
  D_module& tensor_coefficients(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    auto reverse_new_d_module = reverse_view(new_d_module);
    const auto reverse_old_d_module = reverse_view(old_d_module);
    positive_crossing_.tensor_coefficients(
      reverse_new_d_module,
      reverse_old_d_module,
      upper_algebra,
      lower_algebra
    );
    return new_d_module;
  }
  
#ifdef VERBOSE
  friend std::ostream& operator<<(
    std::ostream& os,
    const Negative_crossing& morse_event
  ) {
    os << "negative crossing at position "
      << morse_event.positive_crossing_.position_;
    return os;
  }
#endif  // VERBOSE
  
 private:
  Positive_crossing positive_crossing_;
};

#endif  // TEST_NEGATIVE_CROSSING_H_