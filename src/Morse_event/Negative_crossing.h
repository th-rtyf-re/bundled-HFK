/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *  Bundled HFK - a knot Floer homology calculator                           *
 *                                                                           *
 *  Copyright (C) 2021-2022  Isaac Ren                                       *
 *  For further details, contact Isaac Ren (gopi3.1415@gmail.com).           *
 *                                                                           *
 *  This program is free software: you can redistribute it and/or modify     *
 *  it under the terms of the GNU General Public License as published by     *
 *  the Free Software Foundation, either version 3 of the License, or        *
 *  (at your option) any later version.                                      *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.   *
 *                                                                           *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef NEGATIVE_CROSSING_H_
#define NEGATIVE_CROSSING_H_

#include <string>
#include <vector>

#ifdef BUNDLED_HFK_VERBOSE_
#include <iostream>
#endif  // BUNDLED_HFK_VERBOSE_

#include <boost/any.hpp>

#include "Positive_crossing.h"
#include "Math_tools/Reverse_D_module.h"

/* Morse event for negative crossings.
 * 
 * The DA-bimodule of a negative crossing is the dual of the corresponding
 * positive crossing. Therefore this class just calls the positive crossing
 * class on a reverse view of the old D-module, and then reverses back the new
 * D-module.
 */
template< class D_module, class Morse_event_options >
class Negative_crossing {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Algebra = typename D_module::Bordered_algebra;
  using Coef_bundle = typename D_module::Coef_bundle;
  using Weights = typename D_module::Weights;  // currently, pair of int
  
  using Positive_crossing = Positive_crossing< Reverse_D_module< D_module >, Morse_event_options >;
  
  Negative_crossing(const std::vector< typename Morse_event_options::Parameter_type >& args) :
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
  
#ifdef BUNDLED_HFK_VERBOSE_
  friend std::ostream& operator<<(
    std::ostream& os,
    const Negative_crossing& morse_event
  ) {
    os << "- at "
      << morse_event.positive_crossing_.position_;
    return os;
  }
#endif  // BUNDLED_HFK_VERBOSE_
  
 private:
  Positive_crossing positive_crossing_;
};

#endif  // NEGATIVE_CROSSING_H_