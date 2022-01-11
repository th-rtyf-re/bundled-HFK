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

#ifndef GLOBAL_MINIMUM_H_
#define GLOBAL_MINIMUM_H_

#include <string>
#include <utility>  // pair
#include <vector>

#ifdef BUNDLED_HFK_VERBOSE_
#include <iostream>
#endif  // BUNDLED_HFK_VERBOSE_

/* Morse event for a global minimum.
 */
template< class D_module, class Morse_event_options >
class Global_minimum {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Algebra = typename D_module::Bordered_algebra;
  using Weights = typename D_module::Weights;
  
  Global_minimum(const std::vector< typename Morse_event_options::Parameter_type >&) { }
  
  /* Topological methods */
  
  std::vector< int > lower_matchings(std::vector< int >) const {
    return {0};
  }
  
  std::vector< bool > upper_orientations(
    std::vector< bool >,
    const std::vector< int >&
  ) const {
    return {false, true};  // we choose the trigonometric orientation
  }
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const {
    return {margins.first + 1, margins.second + 1};
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(
    const std::pair< int, int >& margins,
    const std::pair< int, int >&
  ) const {
    return std::string(
      margins.first - 1, '0') + "u" + std::string(margins.second - 1, '0'
    );
  }
  
  /* Algebraic methods */
  
  std::vector< Weights > get_weights(const Algebra&, const Algebra&) const {
    return {{0, 0}};
  }
  
  std::vector< std::string > get_labels(const Algebra&, const Algebra&) const {
    return {"{}"};
  }
  
  D_module tensor_generators(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra&,
    const Algebra&
  ) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      new_d_module.add_gen_bundle(Idem("0"), 0, gen_handle);
    }
    return new_d_module;
  }
  
  D_module& tensor_coefficients(
    D_module& new_d_module,
    const D_module&,
    const Algebra&,
    const Algebra&
  ) const {
    return new_d_module;
  }
  
#ifdef BUNDLED_HFK_VERBOSE_
  friend std::ostream& operator<<(
    std::ostream& os,
    const Global_minimum& morse_event
  ) {
    os << "global min at 0";
    return os;
  }
#endif  // BUNDLED_HFK_VERBOSE_
};

#endif  // GLOBAL_MINIMUM_H_