/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *  Bundled HFK - a knot Floer homology calculator                           *
 *                                                                           *
 *  Copyright (C) 2021  Isaac Ren                                            *
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

#ifndef EVENT_GLOBAL_MINIMUM_H_
#define EVENT_GLOBAL_MINIMUM_H_

#include <string>

#include "../DA_bimodule.h"
#include "../Knot_slice.h"

/* Knot slice for the global minimum
 * 
 * This is where we fix the orientation of the knot.
 */
class Knot_slice_global_minimum : public Knot_slice {
 public:
  Knot_slice_global_minimum(int, int) :
    Knot_slice(2, 0, 0)
  { }
  
  ~Knot_slice_global_minimum() = default;
  
  Morse_event event_id() const override {
    return Morse_event::global_minimum;
  }
  
  std::vector< int > lower_matchings(std::vector< int >) const override {
    return std::vector< int >(0);
  }
  
  std::vector< bool > upper_orientations(std::vector< bool >, const std::vector< int >&) const override {
    return {false, true};  // we choose the trigonometric orientation
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(const std::pair< int, int >& margins) const override {
    return std::string(margins.first - 1, '0') + "u" + std::string(margins.second - 1, '0');
  }
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const override {
    return {margins.first + 1, margins.second + 1};
  }
  
 private:
  using Knot_slice::upper_n_strands_;
  using Knot_slice::lower_n_strands_;
  using Knot_slice::position_;
};

/* Local DA-bimodule for a local minimum.
 * 
 * A slightly formal DA-bimodule that just forgets all coefficients.
 */
template< class D_module >
class DA_bimodule_global_minimum : public DA_bimodule< D_module > {
 public:
  using DA_bimodule = DA_bimodule< D_module >;
  using Idem = typename DA_bimodule::Idem;
  using Gen_type = typename DA_bimodule::Gen_type;
  using Algebra = typename DA_bimodule::Algebra;
  using Monomial = typename DA_bimodule::Monomial;
  
  DA_bimodule_global_minimum(const std::unique_ptr< Knot_slice >&, Algebra, Algebra) { }
  
  ~DA_bimodule_global_minimum() = default;
  
  Monomial alexander_maslov_weights(const Gen_type generator_type) const override {
    return {0, 0};
  }
  
  std::string get_label(const Gen_type generator_type) const override {
    return "{}";
  }
  
  /* Requirement: the old D-module has 2 strands and one agglomerated generator,
   * with lower idempotent 010.
   */
  D_module box_tensor_product(const D_module& old_d_module) const override {
    D_module new_d_module;
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      new_d_module.add_gen_bundle(Idem("0"), 0, gen_handle, "m");
    }
    return new_d_module;
  }
};

#endif  // EVENT_GLOBAL_MINIMUM_H_