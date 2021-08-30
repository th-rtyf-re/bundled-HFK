/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *  Bundled HFK - a knot Floer homology calculator                           *
 *                                                                           *
 *  Copyright (C) 2021  Isaac Ren                                            *
 *  For further details contact Isaac Ren (gopi3.1415@gmail.com)             *
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

#ifndef EVENT_LOCAL_MAXIMUM_H_
#define EVENT_LOCAL_MAXIMUM_H_

#include <string>
#include <tuple>  // tuple, tie
#include <vector>

#include "../DA_bimodule.h"
#include "../Knot_slice.h"

/* Knot slice for a local maximum.
 */
class Knot_slice_local_maximum : public Knot_slice {
 public:
  Knot_slice_local_maximum(int upper_n_strands, int position) :
    Knot_slice(upper_n_strands, upper_n_strands + 2, position)
  { }
  
  ~Knot_slice_local_maximum() = default;
  
  Event event_id() const override {
    return Event::local_maximum;
  }
  
  std::vector< int > lower_matchings(std::vector< int > matchings) const override {
    for (int& match : matchings) {
      if (match >= position_) {
        match += 2;
      }
    }
    matchings.insert(matchings.begin() + position_, {position_ + 1, position_});
    return matchings;
  }
  
  std::vector< bool > upper_orientations(std::vector< bool > orientations, const std::vector< int >&) const override {
    orientations.erase(orientations.begin() + position_,
    orientations.begin() + position_ + 2);
    return orientations;
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(const std::pair< int, int >& margins) const override {
    return std::string(margins.first, '0')
      + std::string(position_, 'r')
      + "a"
      + std::string(lower_n_strands_ - position_ - 2, 'l')
      + std::string(margins.second, '0');
  }
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const override {
    return {margins.first - 1, margins.second - 1};
  }
  
 private:
  using Knot_slice::upper_n_strands_;
  using Knot_slice::lower_n_strands_;
  using Knot_slice::position_;
};

/* DA-bimodule for a local maximum.
 * 
 * Adapted from ComputeHFKv2/Max.cpp
 */
template< class D_module >
class DA_bimodule_local_maximum : public DA_bimodule< D_module > {
 public:
  using DA_bimodule = DA_bimodule< D_module >;
  using Idem = typename DA_bimodule::Idem;
  using Gen_type = typename DA_bimodule::Gen_type;
  using Algebra = typename DA_bimodule::Algebra;
  using Alg_el = typename DA_bimodule::Alg_el;
  using Monomial = typename DA_bimodule::Monomial;
  
  DA_bimodule_local_maximum(const std::unique_ptr< Knot_slice >& slice, Algebra upper_algebra, Algebra lower_algebra) :
    position_(slice->position()),
    DA_bimodule(upper_algebra, lower_algebra)
  { }
  
  ~DA_bimodule_local_maximum() = default;
  
  Monomial alexander_maslov_weights(const Gen_type generator_type) const override {
    return {0, 0};
  }
  
  std::string get_label(const Gen_type generator_type) const override {
    if (lower_algebra_.orientations[position_]) {  // clock
      return std::string(1, "XYZ"[generator_type])
        + "_{" + std::to_string(position_) + "}";
    }
    else {  // trig
      return std::string(1, "xyz"[generator_type])
        + "_{" + std::to_string(position_) + "}";
    }
  }
  
  D_module box_tensor_product(const D_module& old_d_module) const override {
    D_module new_d_module;
    delta_0_and_1_(old_d_module, new_d_module);
    delta_2_(old_d_module, new_d_module);
    return new_d_module;
  }
  
 private:
  enum : Gen_type {
    X,
    Y,
    Z
  };
  
  void delta_0_and_1_(const D_module& old_d_module, D_module& new_d_module) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      const Idem& old_idem = old_d_module.idem(gen_handle);
      
      if (old_idem[position_]) {  // X and Y
        Idem x_idem = extend_(old_idem, X);
        new_d_module.add_gen_bundle(x_idem, X, gen_handle, "X");
        Idem y_idem = extend_(old_idem, Y);
        new_d_module.add_gen_bundle(y_idem, Y, gen_handle, "Y");
        new_d_module.add_coef_bundle(Alg_el(x_idem, y_idem), X, Y, old_idem);
        new_d_module.add_coef_bundle(Alg_el(y_idem, x_idem), Y, X, old_idem);
      }
      
      else {  // Z
        Idem z_idem = old_idem;
        z_idem.insert(position_, {false, true});
        new_d_module.add_gen_bundle(z_idem, Z, gen_handle, "Z");
      }
    }
  }
  
  /* See [OzsvathSzabo2018, Lemma 8.1] */
  void delta_2_(const D_module& old_d_module, D_module& new_d_module) const {
    for (const auto& coef : old_d_module.coef_bundles()) {
      const Alg_el& alg_el = coef.algebra_element();
      const Idem& back_idem = alg_el.source_idem();
      const Idem& front_idem = alg_el.target_idem();
      std::vector< int > new_U_weights = alg_el.U_weights();
      
      int a1, a2;
      std::tie(a1, a2) = get_local_weights_(coef.algebra_element());
      new_U_weights.insert(new_U_weights.begin() + position_, {0, 0});
      
      /* Lambda expression. Captures back_idem, front_idem, new_U_weights,
       * coef, and new_d_module.
       * 
       * Note: I wanted to make this static, but it doesn't work with coef. 
       */
      auto compose = [&](const Gen_type back_marking, const Gen_type front_marking) {
        const Idem& new_source_idem = extend_(back_idem, back_marking);
        const Idem& new_target_idem = extend_(front_idem, front_marking);
        if (new_source_idem.too_far_from(new_target_idem)) { return; }
        const Alg_el new_alg_el(new_source_idem, new_target_idem, new_U_weights);
        if (new_alg_el.is_null()) { return; }
        new_d_module.add_coef_bundle(new_alg_el, back_marking, front_marking, coef);
      };
      
      Gen_type back_marking;
      Gen_type front_marking;
      
      if (a1 == 1 and a2 == 1) {  // R_2R_1
        compose(Y, X);
      }
      else if (a1 == -1 and a2 == -1) {  // L_1L_2
        compose(X, Y);
      }
      else if (a1 == 1 and a2 == 0) {  // R_1
        compose(Z, X);
      }
      else if (a1 == -1 and a2 == 0) {  // L_1
        compose(X, Z);
      }
      else if (a1 == 0 and a2 == 1) {  // R_2
        compose(Y, Z);
      }
      else if (a1 == 0 and a2 == -1) {  // L_2
        compose(Z, Y);
      }
      else if (a1 == 0 and a2 == 0) {
        if (back_idem[position_]) {  // X, Y
          compose(X, X);
          compose(Y, Y);
        }
        else {  // Z
          compose(Z, Z);
        }
      }
    }
  }  // delta_2_
  
  /* Get LR weights at position_ - 1 and position_. This corresponds to the
   * local LR weights of the upper algebra.
   */
  std::tuple< int, int > get_local_weights_(const Alg_el& alg_el) const {
    const Idem& source_idem = alg_el.source_idem();
    const Idem& target_idem = alg_el.target_idem();
    int a1 = 0;  // will only take values in [-1, 1]
    int a2 = 0;  // will only take values in [-1, 1]
    int i = 0;
    for ( ; i < position_; ++i) {
      a1 = a1 + source_idem[i] - target_idem[i];
    }
    a2 = a1 + source_idem[i] - target_idem[i];  // i = position_
    return std::make_tuple(a1, a2);
  }  // get_local_weights_
  
  Idem extend_(Idem idem, const Gen_type marking) const {
    if (marking == X) {  // 1 -> 110
      idem.insert(position_ + 1, {1, 0});
    }
    else if (marking == Y or marking == Z) {  // x -> 01x
      idem.insert(position_, {0, 1});
    }
    return idem;
  }
  
  int position_;
  
  using DA_bimodule::upper_algebra_;
  using DA_bimodule::lower_algebra_;
};

#endif  // EVENT_LOCAL_MAXIMUM_H_