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

#ifndef EVENT_CROSSING_H_
#define EVENT_CROSSING_H_

#include <cstdlib>  // abs
#include <string>
#include <tuple>
#include <utility>  // pair, swap
#include <vector>

#include "../DA_bimodule.h"
#include "../Knot_slice.h"
#include "Event_crossing_helper.h"

/* Knot slice and DA-bimodule for a crossing
 * 
 * Both of these classes are templated by a boolean, for the sign of the
 * crossing: false = negative, true = positive.
 */

template< bool sign >
class Knot_slice_crossing : public Knot_slice {
 public:
  Knot_slice_crossing(int upper_n_strands, int position) :
    Knot_slice(upper_n_strands, upper_n_strands, position)
  { }
  
  ~Knot_slice_crossing() = default;
  
  Morse_event event_id() const override {
    return sign ? Morse_event::positive_crossing : Morse_event::negative_crossing;
  }
  
  std::vector< int > lower_matchings(std::vector< int > matchings) const override {
    std::swap(matchings[position_], matchings[position_ + 1]);
    std::swap(matchings[matchings[position_]], matchings[matchings[position_ + 1]]);
    return matchings;
  }
  
  std::vector< bool > upper_orientations(std::vector< bool > orientations, const std::vector< int >&) const override {
    std::swap(orientations[position_], orientations[position_ + 1]);
    return orientations;
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(const std::pair< int, int >& margins) const override {
    return std::string(margins.first, '0')
      + std::string(position_, '.')
      + (sign ? "+" : "-")
      + std::string(lower_n_strands_ - position_ - 2, '.')
      + std::string(margins.second, '0');
  }
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const override {
    return margins;
  }
  
 private:
  using Knot_slice::upper_n_strands_;
  using Knot_slice::lower_n_strands_;
  using Knot_slice::position_;
};

/* DA bimodule for a crossing.
 * 
 * Adapted from ComputeHFK2/Crossing.cpp.
 */
template< bool sign, class D_module >
class DA_bimodule_crossing : public DA_bimodule< D_module > {
 public:
  using DA_bimodule = DA_bimodule< D_module >;
  using Idem = typename DA_bimodule::Idem;
  using Gen_type = typename DA_bimodule::Gen_type;
  using Algebra = typename DA_bimodule::Algebra;
  using Alg_el = typename DA_bimodule::Alg_el;
  using Monomial = typename DA_bimodule::Monomial;
  
  DA_bimodule_crossing(const std::unique_ptr< Knot_slice >& slice, Algebra upper_algebra, Algebra lower_algebra) :
    position_(slice->position()),
    DA_bimodule(upper_algebra, lower_algebra)
  { }
  
  ~DA_bimodule_crossing() = default;
  
  Monomial alexander_maslov_weights(const Gen_type generator_type) const override {
    bool left_ori = upper_algebra_.orientations[position_];
    bool right_ori = upper_algebra_.orientations[position_ + 1];
    int a = 0;  // double Alexander
    int m = 0;  // Maslov
    
    /* Case study based on the direction in which both strands are pointing */
    if (left_ori) {
      if (right_ori) {  // positive, pointing north
        if (generator_type == north) {
          a = -1;
          m = -1;
        }
        else if (generator_type == south) {
          a = 1;
        }
      }
      else {  // negative, pointing west
        if (generator_type == east) {
          a = -1;
        }
        else if (generator_type == west) {
          a = 1;
          m = 1;
        }
      }
    }
    else {
      if (right_ori) {  // negative, pointing east
        if (generator_type == east) {
          a = 1;
          m = 1;
        }
        else if (generator_type == west) {
          a = -1;
        }
      }
      else {  // positive, pointing south
        if (generator_type == north) {
          a = 1;
        }
        else if (generator_type == south) {
          a = -1;
          m = -1;
        }
      }
    }
    
    if (!sign) {  // flip signs if negative
      a = -a;
      m = -m;
    }
    return {a, m};
  }
  
  std::string get_label(const Gen_type generator_type) const override {
    if (sign) {  // positive crossing
      return std::string(1, "NESW"[generator_type])
        + "_{" + std::to_string(position_) + "}";
    }
    else {  // negative crossing
      return std::string(1, "nesw"[generator_type])
        + "_{" + std::to_string(position_) + "}";
    }
  }
  
  /* If the crossing is negative, we pass to the dual DA-bimodule, which we
   * box tensor product with the dual upper D-module. This process is linear in
   * the number of differential arcs, but priceless in number of eliminated
   * headaches.
   * 
   * Note: we need to const_cast the old D-module to do this. This is ok
   * because the old D-module is not initialized as const.
   */
  D_module box_tensor_product(const D_module& old_d_module) const override {
    D_module new_d_module;
    if (!sign) {  // dualize if negative crossing
      const_cast< D_module& >(old_d_module).dualize();
    }
    delta_0_and_1_(old_d_module, new_d_module);
    delta_2_(old_d_module, new_d_module);
    delta_3_(old_d_module, new_d_module);
    if (!sign) {  // dualize if negative crossing
      new_d_module.dualize();
    }
    return new_d_module;
  }  // box_tensor_product
  
 private:
  enum : Gen_type {
    north = 0,
    east,   // 1
    south,  // 2
    west,   // 3
    null_NEW = 2,
    null_S = 0,
    null = 0
  };
  
  /* \delta functions */
  
  /* \delta_0 and \delta_1.
   * 
   * \delta_0 refers to the tensor product part of the box tensor product, that
   * is, the addition of nodes and edges.
   */
  void delta_0_and_1_(const D_module& old_d_module, D_module& new_d_module) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      const Idem& old_idem = old_d_module.idem(gen_handle);
      
      if (old_idem[position_ + 1]) {  // north marking
        new_d_module.add_gen_bundle(old_idem, north, gen_handle, "N", alexander_maslov_weights(north));  // temporary?
      }
      else {  // east, south, west markings
        new_d_module.add_gen_bundle(old_idem, south, gen_handle, "S", alexander_maslov_weights(south));  // temporary?
        
        for (int we : {0, 1}) {  // we = 0 is west, we = 1 is east
          if (!old_idem[position_ + 2 * we]) { continue; }  // can't add generator
          Gen_type we_marking = we ? east : west;
          Idem we_idem = extend_(old_idem, we_marking);
          new_d_module.add_gen_bundle(we_idem, we_marking, gen_handle, std::string(1, "WE"[we]), alexander_maslov_weights(we_marking));  // temporary?
          
          Alg_el alg_el(we_idem, old_idem);
          new_d_module.add_coef_bundle(alg_el, we_marking, south, old_idem);
          
          if (upper_algebra_.matchings[position_] != position_ + 1) {  // curved
            std::vector< int > U_curved(upper_algebra_.n_strands, 0);
            U_curved[lower_algebra_.matchings[position_ + we]] = 1;
            Alg_el curved_alg_el(old_idem, we_idem, U_curved);
            if (curved_alg_el.is_null()) { continue; }
            new_d_module.add_coef_bundle(curved_alg_el, south, we_marking, old_idem);
          }
        }  // west, east
      }
    }
  }  // delta_0_and_1_
  
  /* \delta_2.
   */
  void delta_2_(const D_module& old_d_module, D_module& new_d_module) const {
    for (const auto& coef : old_d_module.coef_bundles()) {
      /* Calculate preliminary information */
      const Alg_el& alg_el = coef.algebra_element();
      int a1, a2, u1, u2;
      std::tie(a1, a2, u1, u2) = get_local_weights_(alg_el);
      int pre_hash_index = delta_2_table_.pre_hash_index(a1, a2, u1, u2);
      
      for (const Gen_type front_marking : {north, west, south, east}) {
        if (!extendable_(alg_el.target_idem(), front_marking)) { continue; }
        const Gen_type back_marking = delta_2_table_.positive_look_back[pre_hash_index + front_marking];
        if (!extendable_(alg_el.source_idem(), back_marking)) { continue; }
        
        const Idem new_source_idem = extend_(alg_el.source_idem(), back_marking);
        const Idem new_target_idem = extend_(alg_el.target_idem(), front_marking);
        if (new_source_idem.too_far_from(new_target_idem)) { continue; }  // incompatible idems
        
        /* Calculate new algebra element */
        std::vector<int> U_weights = alg_el.U_weights();
        int v1 = 2 * u1 + std::abs(a1);
        int v2 = 2 * u2 + std::abs(a2);
        if (back_marking == east) { --v2; }
        else if (back_marking == west) { --v1; } 
        if (front_marking == east) { ++v2; }
        else if (front_marking == west) { ++v1; }
        U_weights[position_] = v2 / 2;
        U_weights[position_ + 1] = v1 / 2;
        const Alg_el new_alg_el(new_source_idem, new_target_idem, U_weights);
        if (new_alg_el.is_null()) { continue; }  // null algebra element
        
        new_d_module.add_coef_bundle(new_alg_el, back_marking, front_marking, coef);
      }
    }
  }  // delta_2_
  
  /* \delta_3.
   * 
   * Pre-condition: we have calculated the neighbors of differential arcs.
   */
  void delta_3_(const D_module& old_d_module, D_module& new_d_module) const {
    for (const auto& front_coef : old_d_module.coef_bundles()) {
      const Alg_el& front_alg_el = front_coef.algebra_element();
      for (const auto& back_coef : front_coef.others_to_source()) {
        const Alg_el& back_alg_el = back_coef.algebra_element();
        for (const Gen_type front_marking : {north, east, south, west}) {
          if (!extendable_(front_alg_el.target_idem(), front_marking)) { continue; }
          //if (!extendable_(back_alg_el.source_idem(), south)) { continue; }  // check if I need this
          if (!coef_exists_(back_alg_el, front_alg_el, front_marking)) { continue; }
          
          const Idem& new_source_idem = back_alg_el.source_idem();  // back_marking is south
          const Idem& new_target_idem = extend_(front_alg_el.target_idem(), front_marking);
          if (new_source_idem.too_far_from(new_target_idem)) { continue; }  // algebra element is null
          
          /* Calculate new algebra element */
          auto concat_coef = back_coef.light_concatenate(front_coef);
          auto new_U_weights = concat_coef.algebra_element.U_weights();
          int a1, a2, u1, u2, b1, b2, v1, v2;
          std::tie(a1, a2, u1, u2) = get_local_weights_(back_alg_el);
          std::tie(b1, b2, v1, v2) = get_local_weights_(front_alg_el);
          int w1 = 2 * u1 + 2 * v1 + std::abs(a1) + std::abs(b1) - 1;
          int w2 = 2 * u2 + 2 * v2 + std::abs(a2) + std::abs(b2) - 1;
          if (front_marking == east) { ++w2; } // extra weight corresponding to L_2
          else if (front_marking == west) { ++w1; } // extra weight corresponding to R_1
          new_U_weights[position_] = w2 / 2;
          new_U_weights[position_ + 1] = w1 / 2;
          Alg_el new_alg_el(new_source_idem, new_target_idem, new_U_weights);
          if (new_alg_el.is_null()) { continue; }
          
          /* Make and add new differential arc */
          new_d_module.add_coef_bundle(new_alg_el, south, front_marking, concat_coef);
        }
      }
    }
  }  // delta_3_
  
  /* Auxiliary functions */
    
  /* Adapted from ComputeHFKv2/Utility.cpp, LeftRight.
   * For \delta_2 and \delta_3. 
   * Get local LR weights and U weights. This is because the DA-bimodule for a
   * crossing depends on these local weights.
   * 
   * I factored this out, but this means it gets called more than strictly
   * necessary.
   */
  std::tuple< int, int, int, int > get_local_weights_(const Alg_el& alg_el) const {
    const Idem& source_idem = alg_el.source_idem();
    const Idem& target_idem = alg_el.target_idem();
    int a1 = 0;  // will only take values in [-1, 1]
    int a2 = 0;  // will only take values in [-1, 1]
    int i = 0;
    for ( ; i <= position_; ++i) {
      a1 = a1 + source_idem[i] - target_idem[i];
    }
    a2 = a1 + source_idem[i] - target_idem[i];  // i = position_ + 1
    return std::make_tuple(a1, a2, alg_el.U_weight(position_), alg_el.U_weight(position_ + 1));
  }  // get_local_weights_
  
  /* For \delta_3
   *
   * We follow [OzsvathSzabo2018, Lemma 5.5]. Given a front marking Y, we
   * compute the mid marking I(b, Y), the front marking I(a, I(b, Y)), and the
   * product marking I(ab, Y).
   */
  bool coef_exists_(const Alg_el& back_alg_el, const Alg_el& front_alg_el, const Gen_type front_marking) const {
    if (front_marking == south) { return false; }  // front cannot be south
    
    int hash_index;
    int a1, a2, u1, u2, b1, b2, v1, v2;  // local weights
    std::tie(a1, a2, u1, u2) = get_local_weights_(back_alg_el);
    std::tie(b1, b2, v1, v2) = get_local_weights_(front_alg_el);
    
    hash_index = delta_2_table_.hash_index(b1, b2, v1, v2, front_marking);
    Gen_type mid_marking = delta_2_table_.positive_look_back[hash_index];
    if (mid_marking == null_NEW) { return false; }  // check that mid is not null
    
    hash_index = delta_2_table_.hash_index(a1, a2, u1, u2, mid_marking);
    Gen_type back_marking = delta_2_table_.positive_look_back[hash_index];
    
    hash_index = delta_2_table_.hash_index(a1 + b1, a2 + b2,
      u1 + v1 + (std::abs(a1) + std::abs(b1)) / 2,
      u2 + v2 + (std::abs(a2) + std::abs(b2)) / 2, front_marking);
    Gen_type product_marking = delta_2_table_.positive_look_back[hash_index];
    if (back_marking == product_marking) { return false; }  // check that back is not equal to product
    
    /* If product marking is null, we either have (R_1, R_2U_2^t) for east
     * or (L_2, L_1U_1^n) for west. We exclude all other cases.
     * TO DO: Understand why this is here.
     */
    if (product_marking == null_NEW and front_marking == east
      and !(a1 == 1 and a2 == 0 and u1 == 0 and u2 == 0
        and b1 == 0 and b2 == 1 and v1 == 0)) { return false; }
    else if (product_marking == null_NEW and front_marking == west
      and !(a1 == 0 and a2 == -1 and u1 == 0 and u2 == 0
        and b1 == -1 and b2 == 0 and v2 == 0)) { return false; }
    return true;
  }  // back_marking_exists_
    
  /* Taken from ComputeHFKv2/Crossing.cpp, Extendable */
  bool extendable_(const Idem& idem, const Gen_type marking) const {
    switch (marking) {
      case north: return idem[position_ + 1];
      case east: return (!idem[position_ + 1] and idem[position_ + 2]);
      case south: return !idem[position_ + 1];
      case west: return (!idem[position_ + 1] and idem[position_]);
      default: return false;
    }
  }

  /* Taken from ComputeHFKv2/Crossing.cpp, Extend */
  Idem extend_(Idem idem, const Gen_type marking) const {
    if (marking == east) {
      idem.flip(position_ + 1);
      idem.flip(position_ + 2);
    }
    else if (marking == west) {
      idem.flip(position_);
      idem.flip(position_ + 1);
    }
    return idem;  // including cases north, south
  }
  
  int position_;
  
  /* Inherited data members */
  using DA_bimodule::upper_algebra_;
  using DA_bimodule::lower_algebra_;
  
  /* Lookup table */
  static Delta_2_lookup_tables< Two_bit_table<128> > delta_2_table_;
};

template< bool sign, class D_module >
Delta_2_lookup_tables< Two_bit_table<128> > DA_bimodule_crossing< sign, D_module >::delta_2_table_;

#endif  // EVENT_CROSSING_H_