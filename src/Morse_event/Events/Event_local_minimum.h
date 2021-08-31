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

#ifndef EVENT_LOCAL_MINIMUM_H_
#define EVENT_LOCAL_MINIMUM_H_

#include <string>
#include <tuple>  // tuple, tie
#include <vector>

#include "../DA_bimodule.h"
#include "../Knot_slice.h"

/* Knot slice for a local minimum.
 */
class Knot_slice_local_minimum : public Knot_slice {
 public:
  Knot_slice_local_minimum(int upper_n_strands, int position) :
    Knot_slice(upper_n_strands, upper_n_strands - 2, position)
  { }
  
  ~Knot_slice_local_minimum() = default;
  
  Morse_event event_id() const override {
    return Morse_event::local_minimum;
  }
  
  std::vector< int > lower_matchings(std::vector< int > matchings) const override {
    for (int& match : matchings) {
      if (match == position_) {
        match = matchings[position_ + 1];
      }
      else if (match == position_ + 1) {
        match = matchings[position_];
      }
      else if (match >= position_ + 2) {
        match = match - 2;
      }
    }
    matchings.erase(matchings.begin() + position_, matchings.begin() + position_ + 2);
    return matchings;
  }
  
  std::vector< bool > upper_orientations(std::vector< bool > orientations, const std::vector< int >& upper_matchings) const override {
    int left_match = upper_matchings[position_];
    int right_match = upper_matchings[position_ + 1];
    if (left_match >= position_ + 2) {
      left_match -= 2;
    }
    if (right_match >= position_ + 2) {
      right_match -= 2;
    }
    bool left_ori = orientations[left_match];
    bool right_ori = orientations[right_match];
    orientations.insert(orientations.begin() + position_, {right_ori, left_ori});
    return orientations;
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(const std::pair< int, int >& margins) const override {
    return std::string(margins.first - 1, '0')
      + std::string(position_, 'l')
      + "u"
      + std::string(upper_n_strands_ - position_ - 2, 'r')
      + std::string(margins.second - 1, '0');
  }
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const override {
    return {margins.first + 1, margins.second + 1};
  }
  
 private:
  using Knot_slice::upper_n_strands_;
  using Knot_slice::lower_n_strands_;
  using Knot_slice::position_;
};

/* DA-bimodule for a local minimum.
 * 
 * Adapted from ComputeHFKv2/Min.cpp
 * 
 * Following Ozsváth and Szabó, this assumes that the local minimum is in the
 * leftmost position.
 */
template< class D_module >
class DA_bimodule_local_minimum : public DA_bimodule< D_module > {
 public:
  using DA_bimodule = DA_bimodule< D_module >;
  using Idem = typename DA_bimodule::Idem;
  using Gen_type = typename DA_bimodule::Gen_type;
  using Algebra = typename DA_bimodule::Algebra;
  using Alg_el = typename DA_bimodule::Alg_el;
  using Coef_bundle = typename DA_bimodule::Coef_bundle;
  using Monomial = typename DA_bimodule::Monomial;
  
  DA_bimodule_local_minimum(const std::unique_ptr< Knot_slice >& slice, Algebra upper_algebra, Algebra lower_algebra) :
    position_(slice->position()),
    DA_bimodule(upper_algebra, lower_algebra)
  {
    if (position_ != 0) {
      std::cout << "[lm] Warning: local minimum not in position 0."
                << " Behavior is undefined..."
                << std::endl;
    }
  }
  
  ~DA_bimodule_local_minimum() = default;
  
  Monomial alexander_maslov_weights(const Gen_type generator_type) const override {
    return {0, 0};
  }
  
  std::string get_label(const Gen_type generator_type) const override {
    if (!upper_algebra_.orientations[position_]) {  // clock
      return std::string(1, "XY"[generator_type])
        + "_{" + std::to_string(position_) + "}"
        + std::string(1, "LR"[generator_type])
        + "_{" + std::to_string(position_ + generator_type) + "}";
    }
    else {  // trig
      return std::string(1, "xy"[generator_type])
        + "_{" + std::to_string(position_) + "}"
        + std::string(1, "LR"[generator_type])
        + "_{" + std::to_string(position_ + generator_type) + "}";
    }
  }
  
  D_module box_tensor_product(const D_module& old_d_module) const override {
    D_module new_d_module;
    delta_0_(old_d_module, new_d_module);
    
    delta_2_(old_d_module, new_d_module);
    
    delta_geq_4_(old_d_module, new_d_module);
    return new_d_module;
  }
  
 private:
  enum : Gen_type {  // names from [OzsvathSzabo2019, Section 7.2]
    XL1,  // not actually used
    YR2
  };
  
  void delta_0_(const D_module& old_d_module, D_module& new_d_module) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      const Idem& old_idem = old_d_module.idem(gen_handle);
      
      if (extendable_(old_idem, YR2)) {
        Idem yr2_idem = extend_(old_idem, YR2);
        new_d_module.add_gen_bundle(yr2_idem, YR2, gen_handle, "m");
      }
    }
  }
  
  /* \delta_2.
   * The coefficient created corresponds to
   * 
   *     \delta_2(YR_2, U_2^n) = U_\alpha^n \otimes YR_2
   * 
   * Assume that position_ is 0.
   */
  void delta_2_(const D_module& old_d_module, D_module& new_d_module) const {
    for (const auto& coef : old_d_module.coef_bundles()) {
      const Alg_el& alg_el = coef.algebra_element();
      
      if (alg_el.U_weight(0) == 0 and extendable_(alg_el.source_idem(), YR2)
          and extendable_(alg_el.target_idem(), YR2)) {
        /* Add curved weight*/
        std::vector< int > new_U_weights = alg_el.U_weights();
        new_U_weights[upper_algebra_.matchings[0]] += new_U_weights[1];
        
        shorten_and_finish_(new_U_weights, coef, new_d_module);
      }
    }
  }
  
  /* \delta_{\geq 4}
   * Assume position_ is 0.
   */
  void delta_geq_4_(const D_module& old_d_module, D_module& new_d_module) const {
    /* Construct lists of coefficients
     * As is the case everywhere else, algebra indexes start from 0, not 1
     */
    std::vector< Coef_bundle > L_1, R_1, U_0, U_1;
    for (Coef_bundle coef : old_d_module.coef_bundles()) {
      bool s1, t1;  // index 1 of source and target idempotents
      int u0, u1;  // powers of U_0 and U_1
      std::tie(s1, t1, u0, u1) = get_local_info_(coef);
      if (!s1 and t1) {  // L_1
        L_1.push_back(coef);
      }
      else if (s1 and !t1) {  // R_1
        R_1.push_back(coef);
      }
      else if (s1 and t1 and u0 > 0 and u1 == 0) {  // U_0
        U_0.push_back(coef);
      }
      else if (s1 and t1 and u0 == 0 and u1 > 0) {  // U_1
        U_1.push_back(coef);
      }
    }
    
    /* Initialize building blocks of sequences */
    std::vector< Coef_bundle > sequences_back = concatenate_groups_(L_1, U_0);
    const std::vector< Coef_bundle > sequences_mid = concatenate_groups_(U_1, U_0);
    std::vector< Coef_bundle > sequences = concatenate_groups_(sequences_back, R_1);
    
    /* This loop ends for mathematical reasons. If the input is not a knot,
     * there is no guarantee that this ends.
     */
    for (int n_coefs = 1; !sequences.empty(); ++n_coefs) {
      for (const auto& coef : sequences) {
        const Alg_el& alg_el = coef.algebra_element();
        std::vector< int > new_U_weights = alg_el.U_weights();
        new_U_weights[upper_algebra_.matchings[0]] += new_U_weights[1] - n_coefs;
        new_U_weights[upper_algebra_.matchings[1]] += new_U_weights[0] - n_coefs;
        shorten_and_finish_(new_U_weights, coef, new_d_module);
      }
      sequences_back = concatenate_groups_(sequences_back, sequences_mid);
      sequences = concatenate_groups_(sequences_back, R_1);
      if (n_coefs >= 10) { std::cout << "[lm] Giving up!" << std::endl; return; }
    }
  }
  
  /* Concatenate two groups of coefficients.
   * TO DO: sort both sets by... suffix of target (resp. source)? I don't think
   * this works with agglomerated coefs, though
   */
  std::vector< Coef_bundle > concatenate_groups_(const std::vector< Coef_bundle >& back_group, const std::vector< Coef_bundle >& front_group) const {
    std::vector< Coef_bundle > result;
    for (const Coef_bundle& back_coef : back_group) {
      for (const Coef_bundle& front_coef : front_group) {
        if (back_coef.source_idem().too_far_from(front_coef.target_idem())) { continue; }
        if (!back_coef.compatible(front_coef)) { continue; }
        Coef_bundle concat_coef = back_coef * front_coef;
        result.push_back(concat_coef);
      }
    }
    /* modulo 2 */
    std::vector< bool > canceled(result.size(), false);
    for (int i = 0; i < result.size(); ++i) {
      for (int j = i + 1; j < result.size(); ++j) {
        if (!canceled[j] and result[i] == result[j]) {
          canceled[i] = true;
          canceled[j] = true;
          break;
        }
      }
    }
    std::vector< Coef_bundle > new_result;
    for (int i = 0; i < result.size(); ++i) {
      if (!canceled[i]) {
        new_result.push_back(result[i]);
      }
    }
    return new_result;
  }
  
  /* For \delta_{\geq 2}.
   * This auxiliary function does two things, which can't really be separated:
   * shorten idempotents and U weights, and then add composite coefficient.
   */
  void shorten_and_finish_(std::vector< int >& new_U_weights, const Coef_bundle& old_coef, D_module& new_d_module) const {
    /* Shorten idempotents and U weight vector */
    Idem new_source_idem = old_coef.algebra_element().source_idem();
    Idem new_target_idem = old_coef.algebra_element().target_idem();
    new_source_idem.erase(1, 2);
    new_target_idem.erase(1, 2);
    new_U_weights.erase(new_U_weights.begin(), new_U_weights.begin() + 2);
    
    /* Add composite coefficient, if possible */
    if (new_source_idem.too_far_from(new_target_idem)) { return; }
    Alg_el new_alg_el(new_source_idem, new_target_idem, new_U_weights);
    if (new_alg_el.is_null()) { return; }
    new_d_module.add_coef_bundle(new_alg_el, YR2, YR2, old_coef);
  }
  
  std::tuple< bool, bool, int, int > get_local_info_(const Coef_bundle& coef) const {
    const Alg_el& alg_el = coef.algebra_element();
    return std::make_tuple(alg_el.source_idem()[1], alg_el.target_idem()[1],
                           alg_el.U_weight(0), alg_el.U_weight(1));
  }  // get_local_info_
  
  bool extendable_(const Idem& idem, const Gen_type marking) const {
    if (marking == YR2) {
      return (idem[2] and !idem[1] and !idem[0]); // technically don't need !idem[0]
    }
    return false;
  }
  
  Idem extend_(Idem idem, const Gen_type marking) const {
    if (marking == YR2) {
      idem.erase(1, 2);
    }
    return idem;
  }
  
  int position_;  // = 0
  
  using DA_bimodule::upper_algebra_;
  using DA_bimodule::lower_algebra_;
};

#endif  // EVENT_LOCAL_MINIMUM_H_