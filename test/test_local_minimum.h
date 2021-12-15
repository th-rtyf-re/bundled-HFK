#ifndef TEST_LOCAL_MINIMUM_H_
#define TEST_LOCAL_MINIMUM_H_

#include <string>
#include <utility>  // pair
#include <vector>

#ifdef VERBOSE
#include <iostream>
#endif  // VERBOSE

#include <boost/any.hpp>

template< class D_module, class Morse_event_options >
class Local_minimum {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Algebra = typename D_module::Bordered_algebra;
  using Coef_bundle = typename D_module::Coef_bundle;
  using Weights = typename D_module::Weights;
  
  Local_minimum(const std::vector< typename Morse_event_options::Parameter_type >& args) :
    position_(args.empty() ? 0 : Morse_event_options::template parameter_cast< int >(args[0])) 
  {
    if (position_ != 0) {
      std::cout << "[lm] Warning: local minimum not in position 0."
                << " Behavior is undefined..."
                << std::endl;
    }
  }
  
  /* Topological methods */
  
  std::vector< int > lower_matchings(std::vector< int > matchings) const {
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
    matchings.erase(
      matchings.begin() + position_,
      matchings.begin() + position_ + 2
    );
    return matchings;
  }
  
  std::vector< bool > upper_orientations(
    std::vector< bool > orientations,
    const std::vector< int >& upper_matchings
  ) const {
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
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const {
    return {margins.first + 1, margins.second + 1};
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(
    const std::pair< int, int >& margins,
    const std::pair< int, int >& n_strands
  ) const {
    return std::string(margins.first - 1, '0')
      + std::string(position_, 'l')
      + "u"
      + std::string(n_strands.first - position_ - 2, 'r')
      + std::string(margins.second - 1, '0');
  }
  
  /* Algebraic methods */
  
  std::vector< Weights > get_weights(const Algebra&, const Algebra&) const {
    return std::vector< Weights >(2, {0, 0});
  }
  
  std::vector< std::string > get_labels(
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    std::vector< std::string > labels(2);
    std::string symbols;
    if (lower_algebra.orientations[position_]) {  // clock
      symbols = "XY";
    }
    else {  // trig
      symbols = "xyz";
    }
    
    for (auto type : {XL1, YR2}) {
      labels[type] = std::string(1, symbols[type])
        + "_{" + std::to_string(position_) + "}"
        + std::string(1, "LR"[type])
        + "_{" + std::to_string(position_ + type) + "}";
    }
    return labels;
  }
  
  D_module tensor_generators(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra&,
    const Algebra&
  ) const {
    delta_0_(new_d_module, old_d_module);
    return new_d_module;
  }
  
  D_module& tensor_coefficients(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    delta_2_(new_d_module, old_d_module, upper_algebra, lower_algebra);
    delta_geq_4_(new_d_module, old_d_module, upper_algebra, lower_algebra);
    return new_d_module;
  }
    
#ifdef VERBOSE
  friend std::ostream& operator<<(
    std::ostream& os,
    const Local_minimum& morse_event
  ) {
    os << "local minimum at position " << morse_event.position_;
    return os;
  }
#endif  // VERBOSE
  
 private:
  enum {  // names from [OzsvathSzabo2019, Section 7.2]
    XL1,  // not actually used
    YR2
  };
  
  void delta_0_(D_module& new_d_module, const D_module& old_d_module) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      const Idem& old_idem = old_d_module.idem(gen_handle);
      
      if (extendable_(old_idem, YR2)) {
        Idem yr2_idem = extend_(old_idem, YR2);
        new_d_module.add_gen_bundle(yr2_idem, YR2, gen_handle);
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
  void delta_2_(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    for (const auto& coef : old_d_module.coef_bundles()) {
      if (
        old_d_module.U_weights(coef)[0] == 0
        and extendable_(old_d_module.source_idem(coef), YR2)
        and extendable_(old_d_module.target_idem(coef), YR2)
      ) {
        /* Add curved weight*/
        std::vector< int > new_U_weights = old_d_module.U_weights(coef);
        new_U_weights[upper_algebra.matchings[0]] += new_U_weights[1];
        
        shorten_and_finish_(old_d_module, new_d_module, new_U_weights, coef);
      }
    }
  }
  
  /* \delta_{\geq 4}
   * Assume position_ is 0.
   */
  void delta_geq_4_(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    /* Construct lists of coefficients
     * As is the case everywhere else, algebra indexes start from 0, not 1
     */
    std::vector< Coef_bundle > L_1, R_1, U_0, U_1;
    for (Coef_bundle coef : old_d_module.coef_bundles()) {
      bool s1 = old_d_module.source_idem(coef)[1]; // index 1 of source and
      bool t1 = old_d_module.target_idem(coef)[1]; // target idempotents
      int u0 = old_d_module.U_weights(coef)[0]; // powers of U_0 and U_1
      int u1 = old_d_module.U_weights(coef)[1];
      
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
    std::vector< Coef_bundle > sequences_back =
      concatenate_groups_(old_d_module, L_1, U_0);
    const std::vector< Coef_bundle > sequences_mid =
      concatenate_groups_(old_d_module, U_1, U_0);
    std::vector< Coef_bundle > sequences =
      concatenate_groups_(old_d_module, sequences_back, R_1);
    
    /* This loop ends for mathematical reasons. If the input is not a knot,
     * there is no guarantee that this ends.
     */
    for (int n_coefs = 1; !sequences.empty(); ++n_coefs) {
      for (const auto& coef : sequences) {
        std::vector< int > new_U_weights = old_d_module.U_weights(coef);
        new_U_weights[upper_algebra.matchings[0]] += new_U_weights[1] - n_coefs;
        new_U_weights[upper_algebra.matchings[1]] += new_U_weights[0] - n_coefs;
        shorten_and_finish_(old_d_module, new_d_module, new_U_weights, coef);
      }
      sequences_back =
        concatenate_groups_(old_d_module, sequences_back, sequences_mid);
      sequences = concatenate_groups_(old_d_module, sequences_back, R_1);
      // debug, I think
      if (n_coefs >= 10) { std::cout << "[lm] Giving up!" << std::endl; return; }
    }
  }
  
  /* Concatenate two groups of coefficients.
   * TO DO: sort both sets by... suffix of target (resp. source)? I don't think
   * this works with agglomerated coefs, though
   */
  std::vector< Coef_bundle > concatenate_groups_(
    const D_module& old_d_module,
    const std::vector< Coef_bundle >& back_group,
    const std::vector< Coef_bundle >& front_group
  ) const {
    std::vector< Coef_bundle > result;
    for (const Coef_bundle& back_coef : back_group) {
      for (const Coef_bundle& front_coef : front_group) {
        if (
          old_d_module.source_idem(back_coef).too_far_from(
            old_d_module.target_idem(front_coef)
          )
        ) {
          continue;
        }
        if (!old_d_module.compatible(back_coef, front_coef)) { continue; }
        Coef_bundle concat_coef =
          old_d_module.concatenate(back_coef, front_coef);
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
  void shorten_and_finish_(
    const D_module& old_d_module,
    D_module& new_d_module,
    std::vector< int >& new_U_weights,
    const Coef_bundle& old_coef
  ) const {
    /* Shorten idempotents and U weight vector */
    Idem new_source_idem = old_d_module.source_idem(old_coef);
    Idem new_target_idem = old_d_module.target_idem(old_coef);
    new_source_idem.erase(1, 2);
    new_target_idem.erase(1, 2);
    new_U_weights.erase(new_U_weights.begin(), new_U_weights.begin() + 2);
    
    /* Add composite coefficient, if possible */
    if (new_source_idem.too_far_from(new_target_idem)) { return; }
    auto alg_el =
      new_d_module.alg_el(new_source_idem, new_target_idem, new_U_weights);
    new_d_module.add_coef_bundle(alg_el, YR2, YR2, old_coef, old_d_module);
  }
  
  bool extendable_(const Idem& idem, const Gen_type marking) const {
    if (marking == YR2) {
      return (idem[2] and !idem[1] and !idem[0]); // don't need !idem[0]
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
};

#endif  // TEST_LOCAL_MINIMUM_H_