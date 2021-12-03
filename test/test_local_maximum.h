#include <string>
#include <utility>  // pair
#include <vector>

#include <boost/any.hpp>

/* DA-bimodule for a local maximum.
 * 
 * Adapted from ComputeHFKv2/Max.cpp
 */
template< class D_module >
class Local_maximum {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Algebra = typename D_module::Bordered_algebra;
  using Alg_el = typename D_module::Alg_el;
  using Coef_bundle = typename D_module::Coef_bundle;
  using Weights = typename D_module::Weights;  // currently, pair of int
  
  Local_maximum(const std::vector< boost::any >& args) :
    position_(args.empty()? 0 : boost::any_cast< int >(args[0]))
  { }
  
  
  /* Topological methods */
  
  std::vector< int > lower_matchings(std::vector< int > matchings) const {
    for (int match : matchings) {
      if (match >= position_) {
        match += 2;
      }
    }
    matchings.insert(matchings.begin() + position_, {position_ + 1, position_});
    return matchings;
  }
  
  std::vector< bool > upper_orientations(std::vector< bool > orientations, const std::vector< int >&) const {
    orientations.erase(orientations.begin() + position_,
    orientations.begin() + position_ + 2);
    return orientations;
  }
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const {
    return {margins.first - 1, margins.second - 1};
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(const std::pair< int, int >& margins, const std::pair< int, int >& n_strands) const {
    return std::string(margins.first, '0')
      + std::string(position_, 'r')
      + "a"
      + std::string(n_strands.second - position_ - 2, 'l')
      + std::string(margins.second, '0');
  }
  
  /* Algebraic methods */
  
  std::vector< Weights > get_weights(const Algebra&, const Algebra&) const {
    return std::vector< Weights >(3, {0, 0});
  }
  
  std::vector< std::string > get_labels(const Algebra& upper_algebra, const Algebra& lower_algebra) const {
    std::vector< std::string > labels(3);
    std::string symbols;
    if (lower_algebra.orientations[position_]) {  // clock
      symbols = "XYZ";
    }
    else {  // trig
      symbols = "xyz";
    }
    
    for (auto type : {X, Y, Z}) {
      labels[type] = std::string(1, symbols[type])
                     + "_{" + std::to_string(position_) + "}";
    }
    return labels;
  }
  
  D_module tensor_generators(const D_module& old_d_module, const Algebra&, const Algebra&) const {
    D_module new_d_module;
    delta_0_(new_d_module, old_d_module);
    return new_d_module;
  }
  
  D_module& tensor_coefficients(D_module& new_d_module, const D_module& old_d_module, const Algebra&, const Algebra&) const {
    delta_1_(new_d_module, old_d_module);
    delta_2_(new_d_module, old_d_module);
    return new_d_module;
  }
  
 private:
  enum {
    X,
    Y,
    Z
  };
  
  void delta_0_(D_module& new_d_module, const D_module& old_d_module) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      const Idem& old_idem = old_d_module.idem(gen_handle);
      
      if (old_idem[position_]) {  // X and Y
        Idem x_idem = extend_(old_idem, X);
        Idem y_idem = extend_(old_idem, Y);
        new_d_module.add_gen_bundle(x_idem, X, gen_handle);
        new_d_module.add_gen_bundle(y_idem, Y, gen_handle);
      }
      
      else {  // Z
        Idem z_idem = extend_(old_idem, Z);
        new_d_module.add_gen_bundle(z_idem, Z, gen_handle);
      }
    }
  }
  
  void delta_1_(D_module& new_d_module, const D_module& old_d_module) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      const Idem& old_idem = old_d_module.idem(gen_handle);
      
      if (old_idem[position_]) {  // X and Y
        Idem x_idem = extend_(old_idem, X);
        Idem y_idem = extend_(old_idem, Y);
        new_d_module.add_coef_bundle(Alg_el(x_idem, y_idem), X, Y, old_idem);
        new_d_module.add_coef_bundle(Alg_el(y_idem, x_idem), Y, X, old_idem);
      }
    }
  }
  
  /* See [OzsvathSzabo2018, Lemma 8.1] */
  void delta_2_(D_module& new_d_module, const D_module& old_d_module) const {
    for (const auto& coef : old_d_module.coef_bundles()) {
      const Idem& back_idem = old_d_module.source_idem(coef);
      const Idem& front_idem = old_d_module.target_idem(coef);
      std::vector< int > new_U_weights = old_d_module.U_weights(coef);
      
      int a1, a2;
      std::tie(a1, a2) = get_local_weights_(back_idem, front_idem);
      new_U_weights.insert(new_U_weights.begin() + position_, {0, 0});
      
      std::vector< std::pair< Gen_type, Gen_type > > to_compose;  // of length 1 or 2
      
      if (a1 == 1 and a2 == 1) {  // R_2R_1
        to_compose.emplace_back(Y, X);
      }
      else if (a1 == -1 and a2 == -1) {  // L_1L_2
        to_compose.emplace_back(X, Y);
      }
      else if (a1 == 1 and a2 == 0) {  // R_1
        to_compose.emplace_back(Z, X);
      }
      else if (a1 == -1 and a2 == 0) {  // L_1
        to_compose.emplace_back(X, Z);
      }
      else if (a1 == 0 and a2 == 1) {  // R_2
        to_compose.emplace_back(Y, Z);
      }
      else if (a1 == 0 and a2 == -1) {  // L_2
        to_compose.emplace_back(Z, Y);
      }
      else if (a1 == 0 and a2 == 0) {
        if (back_idem[position_]) {  // X, Y
          to_compose.emplace_back(X, X);
          to_compose.emplace_back(Y, Y);
        }
        else {  // Z
          to_compose.emplace_back(Z, Z);
        }
      }
      
      for (const auto& value_pair : to_compose) {
        const Gen_type back_marking = value_pair.first;
        const Gen_type front_marking = value_pair.second;
        const Idem& new_source_idem = extend_(back_idem, back_marking);
        const Idem& new_target_idem = extend_(front_idem, front_marking);
        if (new_source_idem.too_far_from(new_target_idem)) { continue; }
        const Alg_el new_alg_el(new_source_idem, new_target_idem, new_U_weights);
        if (new_alg_el.is_null()) { continue; }
        new_d_module.add_coef_bundle(new_alg_el, back_marking, front_marking, coef, old_d_module);
      }
    }
  }  // delta_2_
  
  /* Get LR weights at position_ - 1 and position_. This corresponds to the
   * local LR weights of the upper algebra.
   */
  std::tuple< int, int > get_local_weights_(const Idem& source_idem, const Idem& target_idem) const {
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
    else if (marking == Y or marking == Z) {  // * -> 01*
      idem.insert(position_, {0, 1});
    }
    return idem;
  }
  
  const int position_;
};