#include <string>
#include <vector>

/* Abstact base class for DA-bimodules.
 *
 * This is the best solution I came up with. Ideally, each type of DA-bimodule
 * would be an object of a DA-bimodule class. However, such an object would
 * have member function objects that need to be templated by the D-module
 * class, which is not possible.
 */
template< class D_module >
class DA_bimodule {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Algebra = typename D_module::Bordered_algebra;
  using Alg_el = typename D_module::Alg_el;
  using Coef_bundle = typename D_module::Coef_bundle;
  using Polynomial = typename D_module::Polynomial;
  using Monomial = typename Polynomial::Monomial;  // currently, pair of int
  
  DA_bimodule() { }
  
  DA_bimodule(Algebra upper_algebra, Algebra lower_algebra) :
    upper_algebra_(upper_algebra),
    lower_algebra_(lower_algebra)
  { }
  
  virtual ~DA_bimodule() = default;
  
  friend D_module box_tensor_product(const D_module& old_d_module, const DA_bimodule& da_bimodule) {
    D_module new_d_module = da_bimodule.tensor_generators_(old_d_module);
    new_d_module.lock_generators(old_d_module, da_bimodule.weights_(), da_bimodule.labels_());
    da_bimodule.tensor_coefficients_(old_d_module, new_d_module);
    new_d_module.lock_coefficients();
    return new_d_module;
  }
  
 protected:
  /* Double Alexander weight and Maslov weight */
  virtual std::vector< Monomial > weights_() const = 0;
  
  /* Given a generator type, give the corresponding LaTeX math label. This is
   * based on notations given by mathematicians.
   */
  virtual std::vector< std::string > labels_() const = 0;
  
  virtual D_module tensor_generators_(const D_module& old_d_module) const = 0;
  virtual D_module& tensor_coefficients_(const D_module& old_d_module, D_module& new_d_module) const = 0;
  
  Algebra upper_algebra_;  // These members really belong to the ABC
  Algebra lower_algebra_;
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
  using Coef_bundle = typename DA_bimodule::Coef_bundle;
  
  template< class Knot_slice_handle >
  DA_bimodule_local_maximum(const Knot_slice_handle& slice, Algebra upper_algebra, Algebra lower_algebra) :
    position_(slice->position()),
    DA_bimodule(upper_algebra, lower_algebra)  // inherited constructor
  { }
  
  ~DA_bimodule_local_maximum() = default;
  
 private:
  enum : Gen_type {
    X,
    Y,
    Z
  };
  
  std::vector< Monomial > weights_() const override {
    return std::vector< Monomial >(3, {0, 0});
  }
  
  std::vector< std::string > labels_() const override {
    std::vector< std::string > labels(3);
    std::string symbols;
    if (lower_algebra_.orientations[position_]) {  // clock
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
  
  
  D_module tensor_generators_(const D_module& old_d_module) const override {
    D_module new_d_module;
    delta_0_(old_d_module, new_d_module);
    return new_d_module;
  }
  
  D_module& tensor_coefficients_(const D_module& old_d_module, D_module& new_d_module) const override{
    delta_1_(old_d_module, new_d_module);
    delta_2_(old_d_module, new_d_module);
    return new_d_module;
  }
  
  
  void delta_0_(const D_module& old_d_module, D_module& new_d_module) const {
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
  
  void delta_1_(const D_module& old_d_module, D_module& new_d_module) const {
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
  void delta_2_(const D_module& old_d_module, D_module& new_d_module) const {
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
    else if (marking == Y or marking == Z) {  // x -> 01x
      idem.insert(position_, {0, 1});
    }
    return idem;
  }
  
  int position_;
  
  using DA_bimodule::upper_algebra_;
  using DA_bimodule::lower_algebra_;
};