#ifndef TEST_D_MODULE_REVERSE_VIEW_
#define TEST_D_MODULE_REVERSE_VIEW_

/* Reverse view of a forest. Should also be a model of the D-module concept
 */
template< class D_module >
class Reverse_D_module {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Bordered_algebra = typename D_module::Bordered_algebra;
  using Alg_el = typename D_module::Alg_el;
  using Weights = typename D_module::Weights;
  
  using Arc = typename D_module::Arc;
  
  using Arc_container = typename D_module::Arc_container;
  
  using Arc_reference = typename D_module::Arc_reference;
  
  using Root_handle_container = typename D_module::Root_handle_container;
  using Gen_bundle_handle = typename D_module::Gen_bundle_handle;
  using Coef_bundle = typename D_module::Coef_bundle;
  
  Reverse_D_module(D_module& forest) :
    d_module_(forest)
  { }
    
  const Root_handle_container& gen_bundle_handles() const {
    return d_module_.gen_bundle_handles();
  }
  
  const Arc_container& coef_bundles() const {
    return d_module_.coef_bundles();
  }
  
  Idem idem(const Gen_bundle_handle& gen_bundle_handle) const {
    return d_module_.idem(gen_bundle_handle);
  }
  
  Idem source_idem(const Arc& arc) const {
    return d_module_.target_idem(arc);
  }
  
  Idem target_idem(const Arc& arc) const {
    return d_module_.source_idem(arc);
  }
  
  std::vector< int > U_weights(const Arc& arc) const {
    return d_module_.U_weights(arc);
  }
  
  void set_as_trivial() {
    d_module_.set_as_trivial();
  }
  
  /* Operations on nodes */
  
  void add_gen_bundle(Idem new_idem, Gen_type new_type, Gen_bundle_handle root_handle) {
    d_module_.add_gen_bundle(new_idem, new_type, root_handle);
  }
  
  void add_gen_bundle(Idem new_idem) {
    d_module_.add_gen_bundle(new_idem);
  }
  
  void lock_generators(const D_module& old_d_module_,
                     const std::vector< Weights >& first_layer_weights,
                     const std::vector< std::string >& first_layer_labels) {
    d_module_.lock_generators(old_d_module_, first_layer_weights, first_layer_labels);
  }
  
  /* Arc creation */
  
  void add_coef_bundle(Alg_el new_value,
                   Gen_type back_marking,
                   Gen_type front_marking,
                   const Arc& old_arc,
                   const D_module& old_d_module_) {
    d_module_.add_coef_bundle(new_value, front_marking, back_marking, old_arc, old_d_module_);
  }
  
  void add_coef_bundle(Alg_el new_value,
                   Gen_type back_marking,
                   Gen_type front_marking,
                   const Idem& old_idem) {
    d_module_.add_coef_bundle(new_value, front_marking, back_marking, old_idem);
  }
  
  void lock_coefficients() {
    d_module_.lock_coefficients();
  }
  
  /* Views on arcs */
  
  std::vector< Arc_reference > arcs_to_source(const Arc& arc) const {
    return d_module_.arcs_from_target(arc);
  }
  
  std::vector< Arc_reference > arcs_from_target(const Arc& arc) const {
    return d_module_.arcs_to_source(arc);
  }
  
  std::vector< Arc_reference > arcs_from_source(const Arc& arc) const {
    return arcs_to_target(arc);
  }
  
  std::vector< Arc_reference > arcs_to_target(const Arc& arc) const {
    return arcs_from_source(arc);
  }
  
  /* Operations on arcs */
  
  bool compatible(const Arc& back_arc, const Arc& front_arc) const {
    return d_module_.compatible(front_arc, back_arc);
  }
  
  Arc concatenate(const Arc& back_arc, const Arc& front_arc) const {
    return d_module_.concatenate(front_arc, back_arc);
  }
  
  void reduce() {
    d_module_.reduce();
  }
  
  template< class Polynomial >
  Polynomial poincare_polynomial(const Idem& idem = Idem("0")) const {
    return d_module_.poincare_polynomial(idem);
  }
  
  /* I/O interface */
  
  friend std::ostream& operator<<(std::ostream& os, const Reverse_D_module& reverse_view_) {
    return os << reverse_view_.d_module_;
  }

#ifdef DRAW
  /* TeXify */
  void TeXify(std::ofstream& write_file) const {
    d_module_.TeXify(write_file);
  }
#endif  // DRAW
  
 private:
  D_module& d_module_;
};

template< class D_module >
Reverse_D_module< D_module > reverse_view(D_module& forest) {
  return Reverse_D_module< D_module >(forest);
}

#endif  // TEST_D_MODULE_REVERSE_VIEW_