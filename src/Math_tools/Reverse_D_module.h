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

#ifndef D_MODULE_REVERSE_VIEW_
#define D_MODULE_REVERSE_VIEW_

/* Reverse view of a D-module.
 * 
 * Reverses the directions of all arcs.
 * 
 * Mathematically speaking, this gives the dual of a D-module.
 * 
 * This class satisfies the D-module concept, and in fact is probably the
 * reference for the concept.
 */
template< class D_module >
class Reverse_D_module {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Bordered_algebra = typename D_module::Bordered_algebra;
  using Alg_el = typename D_module::Alg_el;
  using Weights = typename D_module::Weights;
  
  using Gen_bundle_handle_container =
    typename D_module::Gen_bundle_handle_container;
  using Gen_bundle_handle = typename D_module::Gen_bundle_handle;
  using Coef_bundle = typename D_module::Coef_bundle;
  using Coef_bundle_container = typename D_module::Coef_bundle_container;
  using Coef_bundle_reference = typename D_module::Coef_bundle_reference;
  
  Reverse_D_module(D_module& d_module) :
    d_module_(d_module)
  { }
  
  Reverse_D_module(const D_module& d_module) :
    d_module_(const_cast< D_module& >(d_module))
  { }
    
  const Gen_bundle_handle_container& gen_bundle_handles() const {
    return d_module_.gen_bundle_handles();
  }
  
  const Coef_bundle_container& coef_bundles() const {
    return d_module_.coef_bundles();
  }
  
  Idem idem(const Gen_bundle_handle& gen_bundle_handle) const {
    return d_module_.idem(gen_bundle_handle);
  }
  
  Idem source_idem(const Coef_bundle& coef) const {
    return d_module_.target_idem(coef);
  }
  
  Idem target_idem(const Coef_bundle& coef) const {
    return d_module_.source_idem(coef);
  }
  
  std::vector< int > U_weights(const Coef_bundle& coef) const {
    return d_module_.U_weights(coef);
  }
  
  int U_weight(const Coef_bundle& coef, int position) const {
    return d_module_.U_weight(coef, position);
  }
  
  void set_as_trivial() {
    d_module_.set_as_trivial();
  }
  
  /* Operations on nodes */
  
  void add_gen_bundle(
    Idem new_idem,
    Gen_type new_type,
    Gen_bundle_handle root_handle
  ) {
    d_module_.add_gen_bundle(new_idem, new_type, root_handle);
  }
  
  void add_gen_bundle(Idem new_idem) {
    d_module_.add_gen_bundle(new_idem);
  }
  
  void lock_generators(
    const Reverse_D_module< D_module >& old_d_module,
    const std::vector< Weights >& first_layer_weights,
    const std::vector< std::string >& first_layer_labels
  ) {
    d_module_.lock_generators(
      old_d_module.d_module_,
      first_layer_weights,
      first_layer_labels
    );
  }
  
  /* Coef_bundle creation */
  
  template< class ...Args >
  Alg_el alg_el(
    const Idem source_idem,
    const Idem target_idem,
    Args&&... args
  ) const {
    return Alg_el(target_idem, source_idem, args...);
  }
  
  void add_coef_bundle(
    const Alg_el& new_value,
    const Gen_type back_marking,
    const Gen_type front_marking,
    const Coef_bundle& old_coef,
    const Reverse_D_module< D_module >& old_d_module
  ) {
    d_module_.add_coef_bundle(
      new_value,
      front_marking,
      back_marking,
      old_coef,
      old_d_module.d_module_
    );
  }
  
  void add_coef_bundle(
    const Alg_el& new_value,
    const Gen_type back_marking,
    const Gen_type front_marking,
    const Idem& old_idem
  ) {
    d_module_.add_coef_bundle(
      new_value,
      front_marking,
      back_marking,
      old_idem
    );
  }
  
  void lock_coefficients() {
    d_module_.lock_coefficients();
  }
  
  /* Views on coefs */
  
  const std::vector< Coef_bundle_reference >
  others_to_source(const Coef_bundle& coef) const {
    return d_module_.others_from_target(coef);
  }
  
  const std::vector< Coef_bundle_reference >
  others_from_target(const Coef_bundle& coef) const {
    return d_module_.others_to_source(coef);
  }
  
  const std::vector< Coef_bundle_reference >
  others_from_source(const Coef_bundle& coef) const {
    return others_to_target(coef);
  }
  
  const std::vector< Coef_bundle_reference >
  others_to_target(const Coef_bundle& coef) const {
    return others_from_source(coef);
  }
  
  /* Operations on coefs */
  
  bool compatible(
    const Coef_bundle& back_coef,
    const Coef_bundle& front_coef
  ) const {
    return d_module_.compatible(front_coef, back_coef);
  }
  
  Coef_bundle concatenate(
    const Coef_bundle& back_coef,
    const Coef_bundle& front_coef
  ) const {
    return d_module_.concatenate(front_coef, back_coef);
  }
  
  void reduce() {
    d_module_.reduce();
  }
  
  template< class Polynomial >
  Polynomial poincare_polynomial(const Idem& idem = Idem("0")) const {
    return d_module_.poincare_polynomial(idem);
  }
  
  /* I/O interface */
  
  friend std::ostream& operator<<(
    std::ostream& os,
    const Reverse_D_module& reverse_view_
  ) {
    return os << reverse_view_.d_module_;
  }

#ifdef BUNDLED_HFK_DRAW_
  /* TeXify */
  void TeXify(std::ofstream& write_file) const {
    d_module_.TeXify(write_file);
  }
#endif  // BUNDLED_HFK_DRAW_
  
 private:
  D_module& d_module_;
};

template< class D_module >
Reverse_D_module< D_module > reverse_view(D_module& d_module) {
  return Reverse_D_module< D_module >(d_module);
}

template< class D_module >
const Reverse_D_module< D_module > reverse_view(const D_module& d_module) {
  return Reverse_D_module< D_module >(d_module);
}

#endif  // TEST_D_MODULE_REVERSE_VIEW_