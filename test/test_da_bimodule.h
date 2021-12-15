#ifdef VERBOSE
#include <iostream>
#endif  // VERBOSE

/* The user defines classes that satisfy the Morse event concept, but it's
 * easier for Knot diagram to reason with DA bimodules.
 * These are constructed when we call the knot Floer homology function.
 */
template< class Morse_event, class D_module >
class DA_bimodule {
 public:
  using Algebra = typename D_module::Bordered_algebra;
  
  DA_bimodule(
    const Morse_event morse_event,
    const Algebra upper_algebra,
    const Algebra lower_algebra
  ) :
    morse_event(morse_event),
    upper_algebra(upper_algebra),
    lower_algebra(lower_algebra)
  { }
  
  const Morse_event morse_event;
  const Algebra upper_algebra;
  const Algebra lower_algebra;
  
  friend D_module box_tensor_product(
    const DA_bimodule& da_bimodule,
    const D_module& old_d_module
  ) {
    D_module new_d_module;
    da_bimodule.morse_event.tensor_generators(
      new_d_module,
      old_d_module,
      da_bimodule.upper_algebra,
      da_bimodule.lower_algebra
    );
    const auto weights = da_bimodule.morse_event.get_weights(
      da_bimodule.upper_algebra,
      da_bimodule.lower_algebra
    );
    const auto labels = da_bimodule.morse_event.get_labels(
      da_bimodule.upper_algebra,
      da_bimodule.lower_algebra
    );
    new_d_module.lock_generators(old_d_module, weights, labels);
    da_bimodule.morse_event.tensor_coefficients(
      new_d_module,
      old_d_module,
      da_bimodule.upper_algebra,
      da_bimodule.lower_algebra
    );
    new_d_module.lock_coefficients();
    return new_d_module;
  }
  
#ifdef VERBOSE
  friend std::ostream& operator<<(
    std::ostream& os,
    const DA_bimodule& da_bimodule
  ) {
    os << da_bimodule.morse_event;
    return os;
  }
#endif  // VERBOSE
};