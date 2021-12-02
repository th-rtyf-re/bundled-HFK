/* The user defines classes that satisfy the Morse event concept, but it's
 * easier for Knot diagram to reason with DA bimodules.
 * These are constructed when we call the knot Floer homology function.
 */
template< class Morse_event, class D_module >
struct DA_bimodule {
  using Algebra = typename D_module::Bordered_algebra;
  
  DA_bimodule(const Morse_event morse_event,
              const Algebra upper_algebra,
              const Algebra lower_algebra) :
    morse_event(morse_event),
    upper_algebra(upper_algebra),
    lower_algebra(lower_algebra)
  { }
  
  const Morse_event morse_event;
  const Algebra upper_algebra;
  const Algebra lower_algebra;
  
  friend D_module box_tensor_product(const DA_bimodule& da_bimodule,
                                     const D_module& old_d_module) {
    D_module new_d_module =
      da_bimodule.morse_event.tensor_generators(old_d_module,
                                                da_bimodule.upper_algebra,
                                                da_bimodule.lower_algebra);
    std::cout << "[dab] tensored generators" << std::endl;
    
    auto weights = da_bimodule.morse_event.weights(da_bimodule.upper_algebra,
                                                   da_bimodule.lower_algebra);
    std::cout << "[dab] computed weights" << std::endl;
    auto labels = da_bimodule.morse_event.labels(da_bimodule.upper_algebra,
                                                 da_bimodule.lower_algebra);
    std::cout << "[dab] computed labels" << std::endl;
    new_d_module.lock_generators(old_d_module, weights, labels);
    std::cout << "[dab] locked generators" << std::endl;
    da_bimodule.morse_event.tensor_coefficients(old_d_module,
                                                new_d_module,
                                                da_bimodule.upper_algebra,
                                                da_bimodule.lower_algebra);
    std::cout << "[dab] tensored coefficients" << std::endl;
    new_d_module.lock_coefficients();
    return new_d_module;
  }
};