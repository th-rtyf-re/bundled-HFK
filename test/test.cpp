#define DRAW  // For TeXify

#include <iostream>
#include <fstream>

#include "test_forest.h"
#include "test_local_maximum.h"
#include "test_knot_diagram.h"

#include "Math_tools/Poincare_polynomial.h"

int main() {
  using D_module_short = Forest<>;
  using D_module_long = Forest< Forest_options_default_long >;
  
  using KD = Knot_diagram< D_module_pack< D_module_short, D_module_long >, Morse_event_pack< Local_maximum > >;
  
  KD kd;
  
  kd.import_csv("test_morse_events.csv");
  
  std::cout << "[main] Doing math..." << std::endl;
  
  std::cout << "Poincare polynomial "
            << kd.knot_Floer_homology< Poincare_polynomial, D_module_short >() << std::endl;
  
  

  return 0;
}


// int foo() {
//   using Forest = Forest<>;
//   using Idem = typename Forest::Idem;
//   using Alg_el = typename Forest::Alg_el;
//   
//   std::vector< std::pair< int, int > > monomial_defaults(10, {0, 0});
//   std::vector< std::string > label_defaults(10, "X");
//   
//   Forest f0;
//   
//   f0.add_gen_bundle(Idem("0"));
//   f0.add_gen_bundle(Idem("1"));
//   f0.lock_generators();
//   
//   std::cout << f0;
//   
//   Forest f1;
//   f1.add_gen_bundle(Idem("0"), 1, 0);
//   f1.add_gen_bundle(Idem("0"), 2, 1);
//   f1.add_gen_bundle(Idem("1"), 3, 0);
//   f1.add_gen_bundle(Idem("1"), 4, 1);
//   f1.lock_generators(f0, monomial_defaults, label_defaults);
//   
//   f1.add_coef_bundle(Alg_el(Idem("0"), Idem("1")), 1, 3, Idem("0"));
//   f1.lock_coefficients();
//   
//   std::cout << f1;
//   
//   Forest f2;
//   f2.add_gen_bundle(Idem("0"), 5, 0);
//   f2.add_gen_bundle(Idem("0"), 6, 3);
//   f2.add_gen_bundle(Idem("1"), 7, 0);
//   f2.add_gen_bundle(Idem("1"), 8, 3);
//   f2.lock_generators(f1, monomial_defaults, label_defaults);
//   
//   Alg_el a1(Idem("0"), Idem("0"));
//   Alg_el a2(Idem("1"), Idem("1"));
//   
//   f2.add_coef_bundle(0, 7, a1);
//   f2.add_coef_bundle(3, 10, a1);
//   f2.add_coef_bundle(3, 2, a2);
//   f2.add_coef_bundle(11, 1, a2);
//   f2.add_coef_bundle(4, 8, a1);
//   
//   f2.add_coef_bundle(Alg_el(Idem("0"), Idem("1")), 5, 8, *f1.coef_bundles().begin(), f1);
//   f2.lock_coefficients();
//   
//   std::cout << f2;
//   f2.reduce();
//   
//   std::cout << f2;
//     
// #ifdef DRAW
//   std::ofstream suffix_forest("differential_suffix_forest.tex");
//   f1.TeXify(suffix_forest);
//   f2.TeXify(suffix_forest);
//   suffix_forest.close();
// #endif  // DRAW
//   
//   return 0;
// }