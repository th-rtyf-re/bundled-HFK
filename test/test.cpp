#define DRAW  // define for LaTeX-related functionality
#define VERBOSE  // define for more verbose console

#include <iostream>
#include <fstream>

#include "test_forest.h"
#include "test_knot_diagram.h"

#include "test_positive_crossing.h"
#include "test_negative_crossing.h"
#include "test_local_maximum.h"
#include "test_local_minimum.h"
#include "test_global_minimum.h"

#include "Math_tools/Poincare_polynomial.h"

int main() {
  using D_module_short = Forest<>;
  using D_module_long = Forest< Forest_options_default_long >;
  
  using KD = Knot_diagram< Positive_crossing, Negative_crossing, Local_maximum, Local_minimum, Global_minimum >;
  
  KD kd;
  
  kd.import_csv("test_morse_events.csv");
  
  std::ofstream knot_diagram("knot_diagram.tex");
  kd.TeXify(knot_diagram);
  knot_diagram.close();
  
  auto pp = kd.knot_Floer_homology< Poincare_polynomial, D_module_short >();
  std::cout << "Poincare polynomial: " << pp << std::endl;
  
  return 0;
}


// struct B {
//   template<class T> void f() { }
// };
// 
// template< class B >
// struct D : B {
//   using B::f;
//   void g() { f<int>(); }
// };
// 
// int main() {
//   D< B > d;
//   d.g();
//   return 0;
// }