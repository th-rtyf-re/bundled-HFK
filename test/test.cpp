#define DRAW  // For TeXify

#include <iostream>
#include <fstream>

#include "test_forest.h"
#include "test_reverse_d_module.h"
#include "test_knot_diagram.h"

#include "test_local_maximum.h"
#include "test_local_minimum.h"
#include "test_global_minimum.h"

#include "Math_tools/Poincare_polynomial.h"

int main() {
  using D_module_short = Forest<>;
  using D_module_long = Forest< Forest_options_default_long >;
  
  using KD = Knot_diagram< Local_maximum, Local_minimum, Global_minimum >;
  
  KD kd;
  
  kd.import_csv("test_morse_events.csv");
    
  std::cout << "Poincare polynomial "
            << kd.knot_Floer_homology< Poincare_polynomial, D_module_short >() << std::endl;
  
  
  std::ofstream knot_diagram("knot_diagram.tex");
  kd.TeXify(knot_diagram);
  knot_diagram.close();
  return 0;
}