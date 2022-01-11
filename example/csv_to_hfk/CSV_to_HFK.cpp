#define BUNDLED_HFK_DRAW_  // define for LaTeX-related functionality
//#define BUNDLED_HFK_VERBOSE_  // define for more verbose console

#include <iostream>

#include "Differential_suffix_forest/Differential_suffix_forest.h"
#include "Differential_suffix_forest/Differential_suffix_forest_options.h"
#include "Knot_diagram/Knot_diagram.h"
#include "Math_tools/Poincare_polynomial.h"
#include "Morse_event/Positive_crossing.h"
#include "Morse_event/Negative_crossing.h"
#include "Morse_event/Local_maximum.h"
#include "Morse_event/Local_minimum.h"
#include "Morse_event/Global_minimum.h"

int main(int argc, char* argv[]) {
  using Knot_diagram = Knot_diagram<
    Positive_crossing,
    Negative_crossing,
    Local_maximum,
    Local_minimum,
    Global_minimum
  >;
  
  if (argc <= 1) {
    std::cout << "[main] No file given! Exiting..." << std::endl;
    return 0;
  }
  
  std::cout << "[main] Reading CSV file " << argv[1] << "..." << std::endl;
  
  // Set log file, where developer messages are sent.
  std::ofstream log_file("log.txt");
  std::clog.rdbuf(log_file.rdbuf());
  
  // Other output files
  std::ofstream knot_diagram_out("knot_diagrams.tex");
  std::ofstream polynomial_out("poincare_polynomials.tex");
  
  std::ifstream in_file(argv[1]);
  Knot_diagram knot_diagram;
  knot_diagram.import_csv(in_file);
  
  Poincare_polynomial pp;
  
  // Choose D-module class based on number of strands
  if (knot_diagram.max_n_strands() <= 31) {
    pp = knot_diagram.knot_Floer_homology<
      Poincare_polynomial,
      Differential_suffix_forest< Forest_options_default_short >
    >();
  }
  else {
    pp = knot_diagram.knot_Floer_homology<
      Poincare_polynomial,
      Differential_suffix_forest< Forest_options_default_long >
    >();
  }
  
  // Output knot Floer homology and other information
  std::cout << u8"[main] Poincar\u00E9 polynomial: " << pp << std::endl;
  knot_diagram.TeXify(knot_diagram_out);
  polynomial_out << pp << std::flush;
  
  // Close files
  knot_diagram_out.close();
  polynomial_out.close();
  in_file.close();
  
  return 0;
}