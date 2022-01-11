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
  // Template the Knot diagram class with the Morse events we allow in CSV files
  using Knot_diagram = Knot_diagram<
    Positive_crossing,
    Negative_crossing,
    Local_maximum,
    Local_minimum,
    Global_minimum
  >;
  
  // Program takes filename as argument
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
  
  // Import CSV file to knot diagram
  std::ifstream in_file(argv[1]);
  Knot_diagram knot_diagram;
  knot_diagram.import_csv(in_file);
  
  // Compute knot Floer homology as a Poincaré polynomial, choosing the
  // underlying D-module class based on the number of strands in the diagram:
  // this is done because smaller knot diagrams can represent idempotents using
  // integers.
  Poincare_polynomial pp;
  if (knot_diagram.max_n_strands() <= 30) {
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
#ifdef BUNDLED_HFK_DRAW_
  knot_diagram.TeXify(knot_diagram_out);
#endif  // BUNDLED_HFK_DRAW_
  polynomial_out << pp << std::flush;
  
  // Close files
  knot_diagram_out.close();
  polynomial_out.close();
  in_file.close();
  
  return 0;
}