/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *  Bundled HFK - a knot Floer homology calculator                           *
 *                                                                           *
 *  Copyright (C) 2021  Isaac Ren                                            *
 *  For further details contact Isaac Ren (gopi3.1415@gmail.com)             *
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

#define DRAW
//#define DEBUG

#include <iostream>
#include <fstream>
#include <vector>

#include "Differential_suffix_forest/Differential_suffix_forest.h"
#include "Differential_suffix_forest/Differential_suffix_forest_options.h"
#include "Knot_diagram/Knot_diagram.h"

/* CSV to HFK.
 * 
 * Take as input a .csv file with
 * 
 *     position, event
 * 
 * on each line, and return the knot Floer homology as a Poincaré polynomial.
 * 
 * Parameters:
 * - define DRAW to draw knot diagram and differential suffix trees to .tex
 *   files, which can then be compiled in drawings.tex
 * - define DEBUG to print a lot of undocumented information in console.
 */

template< class Forest_options >
void aux(const Knot_diagram& knot_diagram) {
  using DSF = Differential_suffix_forest< Forest_options >;
  auto pp = knot_diagram.knot_Floer_homology< DSF >();
  std::cout << u8"[hfk] Poincar\u00E9 polynomial: " << pp.to_string() << std::endl;
}

/* First argument is the Morse event CSV file */
int main(int argc, char* argv[]) {
  Knot_diagram knot_diagram;
  if (argc == 1) {
    std::cout << "[hfk] No file given! Exiting..." << std::endl;
    return 1;
  }
  else if (argc >= 2) {
    if (argc > 2) {
      std::cout << "[hfk] Extra options detected. Ignoring them..." << std::endl;
    }
    knot_diagram.import_csv(argv[1]);
  }
  
#ifdef DRAW
  std::ofstream knot_diagram_out("knot_diagram.tex");
  knot_diagram.TeXify(knot_diagram_out);
  knot_diagram_out.close();
#endif  // DRAW
  
  /* Check if we take short or long idempotents */
  if (knot_diagram.max_n_strands() <= 31) {
    aux< Forest_options_default_short >(knot_diagram);
  }
  else {
    aux< Forest_options_default_long >(knot_diagram);
  }
  
  std::cout << "[test] done" << std::endl;
  return 0;
}