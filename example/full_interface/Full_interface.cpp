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

//#define BUNDLED_HFK_DRAW_  // define for LaTeX-related functionality
#define BUNDLED_HFK_VERBOSE_  // define for more verbose console

#include <ctime>
#include <iostream>
#include <fstream>
#include <limits>  // numeric_limits
#include <string>

#include "Knot_interface.h"

void load_and_compute_regina(
  std::ifstream& knots_csv,
  std::ofstream& knot_diagram_out,
  std::ofstream& planar_diagram_out,
  std::ofstream& polynomial_out,
  int start = 1,
  int stop = -1
) {
  --stop;
  while (--start) {
    --stop;
    knots_csv.ignore(std::numeric_limits< std::streamsize >::max(), '\n');
  }
  for (std::string line; std::getline(knots_csv, line, '\n');) {
    int i1 = line.find(',');
    int i2 = line.find(',', i1 + 1);
    std::string name = line.substr(0, i1);
    std::string knot_sig = line.substr(i1 + 1, i2 - i1 - 1);
    
    std::cout << "[main] Loading " << name << "..." << std::endl;
    Knot_interface interface;
    interface.import_regina_signature(knot_sig);
    
#ifdef BUNDLED_HFK_DRAW_
    interface.knot_diagram().TeXify(knot_diagram_out);
    knot_diagram_out << std::endl;
#endif  // BUNDLED_HFK_DRAW_
#ifdef BUNDLED_HFK_VERBOSE_
    planar_diagram_out << interface.planar_diagram_string() << std::endl;
#endif  // BUNDLED_HFK_VERBOSE_
    
    auto pp = interface.knot_Floer_homology();
    std::cout << u8"[main] Poincar\u00E9 polynomial: " << pp << std::endl;
    
    polynomial_out << pp << std::endl;
    if (!--stop) { return; }
  }
}

void print_help() {
  std::cout
    << "Bundled HFK options:\n"
    << "  --gnu [w] [c]\n"
    << "  -h, --help\n"
    << "  -me, --morse-events <csv file>\n"
    << "  -pd, --planar-diagram <txt file> [<start index>] [<end index>]\n"
    << "  -rs, --regina-signature <csv file> [<start index>] [<end index>]\n";
}


void print_license(std::string option = "") {
  if (option == "w") {
    std::cout
      << "This program is distributed in the hope that it will be useful, "
      << "but WITHOUT ANY WARRANTY; without even the implied warranty of "
      << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the "
      << "GNU General Public License for more details.\n";
  }
  else if (option == "c") {
    std::cout
      << "This program is free software: you can redistribute it and/or modify "
      << "it under the terms of the GNU General Public License as published by "
      << "the Free Software Foundation, either version 3 of the License, or "
      << "(at your option) any later version.\n";
  }
}


int main(int argc, char* argv[]) {
  std::cout
    << "\nBundled HFK - a bordered Floer knot homology calculator\n\n"
    << "Copyright (C) 2021-2022  Isaac Ren\n\n"
    << "This program comes with ABSOLUTELY NO WARRANTY; for details run the "
    << "program with option `-gnu w'.\n"
    << "This is free software, and you are welcome to redistribute it under "
    << "certain conditions; run the program with option `-gnu c' for details.\n"
    << std::endl;
  
  if (argc == 1 or std::string(argv[1]) == "-h" or std::string(argv[1]) == "--help") {
    print_help();
    return 0;
  }
  
  if (std::string(argv[1]) == "-gnu") {
    if (argc == 2) {
      print_license();
    }
    else {
      print_license(std::string(argv[2]));
    }
    return 0;
  }
  
  if (argc == 2) {
    print_help();
    return 0;
  }
  
  std::clock_t c_start = std::clock();
  std::ifstream in_file(argv[2]);
  
  if (std::string(argv[1]) == "--morse-event" or std::string(argv[1]) == "-me") {
    std::cout << "[main] Reading morse events from " << argv[2] << "..." << std::endl;
    Knot_interface interface;
    interface.import_morse_events(argv[2]);
    auto pp = interface.knot_Floer_homology();
    std::cout << u8"[main] Poincar\u00E9 polynomial: " << pp << std::endl;
  }
  
  else if (std::string(argv[1]) == "--planar-diagram" or std::string(argv[1]) == "-pd") {
    std::cout << "[main] Reading " << argv[2] << "..." << std::endl;
    Knot_interface interface;
    
    /* Prepare list of planar diagram strings.
     * 
     * This part is taken from ComputeHFKv2/Main.cpp
     */
    std::vector< std::string > pd_strings;
    std::string planar_diagram_string;
    char ch;
    
    // Ignore everything until first planar diagram
    while(in_file >> ch) {
      if(ch == 'P') {
        in_file >> ch;
      }
      if(ch == 'D') {
        break;
      }
	}
	
    while (in_file >> ch) {
      if (ch != ' ' and ch != 'P') {
        planar_diagram_string.push_back(ch);
      }
      else if (ch == 'P') {
        pd_strings.push_back(planar_diagram_string);
        planar_diagram_string.clear();
      }
    }
    pd_strings.push_back(planar_diagram_string);
    
    int start = 0;
    int stop = pd_strings.size();
    
    if (argc >= 4) {
      start = std::stoi(argv[3]);
    }
    if (argc >= 5) {
      stop = std::stoi(argv[4]);
    }
    
    std::cout
      << "[main] Reading planar diagrams #"
      << (start + 1)
      << " to #"
      << stop
      << "...\n";
    
    for (int i = start; i < stop; ++i) {
      std::cout << "[main] Reading planar diagram #" << (i + 1) << "..." << std::endl;
      interface.import_planar_diagram(pd_strings[i]);
      auto pp = interface.knot_Floer_homology();
      std::cout << u8"[main] Poincar\u00E9 polynomial: " << pp << std::endl;
    }
  }
  
  else if (std::string(argv[1]) == "--regina-signature" or std::string(argv[1]) == "-rs") {
    std::ofstream knot_diagram_out("knot_diagrams.tex");
    std::ofstream planar_diagram_out("planar_diagrams.txt");
    std::ofstream polynomial_out("poincare_polynomials.tex");
    
    /* A knot number is provided. */
    if (argc >= 4) {
      int start = std::stoi(argv[3]);
      int stop = start;
      if (argc >= 5) {
        stop = std::stoi(argv[4]);
      }
      std::cout
        << "[main] Reading Regina signatures from "
        << argv[2]
        << " from "
        << start
        << " to "
        << stop
        << "..."
        << std::endl;
      
      in_file.ignore(std::numeric_limits< std::streamsize >::max(), '\n');
      load_and_compute_regina(in_file, knot_diagram_out, planar_diagram_out, polynomial_out, start, stop);
    }
    else {
      std::cout << "[main] Reading Regina signatures from " << argv[2] << "..." << std::endl;
      /* get rid of first line */
      in_file.ignore(std::numeric_limits< std::streamsize >::max(), '\n');
      load_and_compute_regina(in_file, knot_diagram_out, planar_diagram_out, polynomial_out);
    }
    
    knot_diagram_out.close();
    planar_diagram_out.close();
    polynomial_out.close();
  }
  
  in_file.close();
  std::clock_t c_end = std::clock();
  std::cout << "CPU time used: "
    << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC
    << " ms"
    << std::endl;
  
  return 0;
}