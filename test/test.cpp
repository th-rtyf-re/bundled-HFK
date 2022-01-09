//#define DRAW  // define for LaTeX-related functionality
#define VERBOSE  // define for more verbose console

#include <ctime>
#include <iostream>
#include <fstream>
#include <limits>  // numeric_limits
#include <string>

#include "Knot_interface.h"

void aux(std::ifstream& knots_csv,
         std::ofstream& knot_diagram_out,
         std::ofstream& planar_diagram_out,
         std::ofstream& polynomial_out,
         int start = 1,
         int stop = -1) {
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
    
#ifdef DRAW
    interface.knot_diagram().TeXify(knot_diagram_out);
    knot_diagram_out << std::endl;
#endif  // DRAW
#ifdef VERBOSE
    planar_diagram_out << interface.planar_diagram_string() << std::endl;
#endif  // VERBOSE
    
    auto pp = interface.knot_Floer_homology();
    std::cout << u8"[main] Poincar\u00E9 polynomial: " << pp << std::endl;
    
    polynomial_out << pp << std::endl;
    if (!--stop) { return; }
  }
}

int main(int argc, char* argv[]) {
  if (argc <= 2) {
    std::cout
      << "arguments: "
      << "-rs/--regina-signature [csv file] [knot number/start] [knot end], "
      << "-pd/--planar-diagram [txt file], "
      << "-me/--morse-events [csv file]."
      << std::endl;
    return 0;
  }
  
  std::clock_t c_start = std::clock();
  std::ifstream in_file(argv[2]);
  
  if (std::string(argv[1]) == "--planar-diagram" or std::string(argv[1]) == "-pd") {
    std::cout << "[main] Reading planar diagrams from " << argv[2] << "..." << std::endl;
    Knot_interface interface;
    
    std::vector< std::string > pd_strings;
    std::string planar_diagram_string;
    char ch;
    
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
    
    interface.import_planar_diagram(pd_strings[0]);
    auto pp = interface.knot_Floer_homology();
    std::cout << u8"[main] Poincar\u00E9 polynomial: " << pp << std::endl;
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
      aux(in_file, knot_diagram_out, planar_diagram_out, polynomial_out, start, stop);
    }
    else {
      std::cout << "[main] Reading Regina signatures from " << argv[2] << "..." << std::endl;
      /* get rid of first line */
      in_file.ignore(std::numeric_limits< std::streamsize >::max(), '\n');
      aux(in_file, knot_diagram_out, planar_diagram_out, polynomial_out);
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