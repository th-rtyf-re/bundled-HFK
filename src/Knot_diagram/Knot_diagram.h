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

#ifndef KNOT_DIAGRAM_H_
#define KNOT_DIAGRAM_H_

#include <fstream>
#include <string>  // stoi
#include <utility>  // swap
#include <vector>

#include "Morse_event/Event_registers.h"  // Knot_slice_register, DA_bimodule_register

/* Knot diagram.
 * 
 * Class that serves as an interface with the input knot. Currently, we take
 * a CSV file of Morse events. The CSV file should have two integers per line,
 * separated by a comma:
 *
 *   position0,event0
 *   position1,event1
 *   etc.
 *
 * The code for events is:
 * - 0: positive crossing
 * - 1: negative crossing
 * - 2: local maximum
 * - 3: local minimum
 * Comment out a line by adding a '#' at the beginning.
 * 
 * This class also calculates knot Floer homology, provided a class satisfying
 * some "D-module" concept.
 */
class Knot_diagram {
 public:
  using Knot_slice_handle = typename Knot_slice_register::Knot_slice_handle;
  
  Knot_diagram() { }
  
  /* Delete copy constructor and assignment because of unique pointers */
  Knot_diagram(Knot_diagram&) = delete;
  Knot_diagram& operator=(Knot_diagram&) = delete;
  Knot_diagram(Knot_diagram&&) = default;
  Knot_diagram& operator=(Knot_diagram&&) = default;
  ~Knot_diagram() = default;
  
  /* Import from CSV.
   * Import a Morse event list from a CSV file given by filename.
   */
  template< class Filename_type >
  void import_csv(const Filename_type& filename) {
    std::cout << "[kd] Importing morse event list..." << std::endl;
    std::ifstream morse_event_csv(filename);
    
    if (!morse_event_csv.is_open()) {  // file not found
      std::cout << "[kd] Morse event CSV file not found!"
                << " Exiting..."
                << std::endl;
      exit(1);
    }
    
    std::vector< std::pair< int, Morse_event > > morse_events;
    int line_number = 1;
    for (std::string line; std::getline(morse_event_csv, line, '\n'); ++line_number) {
      if (line[0] == '#') { continue; }  // comment line
      else if (line.find_first_of("0123456789") != 0) {  // syntax error
        std::cout << "[kd] Line "
                  << line_number
                  << " is incorrectly formatted."
                  << " Commented lines start with '#'."
                  << " Ignoring line..."
                  << std::endl;
        continue;
      }
      auto i = line.find(',');
      std::string position = line.substr(0, i);
      std::string pre_event = line.substr(i + 1);
      morse_events.push_back({std::stoi(position), static_cast< Morse_event >(std::stoi(pre_event))});
    }
    morse_event_csv.close();
    
    import_morse_events(morse_events);
  }  // import_csv
  
  /* Create knot slices based on a list (position, event) */
  void import_morse_events(std::vector< std::pair< int, Morse_event > > morse_events) {
    if (!morse_events.empty()) {
      morse_events.back().second = Morse_event::global_minimum;  // set global minimum
    }
    
    /* Add knot slices and calculate max_n_strands_ at the same time */
    int n_strands = 0;
    max_n_strands_ = 0;
    for (auto& value : morse_events) {
      int position = value.first;
      auto event_id = value.second;
      knot_slices_.push_back(Knot_slice_register::knot_slice(event_id, position, n_strands));
      n_strands = knot_slices_.back()->lower_n_strands();
      if (n_strands > max_n_strands_) {
        max_n_strands_ = n_strands;
      }
    }
    
    if (n_strands > 0) {
      std::cout << "[kd] Knot diagram not closed. Behavior is undefined." << std::endl;
    }
  }  // import_morse_events
  
  int max_n_strands() const {
    return max_n_strands_;
  }
  
  /* Calculate knot Floer homology using the templated classes for D-modules
   * and DA-bimodules.
   * 
   * D_module must satisfy some "D-module" concept.
   */
  template< class D_module >
  typename D_module::Polynomial knot_Floer_homology() const {
    auto da_bimodules = da_bimodules_< D_module >();
    
    std::ofstream write_file("differential_suffix_tree.tex");  // debug LaTeX file
    D_module d_module;
    d_module.empty_diagram();  // TO DO: move to a Knot_slice?
    
    /* Box tensor product for each slice */
    for (auto da_bimodule_it = da_bimodules.begin();
              da_bimodule_it != da_bimodules.end(); ++da_bimodule_it) {
#ifdef DEBUG
      std::cout << "[hfk] Tensoring... " << std::flush;
#endif  // DEBUG
      d_module = box_tensor_product(d_module, **da_bimodule_it);
      
#ifdef DEBUG
      std::cout << "Reducing... " << std::flush;
#endif  // DEBUG
      while (!d_module.reduce(da_bimodule_it, write_file)) { /* Not done reducing */ }
      
#ifdef DEBUG
      d_module.debug();
      std::cout << "Number of arcs: " << d_module.aggl_coefs().size() << std::flush;
      std::cout << std::endl;
#endif  // DEBUG
    }
        
    write_file.close();
    return d_module.poincare_polynomial(std::prev(da_bimodules.end()));
  }  // knot_Floer_homology
  
 private:
  /* Calculate and return DA-bimodule handles */
  template< class D_module, class DA_bimodule_handle = typename DA_bimodule_register< D_module >::DA_bimodule_handle >
  std::vector< DA_bimodule_handle > da_bimodules_() const {
    auto algebras = bordered_algebras_< typename D_module::Bordered_algebra >();
    std::vector< DA_bimodule_handle > da_bimodules;
    for (int i = 0; i < knot_slices_.size(); ++i) {
      da_bimodules.push_back(DA_bimodule_register< D_module >::da_bimodule(knot_slices_[i], algebras[i], algebras[i + 1]));
    }
    return da_bimodules;
  }
  
  /* Calculate and return bordered algebras */
  template< class Bordered_algebra >
  std::vector< Bordered_algebra > bordered_algebras_() const {
    std::vector< Bordered_algebra > algebras;
    
    /* Place algebras */
    Bordered_algebra inital_algebra;
    inital_algebra.n_strands = 0;
    algebras.push_back(inital_algebra);
    
    for (const auto& slice : knot_slices_) {
      Bordered_algebra lower_algebra;
      lower_algebra.n_strands = slice->lower_n_strands();
      algebras.push_back(lower_algebra);
    }
    
    /* Calculate matchings top-down */
    for (int i = 0; i < knot_slices_.size() - 1; ++i) {
      algebras[i + 1].matchings = knot_slices_[i]->lower_matchings(algebras[i].matchings);
    }
    
    /* Calculate orientations bottom-up */
    for (int i = knot_slices_.size() - 1; i >= 0; --i) {
      algebras[i].orientations = knot_slices_[i]->upper_orientations(algebras[i + 1].orientations, algebras[i].matchings);
    }
    return algebras;
  }
  
 public:
  /* TeXify */
  void TeXify(std::ofstream& write_file) const {
    std::cout << "[kd] Drawing knot..." << std::endl;
    std::vector< std::pair< int, int > > margins = get_margins_();
    
    //std::ofstream write_file;
    //write_file.open(filename);
    
    write_file << "\\KnotDiagram{" << std::endl;
    for (int i = 0; i < knot_slices_.size(); ++i) {
      write_file << knot_slices_[i]->to_string(margins[i]) << std::endl;
	}
    write_file << "}" << std::flush;
    
    //write_file.close();
  }
  
 private:
  std::vector< std::pair< int, int > > get_margins_() const {
    std::vector< std::pair< int, int > > margins;
    std::pair< int, int > current_margins = {max_n_strands_ / 2, max_n_strands_ / 2};
    for (const auto& slice : knot_slices_) {
      current_margins = slice->update_margins(current_margins);
      margins.push_back(current_margins);
    }
    return margins;
  }
  
  std::vector< Knot_slice_handle > knot_slices_;
  int max_n_strands_;
};

#endif  // KNOT_DIAGRAM_H_