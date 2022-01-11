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

#ifndef KNOT_DIAGRAM_H_
#define KNOT_DIAGRAM_H_

#include <fstream>
#include <iostream>
#include <string>  // stoi
#include <utility>  // pair
#include <vector>

//#include <boost/any.hpp>

#include "Differential_suffix_forest/Differential_suffix_forest.h"
#include "Math_tools/DA_bimodule.h"
#include "Morse_event/Morse_event.h"
#include "Morse_event/Morse_event_options.h"

/* Knot diagram.
 * 
 * This class stores the information of a knot through a knot diagram, and can
 * compute its knot Floer homology.
 * 
 * The knot diagram is a list of Morse events, which are determined by the user
 * via the class template.
 */
template< template< class, class > class ...Morse_events >
class Knot_diagram {
 public:
  using D_module_default = Differential_suffix_forest<>;
  using Morse_event_options = Morse_event_options_int;
  using Parameter_type = typename Morse_event_options::Parameter_type;
  
  /* A way to import Morse events. The integer is the index of the Morse event
   * in the template pack, and the vector of parameters are the Morse event
   * constructor arguments.
   */
  using Morse_data_container =
    std::vector< std::pair< int, std::vector< Parameter_type > > >;
  
  Knot_diagram () { }
  
  void import_data(const Morse_data_container& morse_data) {
    morse_data_ = morse_data;
    morse_events_default_ = Detail_<>::get_morse_events(morse_data_);
  }
  
  /* Import from CSV.
   * Import a Morse event list from a CSV file. Currently, we assume that the
   * lines of the file are of the form
   * 
   *      event,position.
   * 
   */
  void import_csv(std::ifstream& morse_event_csv) {
    std::cout << "[kd] Importing Morse event list..." << std::endl;
    
    if (!morse_event_csv.is_open()) {  // file not found
      std::cout << "[kd] Morse event CSV file not found!"
                << " Exiting..."
                << std::endl;
      exit(1);
    }
    
    int line_number = 1;
    for (
      std::string line;
      std::getline(morse_event_csv, line, '\n');
      ++line_number
    ) {
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
      int event = std::stoi(line.substr(0, i));
      int position = std::stoi(line.substr(i + 1));
      morse_data_.push_back({event, {position}});
    }
    
    morse_events_default_ = Detail_<>::get_morse_events(morse_data_);
  }  // import_csv
  
  int max_n_strands() const {
    int n_strands = 0;
    int max_n_strands = 0;
    
    for (const auto morse_event : morse_events_default_) {
      // hack to get lower n strands via the lower matchings vector
      n_strands =
        morse_event.lower_matchings(std::vector< int >(n_strands, 0)).size();
      if (n_strands > max_n_strands) {
        max_n_strands = n_strands;
      }
    }
    return max_n_strands;
  }
  
  template< class Polynomial, class D_module >
  Polynomial knot_Floer_homology() const {
#ifdef BUNDLED_HFK_VERBOSE_
    std::cout << "[kd] Computing knot Floer homology..." << std::endl;
#endif  // BUNDLED_HFK_VERBOSE_
#ifdef BUNDLED_HFK_DRAW_
    std::ofstream suffix_forest("differential_suffix_forest.tex");
#endif  // BUNDLED_HFK_DRAW_
    
    const auto da_bimodules = Detail_< D_module >::get_da_bimodules(morse_data_);
    
    D_module d_module;
    d_module.set_as_trivial();
    
    // box tensor product for each Morse event
    for (int i = 0; i != da_bimodules.size(); ++i) {
#ifdef BUNDLED_HFK_VERBOSE_
      std::cout << "[kd] layer " << i << ": " << da_bimodules[i] << "... " << std::flush;
#endif  // BUNDLED_HFK_VERBOSE_
      d_module = box_tensor_product(da_bimodules[i], d_module);
#ifdef BUNDLED_HFK_DRAW_
      suffix_forest << "Before reduction:\n";
      d_module.TeXify(suffix_forest);
      suffix_forest << "\n" << std::flush;
#endif  // BUNDLED_HFK_DRAW_
#ifdef BUNDLED_HFK_VERBOSE_
      std::cout << "reducing... " << std::flush;
#endif  // BUNDLED_HFK_VERBOSE_
      d_module.reduce();
#ifdef BUNDLED_HFK_DRAW_
      suffix_forest << "After reduction:\n";
      d_module.TeXify(suffix_forest);
      suffix_forest << "\n\n" << std::flush;
#endif  // BUNDLED_HFK_DRAW_
#ifdef BUNDLED_HFK_VERBOSE_
      std::cout << "done." << std::endl;
#endif  // BUNDLED_HFK_VERBOSE_
    }
    
#ifdef BUNDLED_HFK_DRAW_
    suffix_forest.close();
#endif  // BUNDLED_HFK_DRAW_
    
    return d_module.template poincare_polynomial< Polynomial >();
  }
  
 private:
  /* All the private methods depend on a choice of D-module. In order to
   * make things hopefully more readable, I've put all these methods as static
   * methods in a private struct Detail_, templated by D_module.
   */
  template< class D_module = D_module_default >
  struct Detail_ {
    using Bordered_algebra = typename D_module::Bordered_algebra;
    using Morse_event = Morse_event< D_module, Morse_event_options >;
    using DA_bimodule = DA_bimodule< Morse_event >;
    
    static std::vector< DA_bimodule > get_da_bimodules(
      const Morse_data_container& morse_data
    ) {
      std::vector< DA_bimodule > da_bimodules;
      
      const auto morse_events = get_morse_events(morse_data);
      const auto algebras = get_bordered_algebras(morse_events);
      
      for (int i = 0; i != morse_events.size(); ++i) {
        da_bimodules.emplace_back(
          morse_events[i],
          algebras[i],
          algebras[i + 1]
        );
      }
      
      return da_bimodules;
    }
    
    static std::vector< Morse_event > get_morse_events(
      const Morse_data_container& morse_data
    ) {
      std::vector< Morse_event > morse_events;
      
      for (const auto value_pair : morse_data) {
        morse_events.push_back(instance(value_pair.first, value_pair.second));
      }
      
      return morse_events;
    }
    
    static std::vector< Bordered_algebra > get_bordered_algebras(
      const std::vector< Morse_event >& morse_events
    ) {
      std::vector< Bordered_algebra > algebras(morse_events.size() + 1);
      
      // Place algebras
      algebras[0].n_strands = 0;
      
      // Calculate matchings and n_strands top-down
      for (int i = 0; i < morse_events.size() - 1; ++i) {
        algebras[i + 1].matchings =
          morse_events[i].lower_matchings(algebras[i].matchings);
        algebras[i + 1].n_strands = algebras[i + 1].matchings.size();
      }
      
      // Calculate orientations bottom-up
      for (int i = morse_events.size() - 1; i >= 0; --i) {
        algebras[i].orientations = morse_events[i].upper_orientations(
          algebras[i + 1].orientations,
          algebras[i].matchings
        );
      }
      return algebras;
    }
    
    /* Return an instance of the i^th Morse event, as listed in the template,
     * constructed using arguments.
     */
    static Morse_event instance(const int i, const std::vector< Parameter_type >& args) {
      return instance_aux_< 0, Morse_events... >(i, args);
    }
    
   private:
    /* Template-recursive function creating an instance of the i^th Morse event
     * 
     * I need to add a dummy class template because I can't partially
     * specialize for the initialization otherwise.
     * The other possibility is requiring at least one Morse event, but
     * I want to leave the possibility of having no Morse events.
     */
    template<
      int,
      template< class, class > class Morse_events_head,
      template< class, class > class ...Morse_events_tail >
    static Morse_event instance_aux_(
      const int i,
      const std::vector< Parameter_type >& args
    ) {
      if (i == 0) {
        return Morse_event(Morse_events_head< D_module, Morse_event_options >(args));
      }
      else {
        return instance_aux_< 0, Morse_events_tail... >(i - 1, args);
      }
    }
    
    // Initialization: 
    template< int >
    static Morse_event instance_aux_(const int, const std::vector< Parameter_type >&) {
      return Morse_event();
    }
  };  // Detail_
  
 public:
#ifdef BUNDLED_HFK_DRAW_
  /* TeXify */
  void TeXify(std::ofstream& write_file) const {
    std::cout << "[kd] Drawing knot..." << std::flush;
    
    const auto margins = get_margins_();
    const auto algebras = Detail_<>::get_bordered_algebras(morse_events_default_);
    
    write_file << "\\KnotDiagram{" << std::endl;
    for (int i = 0; i < morse_events_default_.size(); ++i) {
      write_file << morse_events_default_[i].to_string(
        margins[i],
        {algebras[i].n_strands, algebras[i + 1].n_strands}
      ) << std::endl;
	}
    write_file << "}" << std::flush;
    
    std::cout << "done." << std::endl;
  }
  
 private:
  std::vector< std::pair< int, int > > get_margins_() const {
    std::vector< std::pair< int, int > > margins;
    
    const int m = max_n_strands();
    
    std::pair< int, int > current_margins = {m / 2, m / 2};
    for (const auto& morse_event : morse_events_default_) {
      current_margins = morse_event.update_margins(current_margins);
      margins.push_back(current_margins);
    }
    return margins;
  }
#endif  // BUNDLED_HFK_DRAW_
  
  Morse_data_container morse_data_;
  
  std::vector< Morse_event< D_module_default, Morse_event_options > > morse_events_default_;
};

#endif  // KNOT_DIAGRAM_H_