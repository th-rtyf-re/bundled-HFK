#include <fstream>
#include <iostream>
#include <string>  // stoi
#include <utility>  // pair
#include <vector>

#include <boost/any.hpp>

#include "test_da_bimodule.h"
#include "test_forest.h"
#include "test_morse_event.h"

/* None of this is optimized, because its cost is minimal in the big picture.
 */
template< template< class > class ...Morse_events >
class Knot_diagram {
 public:
  using D_module_default = Forest<>;
  
  Knot_diagram () { }
  
  /* Import from CSV.
   * Import a Morse event list from a CSV file given by filename.
   * 
   * event, position (different from current version)
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
      int event = std::stoi(line.substr(0, i));
      int position = std::stoi(line.substr(i + 1));
      morse_event_data_.push_back({event, {position}});
    }
    morse_event_csv.close();
    
    morse_events_default_ = Detail_<>::get_morse_events(morse_event_data_);
  }  // import_csv
  
  int max_n_strands() const {
    int n_strands = 0;
    int max_n_strands = 0;
    
    for (const auto morse_event : morse_events_default_) {
      // hack to get lower n strands via the lower matchings vector
      n_strands = morse_event.lower_matchings(std::vector< int >(n_strands, 0)).size();
      if (n_strands > max_n_strands) {
        max_n_strands = n_strands;
      }
    }
    return max_n_strands;
  }
  
  template< class Polynomial, class D_module >
  Polynomial knot_Floer_homology() const {
#ifdef VERBOSE
    std::cout << "[kd] Computing knot Floer homology..." << std::endl;
#endif  // VERBOSE
#ifdef DRAW
    std::ofstream suffix_forest("differential_suffix_forest.tex");
#endif  // DRAW
    
    const auto da_bimodules = Detail_< D_module >::get_da_bimodules(morse_event_data_);
    
    //std::ofstream write_file("differential_suffix_forest.tex");  // debug LaTeX file
    D_module d_module;
    d_module.set_as_trivial();
    
    // box tensor product for each Morse event
    for (const auto& da_bimodule : da_bimodules) {
      std::cout << "[kd] Tensoring... " << std::flush;
      d_module = box_tensor_product(da_bimodule, d_module);
#ifdef DRAW
      suffix_forest << "Before reduction:\n";
      d_module.TeXify(suffix_forest);
      suffix_forest << "\n" << std::flush;
#endif  // DRAW
      std::cout << "reducing... " << std::flush;
      d_module.reduce();
#ifdef DRAW
      suffix_forest << "After reduction:\n";
      d_module.TeXify(suffix_forest);
      suffix_forest << "\n\n" << std::flush;
#endif  // DRAW
      std::cout << "done." << std::endl;
    }
    
#ifdef DRAW
    suffix_forest.close();
#endif  // DRAW
    
    return d_module.template poincare_polynomial< Polynomial >();
  }
  
 private:
  
  /* All the private methods depend on the choice of D-module. In order to
   * make things hopefully more readable, I've put all these methods in a
   * private struct Detail_, templated by D_module, as static methods.
   */
  template< class D_module = D_module_default >
  struct Detail_ {
    using Bordered_algebra = typename D_module::Bordered_algebra;
    using Morse_event = Morse_event< D_module >;
    using DA_bimodule = DA_bimodule< Morse_event, D_module >;
    
    static std::vector< DA_bimodule > get_da_bimodules(const std::vector< std::pair< int, std::vector< boost::any > > >& morse_event_data) {
      std::vector< DA_bimodule > da_bimodules;
      
      const auto morse_events = get_morse_events(morse_event_data);
      const auto algebras = get_bordered_algebras(morse_events);
      
      for (int i = 0; i != morse_events.size(); ++i) {
        da_bimodules.emplace_back(morse_events[i],
                                  algebras[i],
                                  algebras[i + 1]);
      }
      
      return da_bimodules;
    }
    
    static std::vector< Morse_event > get_morse_events(const std::vector< std::pair< int, std::vector< boost::any > > >& morse_event_data) {
      std::vector< Morse_event > morse_events;
      
      for (const auto value_pair : morse_event_data) {
        morse_events.push_back(instance(value_pair.first, value_pair.second));
      }
      
      return morse_events;
    }
    
    static std::vector< Bordered_algebra > get_bordered_algebras(const std::vector< Morse_event >& morse_events) {
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
        algebras[i].orientations =
          morse_events[i].upper_orientations(
            algebras[i + 1].orientations, algebras[i].matchings);
      }
      return algebras;
    }
    
    /* Return an instance of the i^th Morse event, as listed in the template,
     * constructed using arguments.
     */
    static Morse_event instance(const int i, const std::vector< boost::any >& args) {
      return instance_aux_< 0, Morse_events... >(i, args);
    }
    
   private:
    /* I need to add a dummy class template because I can't partially
     * specialize for the initialization otherwise.
     * The other possibility is requiring at least one Morse event, but
     * I want to leave the possibility of having no Morse events.
     */
    template<
      int,
      template< class > class Morse_events_head,
      template< class > class ...Morse_events_tail >
    static Morse_event instance_aux_(const int i, const std::vector< boost::any >& args) {
      if (i == 0) {
        return Morse_event(Morse_events_head< D_module >(args));
      }
      else {
        return instance_aux_< 0, Morse_events_tail... >(i - 1, args);
      }
    }
    
    // Initialization: 
    template< int >
    static Morse_event instance_aux_(const int, const std::vector< boost::any >&) {
      return Morse_event();
    }
  };
  
 public:
  /* TeXify */
  void TeXify(std::ofstream& write_file) const {
    std::cout << "[kd] Drawing knot..." << std::endl;
    
    const auto margins = get_margins_();
    const auto algebras = Detail_<>::get_bordered_algebras(morse_events_default_);
    
    write_file << "\\KnotDiagram{" << std::endl;
    for (int i = 0; i < morse_events_default_.size(); ++i) {
      write_file << morse_events_default_[i].to_string(margins[i], {algebras[i].n_strands, algebras[i + 1].n_strands}) << std::endl;
	}
    write_file << "}" << std::flush;
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
  
  std::vector< std::pair< int, std::vector< boost::any > > > morse_event_data_;
  std::vector< Morse_event< D_module_default > > morse_events_default_;
};