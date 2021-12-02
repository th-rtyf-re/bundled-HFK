#include <fstream>
#include <iostream>
#include <string>  // stoi
#include <utility>  // pair
#include <vector>

#include <boost/any.hpp>

#include "test_morse_event.h"
#include "test_da_bimodule.h"

template< class ...Ts >
struct D_module_pack;

template< template< class > class ...Ts >
struct Morse_event_pack;

template< class ...Ts >
class Knot_diagram;

// Need at least one D_module and one Morse event
template< class D_module_default, class ...D_module_others, template< class > class ...Morse_events >
class Knot_diagram< D_module_pack< D_module_default, D_module_others... >, Morse_event_pack< Morse_events... > > {
 public:
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
  }  // import_csv
  
  int max_n_strands() const {
    int n_strands = 0;
    int max_n_strands = 0;
    
    const auto morse_events = Detail_< D_module_default >::morse_events_();
    
    for (const auto morse_event : morse_events) {
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
#ifdef DRAW
    std::ofstream suffix_forest("differential_suffix_forest.tex");
#endif  // DRAW
    
    const auto da_bimodules = Detail_< D_module >::da_bimodules_(morse_event_data_);
    
    //std::ofstream write_file("differential_suffix_forest.tex");  // debug LaTeX file
    D_module d_module;
    d_module.set_as_trivial();
    
    std::cout << "[kd] D-module is trivial" << std::endl;
    
    // box tensor product for each Morse event
    for (const auto& da_bimodule : da_bimodules) {
      std::cout << "[kd] start loop" << std::endl;
      d_module = box_tensor_product(da_bimodule, d_module);
      std::cout << "[kd] tensor done" << std::endl;
      d_module.reduce();
#ifdef DRAW
      d_module.TeXify(suffix_forest);
#endif  // DRAW
    }
    
#ifdef DRAW
    suffix_forest.close();
#endif  // DRAW
    
    std::cout << "[kd] math done" << std::endl;
    
    return d_module.template poincare_polynomial< Polynomial >();
  }
  
 private:
  
  /* All the private methods depend on the choice of D-modules. In order to
   * make things hopefully more readable, I've put all these methods in a
   * private struct Detail_, templated by D_module, as static methods.
   */
  template< class D_module >
  struct Detail_ {
    using Bordered_algebra = typename D_module::Bordered_algebra;
    using Morse_event = Morse_event< D_module >;
    using DA_bimodule = DA_bimodule< Morse_event, D_module >;
    
    static std::vector< DA_bimodule > da_bimodules_(const std::vector< std::pair< int, std::vector< boost::any > > >& morse_event_data) {
      std::vector< DA_bimodule > da_bimodules;
      
      const auto morse_events = morse_events_(morse_event_data);
      const auto algebras = bordered_algebras_(morse_events);
      
      for (int i = 0; i != morse_events.size(); ++i) {
        da_bimodules.emplace_back(morse_events[i],
                                  algebras[i],
                                  algebras[i + 1]);
      }
      
      return da_bimodules;
    }
    
    static std::vector< Morse_event > morse_events_(const std::vector< std::pair< int, std::vector< boost::any > > >& morse_event_data) {
      std::vector< Morse_event > morse_events;
      
      for (const auto value_pair : morse_event_data) {
        morse_events.push_back(instance_< Morse_events... >(value_pair.first, value_pair.second));
      }
      
      return morse_events;
    }
    
    template<
      template< class > class Morse_events_head,
      template< class > class Morse_events_neck,
      template< class > class ...Morse_events_tail >
    static Morse_event instance_(const int i, const std::vector< boost::any >& args) {
      if (i == 0) {
        return Morse_event(Morse_events_head< D_module >(args));
      }
      else {
        return instance_< Morse_events_neck, Morse_events_tail... >(i - 1, args);
      }
    }
    
    // Initialization
    template< template< class > class Morse_events_head >
    static Morse_event instance_(const int, const std::vector< boost::any >& args) {
      return Morse_event(Morse_events_head< D_module >(args));
    }
    
    static std::vector< Bordered_algebra > bordered_algebras_(const std::vector< Morse_event >& morse_events) {
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
      algebras.back().orientations = {false, false, true, true};  // temporary debug stuff
      for (int i = morse_events.size() - 1; i >= 0; --i) {
        algebras[i].orientations =
          morse_events[i].upper_orientations(
            algebras[i + 1].orientations, algebras[i].matchings);
      }
      return algebras;
    }
  };
  
  std::vector< std::pair< int, std::vector< boost::any > > > morse_event_data_;
};