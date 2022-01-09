#ifndef KNOT_INTERFACE_H_
#define KNOT_INTERFACE_H_

#include <array>
#include <cstdlib>  // strtol
#include <fstream>
#include <iostream>
#include <sstream>  // stringstream
#include <string>
#include <tuple>  // tie
#include <utility>  // pair
#include <vector>

namespace ComputeHFKv2 {
  #include "ComputeHFKv2/Diagrams.h"  // from ComputeHFKv2
  #include "ComputeHFKv2/Diagrams.cpp"
}
#include "link/link.h"  // from Regina

#include "Differential_suffix_forest/Differential_suffix_forest.h"
#include "Knot_diagram/Knot_diagram.h"

#include "Morse_event/Positive_crossing.h"
#include "Morse_event/Negative_crossing.h"
#include "Morse_event/Local_maximum.h"
#include "Morse_event/Local_minimum.h"
#include "Morse_event/Global_minimum.h"

#include "Math_tools/Poincare_polynomial.h"

/* An interface between Regina knots and our knots, passing by ComputeHFKv2
 */
class Knot_interface {
 public:
  using Knot_diagram = Knot_diagram<
    Positive_crossing,
    Negative_crossing,
    Local_maximum,
    Local_minimum,
    Global_minimum
  >;
  using Morse_data_container = typename Knot_diagram::Morse_data_container;
  
  enum : int {
    positive_crossing,
    negative_crossing,
    local_maximum,
    local_minimum,
    global_minimum
  };
  
  Knot_interface() { }
  
  void import_regina_signature(const std::string& knot_sig) {
    {  // swap in
      regina::Link other_knot(knot_sig);
      regina_knot_.swapContents(other_knot);
    }
    planar_diagram_string_ = planar_diagram_to_string_(regina_to_planar_diagram_(regina_knot_));
    
    import_planar_diagram(planar_diagram_string_);
  }
  
  void import_planar_diagram(const std::string& planar_diagram_string) {
    planar_diagram_string_ = planar_diagram_string;
    ComputeHFKv2::PlanarDiagram planar_diagram(planar_diagram_string_);
    compute_hfk_v2_morse_code_ = planar_diagram.GetSmallGirthMorseCode();
    auto morse_data = morse_code_to_data_(compute_hfk_v2_morse_code_);
    knot_diagram_.import_data(morse_data);
  }
  
  regina::Link regina_knot() const {
    return regina_knot_;
  }
  
  std::string planar_diagram_string() const {
    return planar_diagram_string_;
  }
  
  ComputeHFKv2::MorseCode compute_hfk_v2_morse_code() const {
    return compute_hfk_v2_morse_code_;
  }
  
  const Knot_diagram& knot_diagram() const {  // return a viewing reference
    return knot_diagram_;
  }
  
  template< template< class > class D_module = Forest >
  Poincare_polynomial knot_Floer_homology() const {
    if (knot_diagram_.max_n_strands() <= 31) {
      return knot_diagram_.knot_Floer_homology< Poincare_polynomial, D_module< Forest_options_default_short > >();
    }
    else {
      return knot_diagram_.knot_Floer_homology< Poincare_polynomial, D_module< Forest_options_default_long > >();
    }
  }
    
 private:
  /* For planar diagrams, we order the strands as follows:
   * 
   *    2
   *    ^
   *    |   
   * 3 --- 1
   *    |   
   *    0
   * 
   */
  std::vector< std::array< int, 4 > > regina_to_planar_diagram_(const regina::Link& knot) const {
    std::vector< std::array< int, 4 > > planar_diagram;
    std::string oriented_Gauss_code = knot.orientedGauss();
    
    int next_crossing = 0;  // crossings start at 0
    int last_arc = 1;  // arcs start at 1 because mathematicians...
    
    bool negative;
    bool left;
    int crossing;
    int in_direction, out_direction;  // directions of the arcs we are following
    
    for (int pos = 0; pos < oriented_Gauss_code.size(); ) {
      std::tie(negative, left, crossing) = unpack_(oriented_Gauss_code, pos);
      
      if (crossing > next_crossing) {
        std::cout << "[ki] crossing too big" << std::endl;
        continue;
      }
      else if (crossing == next_crossing) {
        planar_diagram.emplace_back();
        ++next_crossing;
      }
      
      if (negative) {
        in_direction = 0;
        out_direction = 2;
      }
      else {
        if (left) {
          in_direction = 3;
          out_direction = 1;
        }
        else {
          in_direction = 1;
          out_direction = 3;
        }
      }
      
      planar_diagram[crossing][in_direction] = last_arc++;
      planar_diagram[crossing][out_direction] = last_arc;
    }
    /* Fix the last arc */
    planar_diagram[crossing][out_direction] = 1;
    return planar_diagram;
  }  // get_planar_diagram_
  
  std::tuple< bool, bool, int > unpack_(const std::string& code, int& pos) const {
    bool negative = (code[pos] == '-');
    bool left = (code[pos + 1] == '<');
    int crossing = strtol(&code[pos + 2], nullptr, 10) - 1;  // strings are contiguous, right?
    pos = code.find_first_not_of("0123456789 ", pos+2);
    if (pos == std::string::npos) { pos = code.size(); }  // end was reached
    return std::make_tuple(negative, left, crossing);
  }  // unpack_
  
  std::string planar_diagram_to_string_(const std::vector< std::array< int, 4 > >& planar_diagram) const {
    std::stringstream ss;
    ss << "PD[";
    for (auto& a : planar_diagram) {
      ss << "X[" << a[0] << "," << a[1] << "," << a[2] << "," << a[3] << "],";
    }
    ss.seekp(-1, ss.cur);
    ss << "]" << std::flush;
    
    return ss.str();
  }  // planar_diagram_to_string_
  
  Morse_data_container morse_code_to_data_(const ComputeHFKv2::MorseCode& morse_code) const {
    std::vector< int > raw_morse_data = morse_code.GetMorseList();
    
    Morse_data_container morse_data;
    for (int i = 0; i < raw_morse_data.size() - 1; ++i) {
      if (raw_morse_data[i] == 1000 or raw_morse_data[i] == 1001) {
        morse_data.push_back({local_maximum, {raw_morse_data[++i] - 1}});
      }
      else if (raw_morse_data[i] == -1000 or raw_morse_data[i] == -1001) {
        morse_data.push_back({local_minimum, {0}});
      }
      else if (raw_morse_data[i] > 0) {
        morse_data.push_back({positive_crossing, {raw_morse_data[i] - 1}});
      }
      else {
        morse_data.push_back({negative_crossing, {-raw_morse_data[i] - 1}});
      }
    }
    
    // last event is always global minimum
    morse_data.push_back({global_minimum, {0}});
    
    return morse_data;
  }
  
  regina::Link regina_knot_;  // from Regina
  std::string planar_diagram_string_;  // i.e. PD[X[a,b,c,d],...]
  ComputeHFKv2::MorseCode compute_hfk_v2_morse_code_;  // from ComputeHFKv2
  Knot_diagram knot_diagram_;  // our knot diagram
};

#endif  // KNOT_INTERFACE_H_