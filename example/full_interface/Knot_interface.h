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

#ifndef KNOT_INTERFACE_H_
#define KNOT_INTERFACE_H_

#include <array>
#include <cstdlib>  // strtol
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>  // tie
#include <vector>

#include "link/link.h"  // from Regina
#include "Planar_diagram.h"

#include "Differential_suffix_forest/Differential_suffix_forest.h"
#include "Differential_suffix_forest/Differential_suffix_forest_options.h"
#include "Knot_diagram/Knot_diagram.h"

#include "Morse_event/Positive_crossing.h"
#include "Morse_event/Negative_crossing.h"
#include "Morse_event/Local_maximum.h"
#include "Morse_event/Local_minimum.h"
#include "Morse_event/Global_minimum.h"

#include "Math_tools/Poincare_polynomial.h"

/* An interface between Regina knots and our knots, passing by planar diagrams
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
    planar_diagram_ = Planar_diagram(regina_to_planar_diagram_(regina_knot_));
    srand(0);
    auto legacy_morse_code = planar_diagram_.get_legacy_morse_code();
    auto morse_data = morse_code_to_data_(legacy_morse_code);
    knot_diagram_.import_data(morse_data);
  }
  
  void import_planar_diagram(const std::string& planar_diagram_string) {
    planar_diagram_ = Planar_diagram(planar_diagram_string);
    srand(0);  // getting Morse code is random; fix the seed for consistency
    auto legacy_morse_code = planar_diagram_.get_legacy_morse_code();
    auto morse_data = morse_code_to_data_(legacy_morse_code);
    knot_diagram_.import_data(morse_data);
  }
  
  void import_morse_events(std::ifstream& in_file) {
    knot_diagram_.import_csv(in_file);
  }
  
  regina::Link regina_knot() const {
    return regina_knot_;
  }
  
  std::string planar_diagram_string() const {
    return planar_diagram_.to_string();
  }  // planar_diagram_to_string_
  
  const Knot_diagram& knot_diagram() const {  // return a viewing reference
    return knot_diagram_;
  }
  
  Poincare_polynomial knot_Floer_homology() const {
    if (knot_diagram_.max_n_strands() <= 31) {
      return knot_diagram_.knot_Floer_homology< Poincare_polynomial, Differential_suffix_forest< Forest_options_default_short > >();
    }
    else {
      return knot_diagram_.knot_Floer_homology< Poincare_polynomial, Differential_suffix_forest< Forest_options_default_long > >();
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
  std::vector< int > regina_to_planar_diagram_(const regina::Link& knot) const {
    std::vector< int > planar_diagram;
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
        planar_diagram.resize(planar_diagram.size() + 4, 0);
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
      
      planar_diagram[4 * crossing + in_direction] = last_arc++;
      planar_diagram[4 * crossing + out_direction] = last_arc;
    }
    /* Fix the last arc */
    planar_diagram[4 * crossing + out_direction] = 1;
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
  
  
  Morse_data_container morse_code_to_data_(const std::vector< int >& legacy_morse_code) const {
    Morse_data_container morse_data;
    for (int i = 0; i < legacy_morse_code.size() - 1; ++i) {
      if (legacy_morse_code[i] == 1000 or legacy_morse_code[i] == 1001) {
        morse_data.push_back({local_maximum, {legacy_morse_code[++i] - 1}});
      }
      else if (legacy_morse_code[i] == -1000 or legacy_morse_code[i] == -1001) {
        morse_data.push_back({local_minimum, {0}});
      }
      else if (legacy_morse_code[i] > 0) {
        morse_data.push_back({positive_crossing, {legacy_morse_code[i] - 1}});
      }
      else {
        morse_data.push_back({negative_crossing, {-legacy_morse_code[i] - 1}});
      }
    }
    
    // last event is always global minimum
    morse_data.push_back({global_minimum, {0}});
    
    return morse_data;
  }
  
  regina::Link regina_knot_;  // from Regina
  Planar_diagram planar_diagram_;
  Knot_diagram knot_diagram_;  // our knot diagram
};

#endif  // KNOT_INTERFACE_H_