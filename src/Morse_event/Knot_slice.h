/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *  Bundled HFK - a knot Floer homology calculator                           *
 *                                                                           *
 *  Copyright (C) 2021  Isaac Ren                                            *
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

#ifndef KNOT_SLICE_H_
#define KNOT_SLICE_H_

#include <string>
#include <utility>  // pair
#include <vector>

#include "Morse_event.h"

/* Knot slice.
 * 
 * This is out implementation of one-dimensional oriented manifolds with
 * boundary, specifically for knot slices.
 * 
 * A knot slice naturally comes with the position parameter, because we are
 * are working with slices made up of many strands, but where events are
 * bounded in width.
 * (To be honest, I'm not satisfied with this justification...)
 * 
 * I think it's better to make the virtual functions pure, so that each event
 * clearly knows what each function does.
 */
class Knot_slice {
 public:
  Knot_slice() { }
  
  Knot_slice(int upper_n_strands, int lower_n_strands, int position) :
    upper_n_strands_(upper_n_strands),
    lower_n_strands_(lower_n_strands),
    position_(position)
  { }
  
  Knot_slice(Knot_slice&) = default;
  Knot_slice& operator=(Knot_slice&) = default;
  Knot_slice(Knot_slice&&) = default;
  Knot_slice& operator=(Knot_slice&&) = default;
  virtual ~Knot_slice() = default;
  
  int upper_n_strands() const {  // never used
    return upper_n_strands_;
  }
  
  int lower_n_strands() const {
    return lower_n_strands_;
  }
  
  int position() const {
    return position_;
  }
  
  /* Return the event id. This must match the Morse_event enum class, as well
   * as the Event registers.
   */
  virtual Morse_event event_id() const = 0;
  /* Given upper matchings for the strands, calculate the lower matchings
   */
  virtual std::vector< int > lower_matchings(std::vector< int > upper_matchings) const = 0;
  /* Given lower orientations and upper matchings, calculate the upper
   * orientations
   */
  virtual std::vector< bool > upper_orientations(std::vector< bool > lower_orientations, const std::vector< int >& upper_matchings) const = 0;
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice
   */
  virtual std::string to_string(const std::pair< int, int >& margins) const = 0;
  /* Horizontal shifts are used when drawing the knot.
   * Update the horizontal shifts for each previous layer, and add the
   * horizontal shift for the current layer. Return the current shift.
   */
  virtual std::pair< int, int > update_margins(std::pair< int, int > prev_margins) const = 0;
  
 protected:
  int upper_n_strands_;
  int lower_n_strands_;
  int position_;
};

#endif  // KNOT_SLICE_H_