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

#ifndef EVENT_REGISTERS_H_
#define EVENT_REGISTERS_H_

#include <memory>  // unique_ptr

#include "Knot_slice.h"
#include "DA_bimodule.h"

#include "Events/Event_crossing.h"
#include "Events/Event_global_minimum.h"
#include "Events/Event_local_maximum.h"
#include "Events/Event_local_minimum.h"

/* Knot slice register.
 * 
 * A struct with only one member, a static function that takes an event ID
 * and produces a unique pointer to a knot slice.
 */
struct Knot_slice_register {
  using Knot_slice_handle = std::unique_ptr< Knot_slice >;
  
  static Knot_slice_handle knot_slice(Event event_id, int position, int n_strands) {
    switch (event_id) {
//    case Event:: /* event name */: return Knot_slice_handle(new /* knot slice class name */ (n_strands, position));
      case Event::positive_crossing: return Knot_slice_handle(new Knot_slice_crossing< false >(n_strands, position));
      case Event::negative_crossing: return Knot_slice_handle(new Knot_slice_crossing< true > (n_strands, position));
      case Event::local_maximum    : return Knot_slice_handle(new Knot_slice_local_maximum    (n_strands, position));
      case Event::local_minimum    : return Knot_slice_handle(new Knot_slice_local_minimum    (n_strands, position));
      case Event::global_minimum   : return Knot_slice_handle(new Knot_slice_global_minimum   (n_strands, position));
    }
    std::cout << "Incorrect event ID. Exiting..." << std::endl;
    exit(1);  // give up
  }
};

/* DA-bimodule register.
 * 
 * A struct with only one member, a static function that takes a knot slice
 * and produces a unique pointer to a DA-bimodule.
 */
template< class D_module >
struct DA_bimodule_register {
  using Base_DA_bimodule = DA_bimodule< D_module >;
  using Knot_slice_handle = Knot_slice_register::Knot_slice_handle;
  using DA_bimodule_handle = std::unique_ptr< Base_DA_bimodule >;
  
  template< class Algebra >
  static DA_bimodule_handle da_bimodule(const Knot_slice_handle& slice, const Algebra& upper_algebra, const Algebra& lower_algebra) {
    switch(slice->event_id()) {
//    case Event:: /* event name */: return DA_bimodule_handle(new /* templated DA-bimodule class name */ (slice, upper_algebra, lower_algebra));
      case Event::positive_crossing: return DA_bimodule_handle(new DA_bimodule_crossing< false, D_module >(slice, upper_algebra, lower_algebra));
      case Event::negative_crossing: return DA_bimodule_handle(new DA_bimodule_crossing< true,  D_module >(slice, upper_algebra, lower_algebra));
      case Event::local_maximum    : return DA_bimodule_handle(new DA_bimodule_local_maximum<   D_module >(slice, upper_algebra, lower_algebra));
      case Event::local_minimum    : return DA_bimodule_handle(new DA_bimodule_local_minimum<   D_module >(slice, upper_algebra, lower_algebra));
      case Event::global_minimum   : return DA_bimodule_handle(new DA_bimodule_global_minimum<  D_module >(slice, upper_algebra, lower_algebra));
    }
    std::cout << "Incorrect event ID. Exiting..." << std::endl;
    exit(1);  // give up
  }
};

#endif  // EVENT_REGISTERS_H_