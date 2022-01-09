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

#ifndef DA_BIMODULE_H_
#define DA_BIMODULE_H_

#ifdef BUNDLED_HFK_VERBOSE_
#include <iostream>
#endif  // BUNDLED_HFK_VERBOSE_

/* DA-bimodule.
 * 
 * The user defines classes that satisfy the Morse event concept, but it's
 * easier for Knot_diagrams to reason with DA bimodules.
 * These are constructed when we call the knot Floer homology function.
 */
template< class Morse_event >
class DA_bimodule {
 public:
  using D_module = typename Morse_event::D_module;
  using Algebra = typename D_module::Bordered_algebra;
  
  DA_bimodule(
    const Morse_event morse_event,
    const Algebra upper_algebra,
    const Algebra lower_algebra
  ) :
    morse_event(morse_event),
    upper_algebra(upper_algebra),
    lower_algebra(lower_algebra)
  { }
  
  const Morse_event morse_event;
  const Algebra upper_algebra;
  const Algebra lower_algebra;
  
  friend D_module box_tensor_product(
    const DA_bimodule& da_bimodule,
    const D_module& old_d_module
  ) {
    D_module new_d_module;
    da_bimodule.morse_event.tensor_generators(
      new_d_module,
      old_d_module,
      da_bimodule.upper_algebra,
      da_bimodule.lower_algebra
    );
    const auto weights = da_bimodule.morse_event.get_weights(
      da_bimodule.upper_algebra,
      da_bimodule.lower_algebra
    );
    const auto labels = da_bimodule.morse_event.get_labels(
      da_bimodule.upper_algebra,
      da_bimodule.lower_algebra
    );
    new_d_module.lock_generators(old_d_module, weights, labels);
    da_bimodule.morse_event.tensor_coefficients(
      new_d_module,
      old_d_module,
      da_bimodule.upper_algebra,
      da_bimodule.lower_algebra
    );
    new_d_module.lock_coefficients();
    return new_d_module;
  }
  
#ifdef BUNDLED_HFK_VERBOSE_
  friend std::ostream& operator<<(
    std::ostream& os,
    const DA_bimodule& da_bimodule
  ) {
    os << da_bimodule.morse_event;
    return os;
  }
#endif  // BUNDLED_HFK_VERBOSE_
};

#endif  // DA_BIMODULE_H_