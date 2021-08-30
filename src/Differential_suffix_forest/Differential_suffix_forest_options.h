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

#ifndef DIFFERENTIAL_SUFFIX_FOREST_OPTIONS_H_
#define DIFFERENTIAL_SUFFIX_FOREST_OPTIONS_H_

#include "Differential_arc.h"
#include "Differential_suffix_forest_node.h"

#include "Bordered_algebra/Bordered_algebra.h"
#include "Bordered_algebra/Idempotent.h"
#include "Math_tools/Poincare_polynomial.h"

/* Differential suffix forest options.
 * 
 * These structs define many of the classes needed for differential suffix
 * forests, to make them accessible to nodes, arcs, and forests, as well as to
 * define options for short idempotents and for long idempotents.
 */
struct Forest_options_default_short {
  using Idem = Idempotent_short;
  using Bordered_algebra = Bordered_algebra< Idem >;
  using Alg_el = typename Bordered_algebra::Element;
  using Polynomial = Poincare_polynomial;
  using Gen_type = unsigned char;  // no need to pass by reference excessively
  
  using Node = Differential_suffix_forest_node< Forest_options_default_short >;
  using Node_container = std::map< const Idem, Node >;
  using Arc = Differential_arc< Forest_options_default_short >;
  using Arc_container = std::list< Arc >;
};

struct Forest_options_default_long {
  using Idem = Idempotent_long< std::vector< bool > >;
  using Bordered_algebra = Bordered_algebra< Idem >;
  using Alg_el = typename Bordered_algebra::Element;
  using Polynomial = Poincare_polynomial;
  using Gen_type = unsigned char;  // no need to pass by reference excessively
  
  using Node = Differential_suffix_forest_node< Forest_options_default_long >;
  using Node_container = std::map< const Idem, Node >;
  using Arc = Differential_arc< Forest_options_default_long >;
  using Arc_container = std::list< Arc >;
};

#endif  // DIFFERENTIAL_SUFFIX_FOREST_OPTIONS_H_