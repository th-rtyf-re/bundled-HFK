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

#ifndef DIFFERENTIAL_SUFFIX_FOREST_OPTIONS_H_
#define DIFFERENTIAL_SUFFIX_FOREST_OPTIONS_H_

#include <vector>
#include <utility>  // pair

#include "Bordered_algebra/Bordered_algebra.h"
#include "Bordered_algebra/Idempotent.h"

struct Forest_options_default_short {
  using Idem = Idempotent_short;
  using Bordered_algebra = Bordered_algebra< Idem >;
  using Alg_el = typename Bordered_algebra::Element;
  using Gen_type = unsigned char;  // no need to pass by reference excessively
  using Weights = std::pair< int, int >;
};

struct Forest_options_default_long {
  using Idem = Idempotent_long< std::vector< bool > >;
  using Bordered_algebra = Bordered_algebra< Idem >;
  using Alg_el = typename Bordered_algebra::Element;
  using Gen_type = unsigned char;  // no need to pass by reference excessively
  using Weights = std::pair< int, int >;
};

#endif  // DIFFERENTIAL_SUFFIX_FOREST_OPTIONS_H_