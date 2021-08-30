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

#ifndef DA_BIMODULE_H_
#define DA_BIMODULE_H_

#include <string>

#include "Math_tools/Poincare_polynomial.h"

/* DA-bimodule.
 * 
 * This is an abstract class implementing DA-bimodules.
 * 
 * DA-bimodules are only implemented in their relation to D-modules, so it
 * makes sense to be templated by the D-module.
 * 
 * Objects of this class are characterized by a member function
 * box_tensor_product, which calculates the box tensor product of the knot
 * slice with a given upper D-module.
 */
template< class D_module >
class DA_bimodule {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Algebra = typename D_module::Bordered_algebra;
  using Alg_el = typename D_module::Alg_el;
  using Coef_bundle = typename D_module::Coef_bundle;
  using Polynomial = typename D_module::Polynomial;
  using Monomial = typename Polynomial::Monomial;  // currently, pair of int
  
  DA_bimodule() { }
  
  DA_bimodule(Algebra upper_algebra, Algebra lower_algebra) :
    upper_algebra_(upper_algebra),
    lower_algebra_(lower_algebra)
  { }
  
  virtual ~DA_bimodule() = default;
  
  /* Double Alexander weight and Maslov weight */
  virtual Monomial alexander_maslov_weights(const Gen_type generator_type) const = 0;
  
  /* Given a generator type, give the corresponding LaTeX math label. This is
   * based on notations given by mathematicians.
   */
  virtual std::string get_label(const Gen_type generator_type) const = 0;
  
  /* This should produce a valid D-module, but currently the output does not
   * have connected nodes and arcs...
   */
  virtual D_module box_tensor_product(const D_module&) const = 0;
  
  D_module operator()(const D_module& d_module) const {
    auto new_d_module = box_tensor_product(d_module);
    new_d_module.post_construction_processing();  // find a better way to do this?
    return new_d_module;
  }
  
 protected:
  Algebra upper_algebra_;  // These members really belong to the ABC
  Algebra lower_algebra_;
};

/* The box tensor product is not a member function. It appears as a standalone
 * function in some namespace (to be defined).
 */
template< class D_module, class DA_bimodule >
D_module box_tensor_product(const D_module& d_module, const DA_bimodule& da_bimodule) {
  return da_bimodule(d_module);
}

#endif  // DA_BIMODULE_H_