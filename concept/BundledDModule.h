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

#include <iostream>

class unspecified;

/* Concept for implementations of bundled D-modules over a path algebra of
 * characteristic 2, up to homotopy equivalence.
 */
template< class BundledDModuleOptions >
class BundledDModule {
 public:
  using IdempotentElement  = typename BundledDModuleOptions::IdempotentElement;
  using GeneratorBundle    = typename BundledDModuleOptions::GeneratorBundle;
  using CoefficientBundle  = typename BundledDModuleOptions::CoefficientBundle;
  using BorderedAlgebra    = typename BundledDModuleOptions::BorderedAlgebra;
  using AlgebraElement     = typename BundledDModuleOptions::AlgebraElement;
  using GeneratorType      = typename BundledDModuleOptions::GeneratorType;
  using PoincarePolynomial = typename BundledDModuleOptions::PoincarePolynomial;
  
  using GeneratorBundleHandle        = unspecified;
  using CoefficientBundleHandle      = unspecified;
  using CoefficientBundleHandleRange = unspecified;
  
  /* SECTION: Constructors */
  BundledDModule() = default;
  BundledDModule(BundledDModule&) = delete;
  BundledDModule& operator=(BundledDModule&) = delete;
  BundledDModule(BundledDModule&&) = default;
  BundledDModule& operator=(BundledDModule&&) = default;
  ~BundledDModule() = default;
  
  /* SECTION: Handle access functions */
  
  /* Get the lower idempotent of the generator bundle */
  const IdempotentElement idem(const GeneratorBundleHandle&) const;
  
  /* Get the lower algebra element of the coefficient bundle
   */
  const AlgebraElement alg_el(const CoefficientBundleHandle&) const;
  
  /* Get the coefficient bundles arriving at a source of the agglomerated
   * coefficient
   */
  const CoefficientBundleHandleRange coefs_to_source(const CoefficientBundleHandle&) const;
  // and the analogous functions, probably
  
  /* SECTION: Generator and coefficient bundle creation functions during
   * box tensor product
   */
  
  /* Create a generator bundle corresponding to all generators of the
   * form new_marking \otimes old_markings, with old_markings in old_gen. The
   * new idempotent for this generator bundle is new_idempotent.
   */
  void add_gen_bundle(const IdempotentElement& new_idempotent,
                      const GeneratorType& new_marking,
                      const GeneratorBundleHandle& old_gen_bundle);
  /* Overloaded function when there is no old generator bundle to copy over
   */
  void add_gen_bundle(const IdempotentElement& new_idempotent);
  
  /* Create an coefficient bundle from an upper coefficient bundle, given the
   * markings of the lowest slice and the new algebra element
   */
  void add_coef_bundle(const AlgebraElement& new_algebra_element,
                       const GeneratorType& back_marking,
                       const GeneratorType& front_marking,
                       const CoefficientBundle& upper_coef_bundle);
  /* Overloaded function when there is no upper coefficient bundle, only an
   * upper idempotent
   */
  void add_coef_bundle(const AlgebraElement& new_algebra_element,
                       const GeneratorType& back_marking,
                       const GeneratorType& front_marking,
                       const IdempotentElement& upper_idempotent);
  
  /* SECTION: Mathematical manipulation of the represented D-module */
  
  /* Set the represented D-module to be trivial (over the trivial algebra, with
   * a single generator).
   */
  void set_as_trivial();
  
  /* Dualize the represented D-module */
  void dualize();
  
  /* Calculate a possibly smaller, homotopically equivalent D-module. Return
   * true if further reduction is possible, false otherwise.
   */
  bool reduce();
  
  /* Return the Poincaré polynomial of the D-module, assuming that the
   * structure morphism is null. Behavior is undefined otherwise.
   */
  PoincarePolynomial poincare_polynomial();
  
  /* Write something to output_file that can be compiled by LaTeX
   */
  void TeXify(std::ostream& output_file);
};