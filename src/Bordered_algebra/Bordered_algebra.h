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

#ifndef BORDERED_ALGEBRA_H_
#define BORDERED_ALGEBRA_H_

#include <algorithm>  // none_of
#include <iostream>
#include <vector>

/* Bordered algebra.
 * 
 * Implementation of bordered algebras for bordered knot Floer homology.
 * Recall that these algebras are path algebras.
 * 
 * A bordered algebra is a struct specifying the type of idempotent, the
 * number of strands, as well as matchings and orientations of the strands.
 * 
 * This struct has a member class, Element, which implements algebra elements.
 * The implementation of bordered algebras does not depend on the implementation
 * of the idempotent ring. However, the implementation of algebra elements does.
 */
template< class Idempotent >
struct Bordered_algebra {
  using Idem = Idempotent;
  
  int n_strands;
  
  std::vector< int > matchings;
  
  std::vector< bool > orientations;
  
  /* Element subclass.
   * 
   * This is our representation for algebra elements, and more specifically
   * algebra monomials. There is no specific way to represent the null monomial;
   * the best we can do is have source and target idempotents that are too far
   * apart, or a generating interval in the U weights.
   */
  class Element {
   public:
    Element(Idem source_idem, Idem target_idem, std::vector< int > U_weights) :
      source_idem_(source_idem),
      target_idem_(target_idem),
      U_weights_(U_weights)
    { }
    
    Element(Idem source_idem, Idem target_idem) :
      source_idem_(source_idem),
      target_idem_(target_idem),
      U_weights_(source_idem.size() - 1, 0)  // default to null U weights
    { }
    
    Idem source_idem() const {
      return source_idem_;
    }
    
    Idem target_idem() const {
      return target_idem_;
    }
    
    std::vector< int > U_weights() const {
      return U_weights_;
    }
    
    void dualize() {
      source_idem_.swap(target_idem_);
    }
    
    int U_weight(int position) const {
      return U_weights_[position];
    }
    
    int& U_weight(int position) {
      return U_weights_[position];
    }
    
    /* Product of algebra elements without checks.
     * This only produces a defined result if the real product is nonzero.
     */
    Element operator*(const Element& other) const {
      std::vector< int > product_U_weights(U_weights_);
      int back_weight = 0;
      int mid_weight = 0;
      int front_weight = 0;
      
      for (int i = 0; i < product_U_weights.size(); ++i) {
        if (source_idem_[i]) { ++back_weight; }
        if (target_idem_[i]) { ++mid_weight; }
        if (other.target_idem_[i]) { ++front_weight; }
        if ((back_weight > mid_weight and mid_weight < front_weight)
            or (back_weight < mid_weight and mid_weight > front_weight)) {
          product_U_weights[i] += other.U_weights_[i] + 1;  // L_iR_i = U_i = R_iL_i
        }
        else {
          product_U_weights[i] += other.U_weights_[i];
        }
      }
      return Element(source_idem_, other.target_idem_, product_U_weights);
    }  // operator*
    
    bool operator==(const Element& other) const {
      return (source_idem_ == other.source_idem_ and target_idem_ == other.target_idem_
              and U_weights_ == other.U_weights_);
    }
    
    /* is_null. Taken from ComputeHFKv2/Utility.cpp, NonZero
     * 
     * Pre-condition: source and target idems are close enough.
     * 
     * The idea is to scan the combinatorial inputs for "generating intervals" as
     * described in [OzsvathSzabo2018, Definition 3.6].
     */
    bool is_null() const {
      bool interval(true);  // flag indicating that we have started an interval
      bool LR_i(false);     // flag indicating that the algebra element has an L_i or R_i
      bool U_i;             // flag indicating that the algebra element has a U_i
      bool source_i;  // i^th bit of the source idem
      bool target_i;  // i^th bit of the target idem
      
      for (int i = 0; i < U_weights_.size(); ++i) {
        U_i = U_weights_[i];
        source_i = source_idem_[i + 1];
        target_i = target_idem_[i + 1];
        
        if (interval and (!source_i or !target_i) and U_i) {  // end of interval
          return true;
        }
        else if ((LR_i and (source_i != target_i))
                 or (!LR_i and !source_i and !target_i)) {  // start of interval
          interval = true;
        } 
        else if (interval and !LR_i and source_i and target_i and U_i) {  // continue the interval
          //interval = true;  // We don't actually need this instruction
        }
        else {  // reset the interval flag
          interval = false;
        }
        LR_i = (LR_i != (source_i != target_i));
      }
      return false;
    }  // is null
    
    bool is_invertible() const {
      return ((source_idem_ == target_idem_)
        and std::none_of(U_weights_.begin(), U_weights_.end(), [](int n) { return n; }));
    }
    
    /* Produce the algebra element associated to the differential cell, in LaTeX
     * math format. Mostly for display and debugging purposes.
     */
    std::string to_string() const {
      int current_difference = 0;
      
      /* Add LR part */
      std::string L_string = "";
      std::string R_string = "";
      for (int i = 0; i < source_idem_.size(); ++i) {
        current_difference = current_difference + source_idem_[i]
                                                - target_idem_[i];
        if (current_difference == 1) {
          R_string = "R_{" + std::to_string(i) + "}" + R_string;
        }
        else if (current_difference == -1) {
          L_string = L_string + "L_{" + std::to_string(i) + "}";
        }
        else if (current_difference == 0) { }
        else {
          std::cout << "[bae] idempotents are too far apart" << std::endl;
        }
      }
      
      /* Add U part */
      std::string LRU_string = L_string + R_string;
      for (int i = 0; i < U_weights_.size(); ++i) {
        if (U_weights_[i] > 0) {
          LRU_string = LRU_string + "U_{" + std::to_string(i) + "}";
          if (U_weights_[i] > 1) {
            LRU_string = LRU_string + "^{" + std::to_string(U_weights_[i]) + "}";
          }
        }
      }
      
      if (LRU_string.size() == 0) {
        LRU_string = "1";
      }
      return LRU_string;
    }  // to_string
    
    friend std::ostream& operator<<(std::ostream& os, const typename Bordered_algebra< Idempotent >::Element& el) {
      os << el.to_string();
      return os;
    }
    
   private:
    Idem source_idem_;
    Idem target_idem_;
    
    std::vector< int > U_weights_;
  };  // Element
};  // Bordered_algebra


#endif  // BORDERED_ALGEBRA_H_