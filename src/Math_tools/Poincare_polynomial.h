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

#ifndef POINCARE_POLYNOMIAL_H_
#define POINCARE_POLYNOMIAL_H_

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>  // pair
#include <vector>

/* Poincaré polynomials.
 * 
 * This is our implementation of two-variable polynomials, say t and q, where
 * the first variable has exponents in 1/2Z, and all coefficients are non-
 * negative
 */
class Poincare_polynomial {
 public:
  using Poly = Poincare_polynomial;  // rename self
  using Monomial = std::pair< int, int >;
  using Coefficient = int;
  using Container = std::vector< std::pair< Monomial, Coefficient > >;
  
  Poincare_polynomial() { }
  
  Poincare_polynomial(Container& terms) : terms_(terms) { }
  
  Poincare_polynomial(const int n) :  // constructor from integer
    terms_(n ? Container({{{0, 0}, n}}) : Container())
  { }
  
  Poincare_polynomial& operator=(const int n) {  // assignment from integer
    if (n == 0) {
      terms_.clear();
    }
    else {
      terms_ = {{{0, 0}, n}};
    }
    return *this;
  }
  
  bool is_null() {
    return (terms_.empty());
  }
  
  Poly operator+(const Poly& other) {
    Container sum;
    auto this_it = terms_.begin();
    auto other_it = other.terms_.begin();
    while (this_it != terms_.end() and other_it != other.terms_.end()) {
      if (this_it->first < other_it->first) {
        sum.push_back(*this_it);
        ++this_it;
      }
      else if (this_it->first > other_it->first) {
        sum.push_back(*other_it);
        ++other_it;
      }
      else {  // equal monomials
        sum.push_back({this_it->first, this_it->second + other_it->second});
        ++this_it;
        ++other_it;
      }
    }
    if (this_it != terms_.end()) {
      sum.insert(sum.end(), this_it, terms_.end());
    }
    if (other_it != other.terms_.end()) {
      sum.insert(sum.end(), other_it, other.terms_.end());
    }
    
    return Poly(sum);
  }
  
  Poly& operator+=(const Poly& other) {
    *this = operator+(other);
    return *this;
  }
  
  /* A monomial-like object has two member coefficients, first and second,
   * where first is in 1/2Z and second is in Z.
   */
  template< class Monomial_like >
  Poly& operator*=(const Monomial_like& monomial) {
    for (auto& term : terms_) {
      term.first.first += monomial.first;
      term.first.second += monomial.second;
    }
    return *this;
  }
  
  // plus_times?
  
  /* TeXify */
  
  /* TeXify as if the first variable is \sqrt{t} */
  std::string to_string() const {
    if (terms_.empty()) {
      return "0";
    }
    else {
      std::stringstream stream;
      for (auto& term : terms_) {
        if (term.second != 1) {
          stream << term.second;
        }
        if (term.first.first != 0) {
          stream << "t";
          if (term.first.first % 2 == 0) {  // even Alexander weight
            if (term.first.first != 2) {
              stream << "^{" << (term.first.first / 2) << "}";
            }
          }
          else if (term.first.first != 2) {  // odd Alexander weight
            stream << "^{\\frac{" << term.first.first << "}{2}}";
          }
        }
        if (term.first.second != 0) {
          stream << "q";
          if (term.first.second != 1) {
            stream << "^{" << term.first.second << "}";
          }
        }
        if (term.second == 1 and term.first.first == 0 and term.first.second == 0) {
          stream << "1";
        }
        stream << " + ";
      }
      std::string result = stream.str();
      result.erase(result.end() - 3, result.end());
      return result;
    }
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Poincare_polynomial& poly) {
    os << poly.to_string();
    return os;
  }
  
 private:
  Container terms_;
};

#endif  // POINCARE_POLYNOMIAL_H_