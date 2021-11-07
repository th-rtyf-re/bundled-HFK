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

#ifndef IDEMPOTENT_H_
#define IDEMPOTENT_H_

#include <cstdint>  // int_fast32_t, I think
#include <iostream>
#include <initializer_list>
#include <string>
#include <vector>
#include <utility>  // swap

/* Idempotents.
 * 
 * This is our implementation of idempotent elements.
 * 
 * The length of an idempotent is equal to the number of strands + 1. If this
 * number is at most 32, then we use 32-bit integers. Otherwise, we use a
 * dynamic object, such as std::vector<bool> or boost::dynamic_bitset.
 * 
 * Subject to change.
 */

/* Idempotent_short.
 * 
 * The type of the idempotent is the fastest integer with at least 32 bits. For
 * my computer, these are just integers, which means I didn't really check
 * integer types. This may cause problems for other users.
 * 
 * Methods should correspond to those of Idempotent_long, which takes from
 * another class.
 * 
 * Behavior is undefined if initialized with an idempotent of length greater
 * than 32.
 */
using Idempotent_short_type = int_fast32_t;
class Idempotent_short {
 public:
  Idempotent_short() : data_(0) { }
  
  /* Convert string to integer...
   * 
   * In practice, this method should be called only when defining the
   * idempotents of aligned maxima or minima.
   */
  Idempotent_short(std::string result) {
    actual_size_ = result.size();
    data_ = 0;
    
    /* Create idempotent by iterating backwards over input string */
    std::string::reverse_iterator rit;
    for (rit = result.rbegin(); rit != result.rend(); ++rit) {
      if (*rit == '0') {
        data_ = data_ << 1;
      }
      else {
        data_ = (data_ << 1) + 1;
      }
    }
  }
  
  bool operator ==(const Idempotent_short& other) const {
    return (data_ == other.data_);
  }
  
  bool operator !=(const Idempotent_short& other) const {
    return (data_ != other.data_);
  }
  
  /* Compare two idempotents. This does not act the same way as Idempotent_long.
   * It's doing lexicographic order of idempotents, in the reverse order, and
   * then reversing the entire order.
   */
  bool operator <(const Idempotent_short& other) const {
    return (data_ < other.data_);
  }
  
  bool operator[](int i) const {
    return (data_ & (1 << i));
  }
  
  size_t size() const {
    return actual_size_;
  }
  
  void flip(int i) {
    data_ = data_ ^ (1 << i);
  }
  
  /* Insert boolean values right before pos */
  void insert(int pos, std::initializer_list< bool > ilist) {
    int power = 1 << pos;
    int q = data_ / power;  // quotient
    int r = data_ % power;  // remainder
    int n_insert = ilist.size();
    for (bool bit : ilist) {
      r += bit << pos;
      ++pos;
    }
    q <<= pos;
    actual_size_ += n_insert;
    data_ = q + r;
  }
  
  /* Erase elements from pos to pos + n_erase - 1 */
  void erase(int pos, int n_erase) {
    int power = 1 << pos;
    int q = data_ / power;
    int r = data_ % power;
    q >>= n_erase;
    q <<= pos;
    actual_size_ -= n_erase;
    data_ = q + r;
  }
  
  bool too_far_from(Idempotent_short other) const {
    /* current_difference only takes values -1, 0, 1 */
    int current_difference = 0;
    for (int i = 0; i < actual_size_; ++i) {
      /* This idempotent is in the lead, but the other one isn't catching up
       */
      if (current_difference > 0 and !other[i]) {
        return true;
      }
      /* Other idempotent is in the lead, but this one isn't catching up
       */
      else if (current_difference < 0 and !Idempotent_short::operator[](i)) {
        return true;
      }
      current_difference = current_difference + Idempotent_short::operator[](i)
                                              - other[i];
    }
    return false;
  }
  
  void swap(Idempotent_short& other) {
    std::swap(data_, other.data_);
    std::swap(actual_size_, other.actual_size_);
  }
  
  /* Convert the idempotent to a string, for display and debugging purposes.
   */
  std::string to_string() const {
    std::string result = "";
    for (int i = 0; i < actual_size_; ++i) {
      // std::cout << "i = " << i << ", actual size = " << actual_size_ << std::endl;
      if (Idempotent_short::operator[](i)) {
        result = result + "1";
      }
      else {
        result = result + "0";
      }
    }
    return result;
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Idempotent_short& idem) {
    os << idem.to_string();
    return os;
  }
  
 private:
  Idempotent_short_type data_;
  
  size_t actual_size_;
};


/* Idempotent_long.
 * 
 * Currently some sort of wrapper for std::vector<bool>, although it might be
 * better to use boost::dynamic_bitset.
 */
template< class Bit_container >
class Idempotent_long {
 public:
  Idempotent_long() { }
  
  /* Convert string to integer...
   * 
   * In practice, this method should be called only when defining the
   * idempotents of aligned maxima or minima.
   */
  Idempotent_long(std::string result) {
    for (int i = 0; i < result.size(); ++i) {
      if (result[i] == '0') {
        data_.push_back(false);
      }
      else {
        data_.push_back(true);
      }
    }
  }
  
  bool operator ==(const Idempotent_long& other) const {
    return (data_ == other.data_);
  }
  
  bool operator !=(const Idempotent_long& other) const {
    return (data_ != other.data_);
  }
  
  bool operator <(const Idempotent_long& other) const {
    return (data_ < other.data_);
  }
  
  bool operator[](int i) const {
    return data_[i];
  }
  
  size_t size() const {
    return data_.size();
  }
  
  /* flip an individual bit
   * Note that there may exist a function flip() which takes no arguments and
   * flips every value.
   */
  void flip(int i) {
    data_[i] = ~data_[i];
  }
  
  /* Insert boolean values right before pos */
  void insert(int pos, std::initializer_list< bool > ilist) {
    data_.insert(data_.begin(), ilist);
  }
  
  void erase(int pos, int n_erase) {
    data_.erase(data_.begin() + pos, data_.begin() + pos + n_erase - 1);
  }
  
  bool too_far_from(Idempotent_long other) const {
    /* current_difference only takes values -1, 0, 1 */
    int current_difference = 0;
    for (int i = 0; i < data_.size(); ++i) {
      /* First idempotent is in the lead, but the second one isn't catching up
       */
      if (current_difference > 0 and !other.data_[i]) {
        return true;
      }
      /* Second idempotent is in the lead, but the first one isn't catching up
       */
      else if (current_difference < 0 and !data_[i]) {
        return true;
      }
      current_difference = current_difference + data_[i]
                                              - other.data_[i];
    }
    return false;
  }
  
  void swap(Idempotent_long& other) {
    data_.swap(other.data_);
  }
  
  /* Convert the idempotent to a string, for display and debugging purposes.
   */  
  std::string to_string() const {
    std::string result = "";
    for (int i = 0; i < data_.size(); ++i) {
      if (data_[i]) {
        result = result + "1";
      }
      else {
        result = result + "0";
      }
    }
    return result;
  }
  
  friend std::ostream& operator<<(std::ostream& os, const Idempotent_long< Bit_container >& idem) {
    os << idem.to_string();
    return os;
  }
  
 private:
  Bit_container data_;
};


#endif  // IDEMPOTENT_H_