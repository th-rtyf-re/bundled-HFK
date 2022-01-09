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

#ifndef TWO_BIT_SET_H_
#define TWO_BIT_SET_H_

#include <bitset>

/* Two bit set.
 *
 * Auxiliary class for positive crossings. Current convention is big-endian.
 */
template< std::size_t size >
class Two_bit_set {
 public:
  using Storage = std::bitset< 2 * size >;
  
  Two_bit_set() { }
  
  int operator[](int i) const {
    return table_[2 * i] + 2 * table_[2 * i + 1];
  }
  
  void set_value(int i, int value) {
    table_[2 * i] = value & 1;
    table_[2 * i + 1] = (value >> 1) & 1;
  }
  
  std::string to_string() const {
    return table_.to_string();
  }
  
 private:
  Storage table_;
};

#endif  // TWO_BIT_SET_H_