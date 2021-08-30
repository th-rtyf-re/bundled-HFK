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

#ifndef EVENT_CROSSING_HELPER_H_
#define EVENT_CROSSING_HELPER_H_

#include <bitset>

/* Slice crossing helper header.
 *
 * This header file contains some generator types. We also construct the
 * lookup table for \delta_2, at compile time.
 *
 * Technically, we could write everything explicitly, but that requires typing
 * a lot...
 */

using Generator_type = unsigned char;  // = Gen_type in Differential_suffix_forest.h
enum : Generator_type {
  north = 0,
  east,   // 1
  south,  // 2
  west,   // 3
  null_NEW = 2,
  null_S = 0,
  null = 0
};

/* Two bit table.
 *
 * Class with two bits. Current convention is big-endian.
 */
template< std::size_t size >
class Two_bit_table {
 public:
  typedef std::bitset< 2 * size > Two_bit_set;
  
  Two_bit_table() { }
  
  int operator[](int i) const {
    return table_[2 * i] + 2 * table_[2 * i + 1];
  }
  
  void update_value(int i, int value) {
    table_[2 * i] = value & 1;
    table_[2 * i + 1] = (value >> 1) & 1;
  }
  
  std::string to_string() const {
    return table_.to_string();
  }
  
 private:
  Two_bit_set table_;
};

template< class Table_with_128_entries_of_at_least_two_bits >
struct Delta_2_lookup_tables {
  typedef Table_with_128_entries_of_at_least_two_bits Table;
  
  Delta_2_lookup_tables() {  // constructor populates the lookup table
    for (auto a1 : {-1, 0, 1}) {
      for (auto a2 : {-1, 0, 1}) {
        for (auto u1 : {-1, 0, 1}) {
          for (auto x : {north, east, south, west}) {
            positive_look_back.update_value(hash_index(a1, a2, u1, 0, x), D2PHSEM(a1, a2, u1, 0, x));
            //negative_look_fore.update_value(hash_index(a1, a2, u1, 0, x), D2PHSEM(-a1, -a2, u1, 0, x));
          }
        }
      }
    }
  }  // constructor
  
  /* Public members*/
  Table positive_look_back;
  
  //Table negative_look_fore;
    
  /* Auxiliary functions */
  
  /* Make a hash index for the data (local LR weight, local U weight, marking).
   *
   * Each index is unique, except when (a1, a2) = (1, -1) or (-1, 1), in which
   * case they share an index
   *
   * a1 a2 i_local_LR
   *  0  0     4
   *  0  1     1
   *  0 -1     7
   *  1  0     5
   *  1  1     3
   *  1 -1     0
   * -1  0     3
   * -1  1     0
   * -1 -1     6
   *
   * The result is an integer between 0 and 127.
   */
  int hash_index(int a1, int a2, int u1, int u2, Generator_type marking) const {
    return pre_hash_index(a1, a2, u1, u2) + (marking << 0);
  }
  
  int pre_hash_index(int a1, int a2, int u1, int u2) const {
    int i_local_LR = 4 + a1 - (3 * a2);  // nine values, with two 'impossible' ones: 0 and 8.
    i_local_LR = i_local_LR & 7;  // mod 8. Not strictly necessary after creation of the table
    
    int i_U_diff;  // sign of u1 - u2, with 3 = negative
    if (u1 > u2) i_U_diff = 1;
    else if (u1 < u2) i_U_diff = 3;
    else i_U_diff = 0;
    
    return (i_local_LR << 4) + (i_U_diff << 2);
  }
  
  /* \delta_2 positive horizontal source edge marking.
   *
   * Adapted from ComputeHFKv2/Crossing.cpp, LookBack.
   *
   * Given the marking of a horizontal target edge with local LR weight and U
   * weight, return the marking of the horizontal source edge.
   */
  Generator_type D2PHSEM(int a1, int a2, int u1, int u2, Generator_type marking) const {
    if (marking == north) {
      if      (             a1 == a2             ) return north;  // 1, L_1L_2, R_2R_1
      else if (a1 ==  1 and a2 ==  0 and u1 <  u2) return east;   // R_1U_2
      else if (a1 ==  1 and a2 ==  0 and u1 >= u2) return west;   // R_1
      else if (a1 ==  0 and a2 == -1 and u1 <= u2) return east;   // L_2
      else if (a1 ==  0 and a2 == -1 and u1 >  u2) return west;   // L_2U_1
      return null_NEW;
    }
    
    else if (marking == east) {
      if      (a1 ==  0 and a2 ==  0 and u1 <= u2) return east;   // 1
      else if (a1 ==  0 and a2 ==  0 and u1 >  u2) return west;   // U1
      else if (a1 == -1 and a2 ==  0             ) return north;  // L_1
      else if (a1 ==  0 and a2 ==  1             ) return north;  // R_2
      return null_NEW;
    }
    
    else if (marking == west) {
      if      (a1 ==  0 and a2 ==  0 and u1 <  u2) return east;   // U2
      else if (a1 ==  0 and a2 ==  0 and u1 >= u2) return west;   // 1
      else if (a1 == -1 and a2 ==  0             ) return north;  // L_1
      else if (a1 ==  0 and a2 ==  1             ) return north;  // R_2
      return null_NEW;
    }
    
    else if (marking == south) {
      if      (             u1 == u2             ) return south;  // 1
      return null_S;
    }
    
    return null;  // this should never happen
  }  
};

#endif  // EVENT_CROSSING_HELPER_H_