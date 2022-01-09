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