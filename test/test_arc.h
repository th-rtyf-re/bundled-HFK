// #include "Bordered_algebra/Bordered_algebra.h"
// 
// /* Potentially a Plain-Old-Data type? */
// template< class Forest_options >
// class Arc {
//  public:
//   using Coefficient_value = typename Forest_options::Alg_el;
//   
//   Arc(int source, int target, Coefficient_value value) :
//     source_(source),
//     target_(target),
//     value_(value)
//   { }
//   
//   int source() const {
//     return source_;
//   }
//   
//   int target() const {
//     return target_;
//   }
//   
//   Coefficient_value value() const {
//     return value_;
//   }
//   
//   bool operator==(const Arc& other) const {
//     return (source_ == other.source_ and target_ == other.target_ and value_ == other.value_);
//   }
//   
//  private:
//   int source_;
//   int target_;
//   
//   Coefficient_value value_;
// };
// 
// template< class Forest_options >
// std::ostream& operator<<(std::ostream& os, const Arc< Forest_options >& arc) {
//   os << "(" << arc.source() << "|" << arc.value() << "|" << arc.target() << ")" << std::flush;
//   return os;
// }