#include <boost/any.hpp>

struct Morse_event_options_default {
  using Parameter_type = boost::any;
  
  template< class T >
  static T parameter_cast(Parameter_type p) {
    return boost::any_cast< T >(p);
  }
};

struct Morse_event_options_int {
  using Parameter_type = int;
  
  template< class T = int >
  static T parameter_cast(Parameter_type p) {
    return p;
  }
};