#ifndef TEST_MORSE_EVENT_H_
#define TEST_MORSE_EVENT_H_

#include <iostream>
#include <string>
#include <utility>  // pair
#include <vector>

#include <boost/type_erasure/any.hpp>
#include <boost/type_erasure/member.hpp>  // BOOST_TYPE_ERASURE_MEMBER
#include <boost/type_erasure/operators.hpp>  // ostreamable

#include <boost/any.hpp>

/* Some metaprogramming nonsense to define Morse events.
 */

BOOST_TYPE_ERASURE_MEMBER(lower_matchings)
BOOST_TYPE_ERASURE_MEMBER(upper_orientations)
BOOST_TYPE_ERASURE_MEMBER(update_margins)
BOOST_TYPE_ERASURE_MEMBER(to_string)

struct Topological_methods :
  boost::mpl::vector<
    has_lower_matchings<
      std::vector< int >(std::vector< int >),
      const boost::type_erasure::_self
    >,
    has_upper_orientations<
      std::vector< bool >(
        std::vector< bool >,
        const std::vector< int >&
      ),
      const boost::type_erasure::_self
    >,
    has_update_margins<
      std::pair< int, int >(std::pair< int, int >),
      const boost::type_erasure::_self
    >,
    has_to_string<
      std::string(const std::pair< int, int >, const std::pair< int, int >&),
      const boost::type_erasure::_self
    >
  >
{ };

BOOST_TYPE_ERASURE_MEMBER(get_weights)
BOOST_TYPE_ERASURE_MEMBER(get_labels)
BOOST_TYPE_ERASURE_MEMBER(tensor_generators)
BOOST_TYPE_ERASURE_MEMBER(tensor_coefficients)

template< class D_module >
struct Algebraic_methods :
  boost::mpl::vector<
    has_get_weights<
      std::vector< typename D_module::Weights >(
        const typename D_module::Bordered_algebra&,
        const typename D_module::Bordered_algebra&
      ),
      const boost::type_erasure::_self
    >,
    has_get_labels<
      std::vector< std::string >(
        const typename D_module::Bordered_algebra&,
        const typename D_module::Bordered_algebra&
      ),
      const boost::type_erasure::_self
    >,
    has_tensor_generators<
      D_module(
        D_module&,
        const D_module&,
        const typename D_module::Bordered_algebra&,
        const typename D_module::Bordered_algebra&
      ),
      const boost::type_erasure::_self
    >,
    has_tensor_coefficients<
      D_module(
        D_module&,
        const D_module&,
        const typename D_module::Bordered_algebra&,
        const typename D_module::Bordered_algebra&
      ),
      const boost::type_erasure::_self
    >
  >
{ };


template< class D_module, class Morse_event_options >
using Morse_event_base = boost::type_erasure::any<
  boost::mpl::vector<
    boost::type_erasure::copy_constructible<>,
    boost::type_erasure::constructible<
      boost::type_erasure::_self(
        const std::vector< typename Morse_event_options::Parameter_type >&
      )
    >,
    Topological_methods,
    Algebraic_methods< D_module >,
#ifdef VERBOSE
    boost::type_erasure::ostreamable<>,
#endif  // VERBOSE
    boost::type_erasure::relaxed
  >
>;

template< class D_module_template, class Morse_event_options >
struct Morse_event : Morse_event_base< D_module_template, Morse_event_options > {
  using Base = Morse_event_base< D_module_template, Morse_event_options >;
  
  using Base::Base;
  using D_module = D_module_template;
};

#endif  // TEST_MORSE_EVENT_H_