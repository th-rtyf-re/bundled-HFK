#ifndef TEST_POSITIVE_CROSSING_H_
#define TEST_POSITIVE_CROSSING_H_

#include <bitset>
#include <cstdlib>  // abs
#include <string>
#include <tuple>
#include <utility>  // pair, swap
#include <vector>

#ifdef VERBOSE
#include <iostream>
#endif  // VERBOSE

#include <boost/any.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/for_each.hpp>

/* Two bit set.
 *
 * Class with two bits. Current convention is big-endian.
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

/* DA bimodule for a crossing.
 * 
 * Adapted from ComputeHFK2/Crossing.cpp.
 */
template< class D_module, class Morse_event_options >
class Positive_crossing {
 public:
  using Idem = typename D_module::Idem;
  using Gen_type = typename D_module::Gen_type;
  using Algebra = typename D_module::Bordered_algebra;
  using Coef_bundle = typename D_module::Coef_bundle;
  using Weights = typename D_module::Weights;
  
  using Table = Two_bit_set< 128 >;
  
  Positive_crossing(const std::vector< typename Morse_event_options::Parameter_type >& args) :
    position_(args.empty() ? 0 : Morse_event_options::template parameter_cast< int >(args[0]))
  { }
  
  /* Topological methods */
  
  std::vector< int > lower_matchings(std::vector< int > matchings) const {
    std::swap(matchings[position_], matchings[position_ + 1]);
    std::swap(
      matchings[matchings[position_]],
      matchings[matchings[position_ + 1]]
    );
    return matchings;
  }
  
  std::vector< bool > upper_orientations(
    std::vector< bool > orientations,
    const std::vector< int >&
  ) const {
    std::swap(orientations[position_], orientations[position_ + 1]);
    return orientations;
  }
  
  std::pair< int, int > update_margins(std::pair< int, int > margins) const {
    return margins;
  }
  
  /* Return the LaTeX KnotDiagram2ASCII string for the knot slice.
   */
  std::string to_string(
    const std::pair< int, int >& margins,
    const std::pair< int, int >& n_strands,
    std::string symbol = "+"
  ) const {
    return std::string(margins.first, '0')
      + std::string(position_, '.')
      + symbol
      + std::string(n_strands.second - position_ - 2, '.')
      + std::string(margins.second, '0');
  }
  
  
  /* Algebraic methods */
  
  std::vector< Weights > get_weights(
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    std::vector< Weights > weights;
    for (auto type : {N, E, S, W}) {
      weights.push_back(
        alexander_maslov_weights_(type, upper_algebra, lower_algebra)
      );
    }
    return weights;
  }
  
  std::vector< std::string > get_labels(
    const Algebra& upper_algebra,
    const Algebra& lower_algebra,
    std::string symbols = "NESW"
  ) const {
    std::vector< std::string > labels;
    for (auto type : {N, E, S, W}) {
      labels.push_back(
        std::string(1, symbols[type]) + "_{" + std::to_string(position_) + "}"
      );
    }
    return labels;
  }
  
  D_module tensor_generators(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra&,
    const Algebra&
  ) const {
    delta_0_(new_d_module, old_d_module);
    return new_d_module;
  }
  
  D_module& tensor_coefficients(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    delta_1_(new_d_module, old_d_module, upper_algebra, lower_algebra);
    delta_2_(new_d_module, old_d_module);
    delta_3_(new_d_module, old_d_module);
    return new_d_module;
  }
  
#ifdef VERBOSE
  template< class, class >
  friend class Negative_crossing;
  
  friend std::ostream& operator<<(
    std::ostream& os,
    const Positive_crossing& morse_event
  ) {
    os << "positive crossing at position " << morse_event.position_;
    return os;
  }
#endif  // VERBOSE
  
 private:
  enum {
    N = 0,
    E,  // 1
    S,  // 2
    W,  // 3
    null_NEW = 2,
    null_S = 0,
    null = 0
  };
  
  Weights alexander_maslov_weights_(
    const int generator_type,
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    bool left_ori = upper_algebra.orientations[position_];
    bool right_ori = upper_algebra.orientations[position_ + 1];
    int a = 0;  // double Alexander
    int m = 0;  // Maslov
    
    /* Case study based on the direction in which both strands are pointing */
    if (left_ori) {
      if (right_ori) {  // positive, pointing N
        if (generator_type == N) {
          a = -1;
          m = -1;
        }
        else if (generator_type == S) {
          a = 1;
        }
      }
      else {  // negative, pointing W
        if (generator_type == E) {
          a = -1;
        }
        else if (generator_type == W) {
          a = 1;
          m = 1;
        }
      }
    }
    else {
      if (right_ori) {  // negative, pointing E
        if (generator_type == E) {
          a = 1;
          m = 1;
        }
        else if (generator_type == W) {
          a = -1;
        }
      }
      else {  // positive, pointing S
        if (generator_type == N) {
          a = 1;
        }
        else if (generator_type == S) {
          a = -1;
          m = -1;
        }
      }
    }
    
    return {a, m};
  }
  
  /* \delta functions */
  
  /* \delta_0 and \delta_1.
   * 
   * \delta_0 refers to the tensor product part of the box tensor product, that
   * is, the addition of nodes and edges.
   */
  void delta_0_(D_module& new_d_module, const D_module& old_d_module) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      const Idem& old_idem = old_d_module.idem(gen_handle);
      
      if (old_idem[position_ + 1]) {  // N marking
        new_d_module.add_gen_bundle(old_idem, N, gen_handle);
      }
      else {  // E, S, W markings
        new_d_module.add_gen_bundle(old_idem, S, gen_handle);
        
        for (int we : {0, 1}) {  // we = 0 is W, we = 1 is E
          if (!old_idem[position_ + 2 * we]) {  // can't add generator
            continue;
          }
          Gen_type we_marking = we ? E : W;
          Idem we_idem = extend_(old_idem, we_marking);
          new_d_module.add_gen_bundle(we_idem, we_marking, gen_handle);
        }  // W, E
      }
    }
  }  // delta_0_and_1_
  
  void delta_1_(
    D_module& new_d_module,
    const D_module& old_d_module,
    const Algebra& upper_algebra,
    const Algebra& lower_algebra
  ) const {
    for (const auto& gen_handle : old_d_module.gen_bundle_handles()) {
      const Idem& old_idem = old_d_module.idem(gen_handle);
      
      if (!old_idem[position_ + 1]) {
        for (int we : {0, 1}) {
          if (!old_idem[position_ + 2 * we]) { continue; }
          Gen_type we_marking = we ? E : W;
          Idem we_idem = extend_(old_idem, we_marking);
          auto alg_el = new_d_module.alg_el(we_idem, old_idem);
          new_d_module.add_coef_bundle(alg_el, we_marking, S, old_idem);
          
          if (upper_algebra.matchings[position_] != position_ + 1) {  // curved
            std::vector< int > U_curved(upper_algebra.n_strands, 0);
            U_curved[lower_algebra.matchings[position_ + we]] = 1;
            auto alg_el = new_d_module.alg_el(old_idem, we_idem, U_curved);
            new_d_module.add_coef_bundle(alg_el, S, we_marking, old_idem);
          }
        }
      }
    }
  }
  
  /* \delta_2.
   */
  void delta_2_(D_module& new_d_module, const D_module& old_d_module) const {
    for (const auto& coef : old_d_module.coef_bundles()) {
      /* Calculate preliminary information */
      int a1, a2, u1, u2;
      std::tie(a1, a2, u1, u2) = get_local_weights_(coef, old_d_module);
      int pre_hash_index = pre_hash_index_(a1, a2, u1, u2);
      
      for (const Gen_type front_marking : {N, W, S, E}) {
        if (!extendable_(old_d_module.target_idem(coef), front_marking)) {
          continue;
        }
        const Gen_type back_marking =
          positive_look_back_[pre_hash_index + front_marking];
        if (!extendable_(old_d_module.source_idem(coef), back_marking)) {
          continue;
        }
        
        const Idem new_source_idem =
          extend_(old_d_module.source_idem(coef), back_marking);
        const Idem new_target_idem =
          extend_(old_d_module.target_idem(coef), front_marking);
        if (new_source_idem.too_far_from(new_target_idem)) {
          continue;  // incompatible idems
        }
        
        /* Calculate new algebra element */
        std::vector<int> U_weights = old_d_module.U_weights(coef);
        int v1 = 2 * u1 + std::abs(a1);
        int v2 = 2 * u2 + std::abs(a2);
        if (back_marking == E) { --v2; }
        else if (back_marking == W) { --v1; } 
        if (front_marking == E) { ++v2; }
        else if (front_marking == W) { ++v1; }
        U_weights[position_] = v2 / 2;
        U_weights[position_ + 1] = v1 / 2;
        auto alg_el =
          new_d_module.alg_el(new_source_idem, new_target_idem, U_weights);
        new_d_module.add_coef_bundle(
          alg_el,
          back_marking,
          front_marking,
          coef,
          old_d_module
        );
      }
    }
  }  // delta_2_
  
  /* \delta_3.
   * 
   * Pre-condition: we have calculated the neighbors of differential arcs.
   */
  void delta_3_(D_module& new_d_module, const D_module& old_d_module) const {
    for (const auto& front_coef : old_d_module.coef_bundles()) {
      for (const auto& back_coef : old_d_module.others_to_source(front_coef)) {
        for (const Gen_type front_marking : {N, E, S, W}) {
          if (
            !extendable_(old_d_module.target_idem(front_coef), front_marking)
          ) {
            continue;
          }
          // check if I need this
          //if (!extendable_(back_alg_el.source_idem(), S)) { continue; }
          if (
            !coef_exists_(back_coef, front_coef, front_marking, old_d_module)
          ) {
            continue;
          }
          
          // back_marking is S
          const Idem& new_source_idem = old_d_module.source_idem(back_coef);
          const Idem& new_target_idem =
            extend_(old_d_module.target_idem(front_coef), front_marking);
          if (new_source_idem.too_far_from(new_target_idem)) {
            continue;  // algebra element is null
          }
          
          /* Calculate new algebra element */
          auto concat_coef = old_d_module.concatenate(back_coef, front_coef);
          auto new_U_weights = old_d_module.U_weights(concat_coef);
          int a1, a2, u1, u2, b1, b2, v1, v2;
          std::tie(a1, a2, u1, u2) =
            get_local_weights_(back_coef, old_d_module);
          std::tie(b1, b2, v1, v2) =
            get_local_weights_(front_coef, old_d_module);
          int w1 = 2 * u1 + 2 * v1 + std::abs(a1) + std::abs(b1) - 1;
          int w2 = 2 * u2 + 2 * v2 + std::abs(a2) + std::abs(b2) - 1;
          if (front_marking == E) { ++w2; } // extra weight for L_2
          else if (front_marking == W) { ++w1; } // extra weight for R_1
          new_U_weights[position_] = w2 / 2;
          new_U_weights[position_ + 1] = w1 / 2;
          /* Make and add new differential arc */
          auto alg_el = new_d_module.alg_el(
            new_source_idem,
            new_target_idem,
            new_U_weights
          );
          new_d_module.add_coef_bundle(
            alg_el,
            S,
            front_marking,
            concat_coef,
            old_d_module
          );
        }
      }
    }
  }  // delta_3_
  
  /* Auxiliary functions */
    
  /* Adapted from ComputeHFKv2/Utility.cpp, LeftRight.
   * For \delta_2 and \delta_3. 
   * Get local LR weights and U weights. This is because the DA-bimodule for a
   * crossing depends on these local weights.
   * 
   * I factored this out, but this means it gets called more than strictly
   * necessary.
   */
  std::tuple< int, int, int, int > get_local_weights_(
    Coef_bundle coef,
    const D_module& old_d_module
  ) const {
    const Idem& source_idem = old_d_module.source_idem(coef);
    const Idem& target_idem = old_d_module.target_idem(coef);
    int a1 = 0;  // will only take values in [-1, 1]
    int a2 = 0;  // will only take values in [-1, 1]
    int i = 0;
    for ( ; i <= position_; ++i) {
      a1 = a1 + source_idem[i] - target_idem[i];
    }
    a2 = a1 + source_idem[i] - target_idem[i];  // i = position_ + 1
    return std::make_tuple(
      a1,
      a2,
      old_d_module.U_weights(coef)[position_],
      old_d_module.U_weights(coef)[position_ + 1]
    );
  }  // get_local_weights_
  
  /* For \delta_3
   *
   * We follow [OzsvathSzabo2018, Lemma 5.5]. Given a front marking Y, we
   * compute the mid marking I(b, Y), the front marking I(a, I(b, Y)), and the
   * product marking I(ab, Y).
   */
  bool coef_exists_(
    const Coef_bundle& back_coef,
    const Coef_bundle& front_coef,
    const Gen_type front_marking,
    const D_module& old_d_module
  ) const {
    if (front_marking == S) { return false; }  // front cannot be S
    
    int hash_index;
    int a1, a2, u1, u2, b1, b2, v1, v2;  // local weights
    std::tie(a1, a2, u1, u2) = get_local_weights_(back_coef, old_d_module);
    std::tie(b1, b2, v1, v2) = get_local_weights_(front_coef, old_d_module);
    
    hash_index = hash_index_(b1, b2, v1, v2, front_marking);
    Gen_type mid_marking = positive_look_back_[hash_index];
    if (mid_marking == null_NEW) { return false; }  // check that mid isn't null
    
    hash_index = hash_index_(a1, a2, u1, u2, mid_marking);
    Gen_type back_marking = positive_look_back_[hash_index];
    
    hash_index = hash_index_(a1 + b1, a2 + b2,
      u1 + v1 + (std::abs(a1) + std::abs(b1)) / 2,
      u2 + v2 + (std::abs(a2) + std::abs(b2)) / 2, front_marking);
    Gen_type product_marking = positive_look_back_[hash_index];
    // check that back is not equal to product
    if (back_marking == product_marking) { return false; }
    
    /* If product marking is null, we either have (R_1, R_2U_2^t) for E
     * or (L_2, L_1U_1^n) for W. We exclude all other cases.
     * TO DO: Understand why this is here.
     */
    if (product_marking == null_NEW and front_marking == E
      and !(a1 == 1 and a2 == 0 and u1 == 0 and u2 == 0
        and b1 == 0 and b2 == 1 and v1 == 0)) { return false; }
    else if (product_marking == null_NEW and front_marking == W
      and !(a1 == 0 and a2 == -1 and u1 == 0 and u2 == 0
        and b1 == -1 and b2 == 0 and v2 == 0)) { return false; }
    return true;
  }  // back_marking_exists_
    
  /* Taken from ComputeHFKv2/Crossing.cpp, Extendable */
  bool extendable_(const Idem& idem, const Gen_type marking) const {
    switch (marking) {
      case N: return idem[position_ + 1];
      case E: return (!idem[position_ + 1] and idem[position_ + 2]);
      case S: return !idem[position_ + 1];
      case W: return (!idem[position_ + 1] and idem[position_]);
      default: return false;
    }
  }

  /* Taken from ComputeHFKv2/Crossing.cpp, Extend */
  Idem extend_(Idem idem, const Gen_type marking) const {
    if (marking == E) {
      idem.flip(position_ + 1);
      idem.flip(position_ + 2);
    }
    else if (marking == W) {
      idem.flip(position_);
      idem.flip(position_ + 1);
    }
    return idem;  // including cases N, S
  }
  
  /* Static methods for reasons */
  
  static Table make_positive_look_back_() {
    Table table;
    boost::mpl::for_each< boost::mpl::vector_c< int, -1, 0, 1 > >(
      [&table](int a1) {
        boost::mpl::for_each< boost::mpl::vector_c< int, -1, 0, 1 > >(
          [=, &table](int a2) {
            boost::mpl::for_each< boost::mpl::vector_c< int, -1, 0, 1 > >(
              [=, &table](int u1) {
                boost::mpl::for_each< boost::mpl::vector_c< int, N, E, S, W > >(
                  [=, &table](int x) {
                    table.set_value(
                      hash_index_(a1, a2, u1, 0, x),
                      D2PHSEM_(a1, a2, u1, 0, x)
                    );
                  }
                );
              }
            );
          }
        );
      }
    );
    return table;
  }
   
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
  static int hash_index_(int a1, int a2, int u1, int u2, int marking) {
    return pre_hash_index_(a1, a2, u1, u2) + (marking << 0);
  }
  
  static int pre_hash_index_(int a1, int a2, int u1, int u2) {
    int i_local_LR = 4 + a1 - (3 * a2);
    // mod 8. Not strictly necessary after creation of the table
    i_local_LR = i_local_LR & 7;
    
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
  static int D2PHSEM_(int a1, int a2, int u1, int u2, int marking) {
    if (marking == N) {
      if      (             a1 == a2             ) return N;  // 1, L_1L_2, R_2R_1
      else if (a1 ==  1 and a2 ==  0 and u1 <  u2) return E;  // R_1U_2
      else if (a1 ==  1 and a2 ==  0 and u1 >= u2) return W;  // R_1
      else if (a1 ==  0 and a2 == -1 and u1 <= u2) return E;  // L_2
      else if (a1 ==  0 and a2 == -1 and u1 >  u2) return W;  // L_2U_1
      return null_NEW;
    }
    
    else if (marking == E) {
      if      (a1 ==  0 and a2 ==  0 and u1 <= u2) return E;  // 1
      else if (a1 ==  0 and a2 ==  0 and u1 >  u2) return W;  // U1
      else if (a1 == -1 and a2 ==  0             ) return N;  // L_1
      else if (a1 ==  0 and a2 ==  1             ) return N;  // R_2
      return null_NEW;
    }
    
    else if (marking == W) {
      if      (a1 ==  0 and a2 ==  0 and u1 <  u2) return E;  // U2
      else if (a1 ==  0 and a2 ==  0 and u1 >= u2) return W;  // 1
      else if (a1 == -1 and a2 ==  0             ) return N;  // L_1
      else if (a1 ==  0 and a2 ==  1             ) return N;  // R_2
      return null_NEW;
    }
    
    else if (marking == S) {
      if      (             u1 == u2             ) return S;  // 1
      return null_S;
    }
    
    return null;  // this should never happen
  }
  
  int position_;
  
  static const Table positive_look_back_;
};

template< class D_module, class Morse_event_options >
const typename Positive_crossing< D_module, Morse_event_options >::Table
  Positive_crossing< D_module, Morse_event_options >::positive_look_back_ =
    Positive_crossing< D_module, Morse_event_options >::make_positive_look_back_();

#endif  // TEST_POSITIVE_CROSSING_H_