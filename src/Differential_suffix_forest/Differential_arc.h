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

#ifndef DIFFERENTIAL_ARC_H_
#define DIFFERENTIAL_ARC_H_

#include <string>
#include <utility>  // pair, swap
#include <vector>

#include "../Utility/Iterators.h"

/* Differential arc.
 * 
 * This implements coefficient bundles for differential suffix trees.
 * 
 * This class contains a member struct, Light_arc, which is an arc that does
 * not know its neighbors. It is a lighter structure, useful for temporary
 * arcs.
 */
template< class Forest_options >
class Differential_arc {
 public:
  using Arc = Differential_arc;  // rename self
  using Idem = typename Forest_options::Idem;
  using Node = typename Forest_options::Node;
  using Arc_iterator = typename Forest_options::Arc_container::iterator;
  using Algebra_element = typename Forest_options::Alg_el;
  
  using Neighbor_container = std::vector< Arc_iterator >;
  using Neighbor_iterator = Double_star_iterator< Neighbor_container >;
  using Neighbor_range = Custom_range< Neighbor_iterator >;
  using Neighbor_const_iterator = Double_star_const_iterator< Neighbor_container >;
  using Neighbor_const_range = Custom_const_range< Neighbor_const_iterator >;
  
  /* Constructors, destructors, access, assignemnt.
   */
  Differential_arc(Node* source, Node* target,
                   Algebra_element algebra_element, int ht) :
    source_(source),
    target_(target),
    algebra_element_(algebra_element),
    ht_(ht),
    // Discrete Morse theory
    to_be_reversed_(false),
    not_reversed_(false),
    gradient(false)
  { }
  
  Node* source() const {
    return source_;
  }
  
  Node* target() const {
    return target_;
  }
  
  Idem source_idem() const {
    return algebra_element_.source_idem();
  }
  
  Idem target_idem() const {
    return algebra_element_.target_idem();
  }
  
  Algebra_element algebra_element() const {
    return algebra_element_;
  }
  
  int ht() const {
    return ht_;
  }
  
  Neighbor_range others_to_source() {
    return Neighbor_range(others_to_source_);
  }
  
  Neighbor_range others_from_source() {
    return Neighbor_range(others_from_source_);
  }
  
  Neighbor_range others_to_target() {
    return Neighbor_range(others_to_target_);
  }
  
  Neighbor_range others_from_target() {
    return Neighbor_range(others_from_target_);
  }
  
  /* This list happens to be called in a const context, so we overload with a
   * const range
   */
  Neighbor_const_range others_to_source() const {
    return Neighbor_const_range(others_to_source_);
  }
  
  void add_source_incoming_arc_it(Arc_iterator other_it) {
    others_to_source_.push_back(other_it);
  }
  
  void add_source_outgoing_arc_it(Arc_iterator other_it) {
    others_from_source_.push_back(other_it);
  }
  
  void add_target_incoming_arc_it(Arc_iterator other_it) {
    others_to_target_.push_back(other_it);
  }
  
  void add_target_outgoing_arc_it(Arc_iterator other_it) {
    others_from_target_.push_back(other_it);
  }
  
  /* General functions.
   */
  
  /* Two arcs are equal if they start at the same node, arrive at the same node,
   * and have the same algebra element (at this point, it can only differ in
   * U weight).
   */
  bool operator==(const Arc& other) const {
    return (source_ == other.source_ and target_ == other.target_
            and algebra_element_ == other.algebra_element_);
  }
  
  void clear_neighbors() {
    others_to_source_.clear();
    others_from_source_.clear();
    others_to_target_.clear();
    others_from_target_.clear();
  }
  
  bool is_closed() const {
    return algebra_element_.is_invertible();
  }
  
  /* Produce the differential cells of the dual D-module. */
  void dualize() {
    algebra_element_.dualize();
    std::swap(source_, target_);
    others_to_source_.swap(others_from_target_);
    others_from_source_.swap(others_to_target_);  // not necessary for our use
  }
  
  std::string algebra_string() const {
    return algebra_element_.to_string();
  }
  
  /* Concatenating arcs involves comparing the target path of the first arc and
   * the source path of the second arc.
   */
  
  /* Check if this arc is compatible with another arc above a certain minimum
   * height. For example, in the case where this height < other height:
   *  
   *   ?  ? |  |
   *   ?  ? |  |                     ?   |   |
   *   .  . .  . <- lower_ht         ?   |   |
   *   \  \ /  /                     .   .   .
   *    \  .  /  <- min_ht            \  |  /
   *     \ | /                         \ | /
   *      \|/                           \|/
   *       .                             .
   *   arc   other           implied case for min_ht = -1
   * 
   * First, we check if the arc source can be extended like the other source.
   * If min_ht = -1, this is all that we do.
   * Next, calculate lower_ht. If lower_ht > min_ht, then we also need to check
   * that the arc target can be extended like the other source. Finally, we
   * check that the arc target and the other source coincide from lower_ht down
   * to min_ht.
   * 
   * We proceed similarly in other case.
   */
  bool compatible_above(const Arc& other, int min_ht = -1) const {
    bool extendable;
    if (ht_ > other.ht_) {
      extendable = other.target_->descendant_like(target_, ht_ - other.ht_).second;
    }
    else {
      extendable = source_->descendant_like(other.source_, other.ht_ - ht_).second;
    }
    if (min_ht < 0) { return extendable; }
    else if (!extendable) { return false; }
    
    Node* node = target_;
    Node* other_node = other.source_;
    
    /* Bring both nodes to the same height by lowering the higher one */
    int lower_ht;
    if (ht_ >= other.ht_) {
      lower_ht = ht_;
      for (; lower_ht > other.ht_; --lower_ht) {
        node = node->parent();
      }
    }
    else {
      lower_ht = other.ht_;
      for (; lower_ht > ht_; --lower_ht) {
        other_node = other_node->parent();
      }
    }
    
    if (lower_ht > min_ht) {
      if (ht_ > other.ht_) {
        extendable = other.source_->descendant_like(target_, ht_ - other.ht_).second;
      }
      else {
        extendable = target_->descendant_like(other.source_, other.ht_ - ht_).second;
      }
    }
    if (!extendable) { return false; }
    
    for (; lower_ht > min_ht; --lower_ht) {
      if (node->type() != other_node->type()) {
        return false;
      }
      node = node->parent();
      other_node = other_node->parent();
    }
    return true;
  }
    
  /* Check that both nodes coincide below as well */
  bool compatible(const Arc& other) const {
    if (algebra_element_.target_idem() != other.algebra_element_.source_idem()) {
      return false;
    }
    
    Node* node = target_;
    Node* other_node = other.source_;
    
    /* Bring both nodes to the same height by lowering the higher one */
    if (ht_ >= other.ht_) {
      int h = ht_;
      for (; h > other.ht_; --h) {
        node = node->parent();
      }
    }
    else {
      int h = other.ht_;
      for (; h > ht_; --h) {
        other_node = other_node->parent();
      }
    }
    
    /* Scan both paths to root */
    while (node->type() == other_node->type()) {
      node = node->parent();
      other_node = other_node->parent();
      if (!node->has_parent()) {  // reached root
        return compatible_above(other);  // check above
      }
    }
    return false;  // exited while loop early
  }
  
  /* A light arc contains the information of a differential arc, but without
   * knowledge of neighbors. This is an optimization to avoid excessive
   * creation and deletion of vectors
   */
  struct Light_arc {
    Light_arc(Node* s, Node* t, Algebra_element a, int h) :
      source(s),
      target(t),
      algebra_element(a),
      ht(h)
    { }
    
    Node* source;
    Node* target;
    
    Algebra_element algebra_element;
    int ht;
  };
  using Light = Light_arc;
  
  /* Concatenate this arc with another arc, provided that they are compatible
   */
  Arc concatenate(const Arc& other) const {
    Node* new_source = source_;
    Node* new_target = other.target_;
    Algebra_element new_alg_el = algebra_element_ * other.algebra_element_;
    int new_ht = ht_;
    if (ht_ > other.ht_) {
      new_target = new_target->descendant_like(target_, ht_ - other.ht_).first;
    }
    else {
      new_source = new_source->descendant_like(other.source_, other.ht_ - ht_).first;
      new_ht = other.ht_;
    }
    return Arc(new_source, new_target, new_alg_el, new_ht);
  }
  
  Arc operator*(const Arc& other) const {
    return this->concatenate(other);
  }
  
  /* Concatenate, but make a light arc */
  Light_arc light_concatenate(const Arc& other) const {
    Node* new_source = source_;
    Node* new_target = other.target_;
    Algebra_element new_alg_el = algebra_element_ * other.algebra_element_;
    int new_ht = ht_;
    if (ht_ > other.ht_) {
      new_target = new_target->descendant_like(target_, ht_ - other.ht_).first;
    }
    else {
      new_source = new_source->descendant_like(other.source_, other.ht_ - ht_).first;
      new_ht = other.ht_;
    }
    return Light_arc(new_source, new_target, new_alg_el, new_ht);
  }
  
  /* Discrete Morse theory.
   */
  bool to_be_reversed() const {
    return to_be_reversed_;
  }
  
  bool not_reversed() const {
    return not_reversed_;
  }
  
  void set_to_be_reversed() {
    to_be_reversed_ = true;
    not_reversed_ = false;
  }
  
  void set_not_reversed() {
    not_reversed_ = true;
    to_be_reversed_ = false;
  }
  
  void reset_flags() {
    not_reversed_ = false;
    to_be_reversed_ = false;
  }
  
  typename std::list< Arc >::iterator gradient_paths_begin() {
    return incoming_partial_gradient_paths_.begin();
  }
  
  typename std::list< Arc >::iterator gradient_paths_end() {
    return incoming_partial_gradient_paths_.end();
  }
  
  std::list< Arc >& gradient_paths() {  // rvalue
    return incoming_partial_gradient_paths_;
  }
  
  void add_incoming_partial_gradient_path(Arc partial_gradient_path_ptr) {
    incoming_partial_gradient_paths_.push_back(partial_gradient_path_ptr);
  }
  
  std::vector< Arc > raised_arcs(int rel_ht) {
    std::vector< Arc > new_arcs;
    if (rel_ht == 0) {  // just returns this arc, shouldn't be used
      new_arcs.emplace_back(source_, target_, algebra_element_, ht_);
    }
    else {
      for (auto source_child_it = source_->children_begin();
                source_child_it != source_->children_end(); ++source_child_it) {
        for (auto target_child_it = target_->children_begin();
                  target_child_it != target_->children_end(); ++target_child_it) {
          if (source_child_it->type() == target_child_it->type()) {
            raised_arcs_aux(rel_ht - 1, new_arcs, &*source_child_it, &*target_child_it, algebra_element_, ht_ + 1);
          }
        }
      }
    }
    return new_arcs;
  }
  
  void raised_arcs_aux(int rel_ht, std::vector< Arc >& new_arcs,
    Node* current_source, Node* current_target, Algebra_element& alg_el, int current_ht) {
    if (rel_ht == 0) {
      new_arcs.emplace_back(current_source, current_target, alg_el, current_ht);
    }
    else {
      for (auto source_child_it = current_source->children_begin();
                source_child_it != current_source->children_end(); ++source_child_it) {
        for (auto target_child_it = current_target->children_begin();
                  target_child_it != current_target->children_end(); ++target_child_it) {
          if (source_child_it->type() == target_child_it->type()) {
            raised_arcs_aux(rel_ht - 1, new_arcs, &*source_child_it, &*target_child_it, alg_el, current_ht + 1);
          }
        }
      }
    }
  }
  
  /* Check if the arc is valid, that is, its source and target nodes can be
   * extended along the same path.
   */
  bool same_paths_exist() const {
    if (!source_->has_children() or !target_->has_children()) {
      return true;
    }
    bool result = false;
    for (auto source_child_it = source_->children_begin();
              source_child_it != source_->children_end(); ++source_child_it) {
      for (auto target_child_it = target_->children_begin();
                target_child_it != target_->children_end(); ++target_child_it) {
        if (source_child_it->type() == target_child_it->type()) {
          result = same_paths_exist(&*source_child_it, &*target_child_it) or result;
        }
      }
    }
    return true;
  }
  
  bool same_paths_exist(Node* current_source, Node* current_target) const {
    if (!source_->has_parent() or !target_->has_parent()) {
      return true;
    }
    bool result = false;
    for (auto source_child_it = current_source->children_begin();
              source_child_it != current_source->children_end(); ++source_child_it) {
      for (auto target_child_it = current_target->children_begin();
                target_child_it != current_target->children_end(); ++target_child_it) {
        if (source_child_it->type() == target_child_it->type()) {
          result = same_paths_exist(&*source_child_it, &*target_child_it) or result;
        }
      }
    }
    return true;
  }
  
  void reverse() {  // not dualize
    algebra_element_.dualize();
    std::swap(source_, target_);
    others_to_source_.swap(others_to_target_);
    others_from_source_.swap(others_from_target_);
  }
  
  bool gradient;
  
  void debug_out(std::string prefix = "[da]") const {
    std::cout << prefix << " arc: " << source_->chain << std::flush;
    std::cout << " (" << source_->n_gen() << std::flush;
    std::cout << ") " << target_->chain << std::flush;
    std::cout << " (" << target_->n_gen() << std::flush;
    std::cout << ") " << algebra_element_.to_string() << std::endl;
  }
  
 private:
  /* Member declarations.
   */
  Node* source_;  // note that these are pointers
  Node* target_;
  
  Algebra_element algebra_element_;
  
  int ht_;  // height
  
  /* I'm not sure if this is the best way to keep track of arc neighbors,
   * because this needs to be handled explicitly. As a consequence, the D-module
   * can be in an invalid state (arcs not connected)
   *
   */
  Neighbor_container others_to_source_;
  Neighbor_container others_from_source_;
  Neighbor_container others_to_target_;
  Neighbor_container others_from_target_;
  
  /* Discrete Morse theory */
  bool to_be_reversed_;
  bool not_reversed_;
  
  std::list< Arc > incoming_partial_gradient_paths_;
};

#endif  // DIFFERENTIAL_ARC_H_