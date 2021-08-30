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

#ifndef DIFFERENTIAL_SUFFIX_FOREST_NODE_
#define DIFFERENTIAL_SUFFIX_FOREST_NODE_

#include <list>
#include <string>
#include <tuple>
#include <utility>  // pair

/* Differential suffix forest node.
 * 
 * This class handles all nodes of a differential suffix forest above the root.
 * We actually could call this class 'tree'...
 */
template< class Forest_options >
class Differential_suffix_forest_node {
 public:
  using Node = Differential_suffix_forest_node;  // rename itself
  using Idem = typename Forest_options::Idem;
  using Polynomial = typename Forest_options::Polynomial;
  using Gen_type = typename Forest_options::Gen_type;
  using Arc_iterator = typename Forest_options::Arc_container::iterator;  // handle?
    
  using Child_container = std::list< Node >;
  using Child_iterator = typename Child_container::iterator;
  using Child_const_iterator = typename Child_container::const_iterator;
  
  Differential_suffix_forest_node(Idem idem) :
    idem_(idem),
    type_(0),
    parent_(nullptr),
    criticalness_(0),
    polynomial_(1),
    weights_(0, 0),
    chain("")
  { }
  
  /* Copy constructor.
   * 
   * Copy a suffix forest, except for the differential arcs. Add the correct
   * parent pointers to children.
   */
  Differential_suffix_forest_node(const Differential_suffix_forest_node& other) :
    idem_(other.idem_),
    children_(other.children_),
    parent_(other.parent_),
    type_(other.type_),
    polynomial_(other.polynomial_),
    weights_(other.weights_),
    chain(other.chain),
    criticalness_(0)  // criticalness is not transferrable
  {
    for (auto& child : children_) {
      child.parent_ = this;
    }
  }
  
  /* Copy assignment.
   * 
   * Same idea
   */
  Node& operator=(const Differential_suffix_forest_node& other) {
    idem_ = other.idem_;
    children_ = other.children_;
    for (auto& child : children_) {
      child.parent_ = this;
    }
    parent_ = other.parent_;
    type_ = other.type_;
    polynomial_ = other.polynomial_;
    weights_ = other.weights_;
    chain = other.chain;
    criticalness_ = 0;  // criticalness is not transferrable
  }
  
  Differential_suffix_forest_node(Differential_suffix_forest_node&& other) = default;
  Node& operator=(Differential_suffix_forest_node&& other) = default;
  ~Differential_suffix_forest_node() = default;
  
  using Arc_iterator_const_range = Normal_const_range< std::list< Arc_iterator > >;
  
//   std::list< Arc_iterator > incoming_arcs() const {
//     return incoming_arcs_;
//   }
  
  Arc_iterator_const_range outgoing_arcs() {
    return Arc_iterator_const_range(outgoing_arcs_);
  }
  
  Child_iterator children_begin() {
    return children_.begin();
  }
  
  Child_iterator children_end() {
    return children_.end();
  }
  
  Child_const_iterator children_cbegin() const {
    return children_.cbegin();
  }
  
  Child_const_iterator children_cend() const {
    return children_.cend();
  }
  
  Child_container children() const {
    return children_;
  }
  
  Idem idem() const {
    return idem_;
  }
  
  bool has_parent() const {
    return (parent_ != nullptr);
  }
  
  bool has_children() const {
    return (!children_.empty());
  }
  
  Node* parent() const {
    return parent_;
  }
  
  Gen_type type() const {
    return type_;
  }
  
  /* TO DO: Check for cancelation?
   */
  void add_incoming_arc_it(Arc_iterator arc_it) {
    incoming_arcs_.push_back(arc_it);
  }
  
  void add_outgoing_arc_it(Arc_iterator arc_it) {
    outgoing_arcs_.push_back(arc_it);
  }
  
  /* Add child with an edge labeled by type. Return a pointer to the child
   * node.
   */
  Node* add_child(Gen_type type, const Node& child, std::string new_chain, int a=0, int m=0) {
    children_.emplace_back(child);
    children_.back().parent_ = this;
    children_.back().type_ = type;
    children_.back().add_weight(a, m);  // temporary?
    children_.back().add_to_chain(new_chain);
    return &(children_.back());
  }
  
  /* Requirement: the other descendant node is rel_ht layers higher
   * than this node */
  std::pair< Node*, bool > descendant_like(Node* other_descendant, int rel_ht) {
    if (rel_ht <= 0) {
      return std::make_pair(this, true);
    }
    else {
      auto result = descendant_like(other_descendant->parent(), rel_ht - 1);
      if (result.second) {
        Node* desc_parent = result.first;
        for (Child_iterator child_it = desc_parent->children_begin();
             child_it != desc_parent->children_end(); ++child_it) {
          if (child_it->type() == other_descendant->type()) {
            return std::make_pair(&*child_it, true);
          }
        }
      }
    }
    return std::make_pair(this, false);  // descendant not found
  }
  
  /* Note: the following two functions could be split into a public function
   * and an auxiliary private function.
   */
  /* Calculate arcs arriving at given differential arc, assuming that it is
   * higher than this node.
   */
  void calculate_arcs_neighboring_source(Arc_iterator arc_it, int rel_ht = 0) {
    for (auto incoming_arc_it : incoming_arcs_) {
      /* Check that the neighbor candidate has a compatible source path. */
      if ((arc_it->ht() <= incoming_arc_it->ht())
        or incoming_arc_it->source()->descendant_like(arc_it->source(), rel_ht).second) {
        arc_it->add_source_incoming_arc_it(incoming_arc_it);
        if (rel_ht > 0) {
          incoming_arc_it->add_target_outgoing_arc_it(arc_it);
        }
      }
    }
    
    for (auto outgoing_arc_it : outgoing_arcs_) {
      if (outgoing_arc_it == arc_it) { continue; }
      /* neighbor candidate is guaranteed to have a compatible source path */
      arc_it->add_source_outgoing_arc_it(outgoing_arc_it);
      if (rel_ht > 0) {
        outgoing_arc_it->add_source_outgoing_arc_it(arc_it);
      }
    }
    
    if (parent_ != nullptr) {
      parent_->calculate_arcs_neighboring_source(arc_it, rel_ht + 1);
    }
  }
  
  /* Calculate arcs leaving from given differential arc, assuming that it is
   * higher than this node.
   */
  void calculate_arcs_neighboring_target(Arc_iterator arc_it, int rel_ht = 0) {
    for (auto outgoing_arc_it : outgoing_arcs_) {
      /* Check that the neighbor candidate has a compatible target path. */
      if ((arc_it->ht() <= outgoing_arc_it->ht())
        or outgoing_arc_it->target()->descendant_like(arc_it->target(), rel_ht).second) {
        arc_it->add_target_outgoing_arc_it(outgoing_arc_it);
        if (rel_ht > 0) {
          outgoing_arc_it->add_source_incoming_arc_it(arc_it);
        }
      }
    }
    
    for (auto incoming_arc_it : incoming_arcs_) {
      if (incoming_arc_it == arc_it) { continue; }
      /* neighbor candidate is guaranteed to have a compatible target path */
      arc_it->add_target_incoming_arc_it(incoming_arc_it);
      if (rel_ht > 0) {
        incoming_arc_it->add_target_incoming_arc_it(arc_it);
      }
    }
    
    if (parent_ != nullptr) {
      parent_->calculate_arcs_neighboring_target(arc_it, rel_ht + 1);
    }
  }
  
  void clear_neighbors() {
    for (auto& child : children_) {
      child.clear_neighbors();
    }
    incoming_arcs_.clear();
    outgoing_arcs_.clear();
  }
  
  template< class DA_bimodule_it >
  void trim_leaves(int rel_ht, DA_bimodule_it da_bimodule_it) {
    if (rel_ht > 0) {
      --da_bimodule_it;
      for (auto& child : children_) {
        child.trim_leaves(rel_ht - 1, da_bimodule_it);
      }
    }
    else {
      poincare_polynomial(da_bimodule_it);  // calculate the Poincaré polynomial
      children_.clear();
    }
  }
  
  /* Calculate and return the Poincaré polynomial associated to the tree. That
   * is, we assume that there are no more arcs above this node.
   */
  template< class DA_bimodule_it >
  Polynomial poincare_polynomial(DA_bimodule_it da_bimodule_it) const {
    if (children_.empty()) {
      return 1;
    }
    else {
      Polynomial polynomial = 0;
      for (auto& child : children_) {
        auto child_da_it = da_bimodule_it;
        --child_da_it;
        auto child_poly = child.poincare_polynomial(child_da_it);
        child_poly *= (*da_bimodule_it)->alexander_maslov_weights(child.type());
        polynomial += child_poly;
      }
      return polynomial;
    }
  }
  
  /* When we know the polynomial is calculated */
  Polynomial poincare_polynomial() const {
    return polynomial_;
  }
  
  std::pair< int, int > weights() const {
    return weights_;
  }
  
  /* Discrete Morse theory */
  
  int criticalness() const {
    return criticalness_;
  }
  
  void set_tbd() {
    criticalness_ = 2;
  }
  
  void set_leaves_tbd() {
    if (children_.empty()) {
      criticalness_ = 2;
    }
    else {
      for (auto& child : children_) {
        child.set_leaves_tbd();
      }
    }
  }
  
  /* A node is:
   * - critical (0) if every child is critical,
   * - to-be-deleted (2) if every child is tbd,
   * - partially critical (1) otherwise.
   * */
  void set_criticalness() {
    if (children_.empty()) { return; }  // already marked
    
    criticalness_ = 0;
    for (auto& child : children_) {  // mark children
      child.set_criticalness();
    }
    
    for (auto& child : children_) {  // determine if critical
      if (child.criticalness() != 0) {
        criticalness_ = 2;
      }
    }
    for (auto& child : children_) {
      if (child.criticalness() != 2 and criticalness_ == 2) {  // determine if tbd
        criticalness_ = 1;
      }
    }
  }
  
  void erase_children() {  // only called if criticalness_ != 2
    for (Child_iterator child_it = children_.begin(); child_it != children_.end(); ) {
      if (child_it->criticalness() == 2) {
        child_it = children_.erase(child_it);
      }
      else {
        child_it->erase_children();
        ++child_it;
      }
    }
  }
  
  void add_to_chain(std::string new_chain) {
    chain += new_chain;
    for (auto& child : children_) {
      child.add_to_chain(new_chain);
    }
  }
  
  void add_weight(int a, int m) {
    weights_.first += a;
    weights_.second += m;
    for (auto& child : children_) {
      child.add_weight(a, m);
    }
  }
  
  int n_gen() {
    if (children_.empty()) {
      return 1;
    }
    else {
      int result = 0;
      for (auto& child : children_) {
        result += child.n_gen();
      }
      return result;
    }
  }
  
  void leaves() const {
    if (children_.empty()) {
      std::cout << "[dsfn] generator: " << chain << " " << weights_.first << " " << weights_.second << std::endl;
    }
    else {
      for (auto& child : children_) {
        child.leaves();
      }
    }
  }
  
  std::string chain;
  
 private:
  Idem idem_;  // at some point, this should be removed.
  
  Child_container children_;
    
  Node* parent_;
  Gen_type type_;
    
  std::list< Arc_iterator > incoming_arcs_;
  std::list< Arc_iterator > outgoing_arcs_;
  
  Polynomial polynomial_;
  std::pair< int, int > weights_;
 
  /* Discrete Morse theory */
  int criticalness_;  // 0 = critical, 1 = partially critical, 2 = to-be-deleted
};

#endif  // DIFFERENTIAL_SUFFIX_FOREST_NODE_