#ifndef TEST_FOREST_H_
#define TEST_FOREST_H_

#include <fstream>
#include <functional>  // reference_wrapper
#include <iostream>
#include <map>
#include <string>
#include <utility>  // pair
#include <vector>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

// For default template values
#include "test_forest_options.h"
#include "test_node_container.h"
#include "test_arc_container.h"

/* Differential suffix forest.
 * 
 * Arcs are stored in a boost multi-index container, indexed by source and
 * target nodes. Thus, we can iterate over arcs in order of source or target.
 * This is useful for finding neighbors without too much overhead (maybe).
 */
template<
  class Forest_options = Forest_options_default_short,
  template< class > class Node_container_template = Node_container,
  template< class, template< class > class > class Arc_container_template = Arc_container
>
class Forest :
  private Arc_container_template< Forest_options, Node_container_template >
{
 public:
  using Idem = typename Forest_options::Idem;
  using Gen_type = typename Forest_options::Gen_type;
  using Bordered_algebra = typename Forest_options::Bordered_algebra;
  using Alg_el = typename Forest_options::Alg_el;
  using Weights = typename Forest_options::Weights;
  
  using Node_container = Node_container_template< Forest_options >;
  using typename Node_container::Root_handle;
  using typename Node_container::Root_handle_container;
  
  using Arc_container = Arc_container_template< Forest_options, Node_container_template >;
  using Arc = typename Arc_container::Arc;
  using Arc_view = typename Arc_container::Arc_view;
  using Arc_iterator = typename Arc_container::Arc_iterator;
  using Arc_reference = typename Arc_container::Arc_reference;
  using Source = typename Arc_container::Source;
  using Target = typename Arc_container::Target;
  
  /* Aliases for external usage */
  using Gen_bundle_handle_container = Root_handle_container;
  using Gen_bundle_handle = Root_handle;
  using Coef_bundle = Arc;
  using Coef_bundle_container = Arc_view;
  using Coef_bundle_iterator = Arc_iterator;
  using Coef_bundle_reference = Arc_reference;
  
  /* Inherited member functions from Node_container and Arc_container */
  
  using Node_container::idem;
  using Node_container::poincare_polynomial;
  
  using Arc_container::source_idem;
  using Arc_container::target_idem;
  using Arc_container::U_weights;
  
  using Arc_container::others_to_source;
  using Arc_container::others_from_target;
  using Arc_container::others_from_source;
  using Arc_container::others_to_target;
  
  using Arc_container::compatible;
  using Arc_container::concatenate;
  
 private:
  // Modifiers
  using Node_container::clear_nodes;
  using Node_container::push_back_root;
  using Node_container::push_back_subtree;
  
  // Relative distances
  using Node_container::to_root;
  using Node_container::to_parent;
  using Node_container::to_next;
  
  // Absolute node positions
  using Node_container::descendants_begin;
  using Node_container::descendants_end;
  using Node_container::descendants_size;
  
  // Bool functions
  using Node_container::is_root;
  using Node_container::has_children;
  
  // Functions for pruning
  using Node_container::node_offsets;
  using Node_container::prune_nodes;
  using Node_container::erase_subtree_nodes;
  
  // Arc iterators
  using Arc_container::arcs_begin;
  using Arc_container::arcs_end;
  
  // Modifiers
  using Arc_container::clear_arcs;
  using Arc_container::insert_arcs;
  using Arc_container::insert_arc;
  
  using Arc_container::resolve_overlaps_after;
  using Arc_container::resolve_overlaps_before;
  
  // CWG109 says that being unable to do this is not a problem :(
  // using Arc_container::erase_arcs_above_node;
  // using Arc_container::raise_arcs_below_node;
  
  using Arc_container::update_arc_endpoints;
  
  
 public:
  const Root_handle_container& gen_bundle_handles() const {
    return this->root_idems_;
  }
  
  const Arc_view& coef_bundles() const {
    return this->arcs_;
  }
  
  void set_as_trivial() {
    clear_nodes();
    clear_arcs();
    add_gen_bundle(Idem("0"));
    lock_generators();
    lock_coefficients();
  }
  
  /* Operations on nodes */
  
  /* Declare subtree to a specified idempotent.
   * Pre-condition: subtrees unlocked.
   */
  void add_gen_bundle(
    Idem new_idem,
    Gen_type new_type,
    Gen_bundle_handle root_handle
  ) {
    declared_subtrees_[new_idem].emplace_back(new_type, root_handle.first);
  }
  
  void add_gen_bundle(Idem new_idem) {
    declared_subtrees_[new_idem];
  }
  
  /* Lock subtrees using another forest.
   * Post-condition: subtrees locked.
   * 
   * For now, first_layer_weights is a vector, indexed by generator type. This
   * imposes that Gen_type can be converted to int.
   */
  void lock_generators(
    const Forest& old_forest,
    const std::vector< Weights >& first_layer_weights,
    const std::vector< std::string >& first_layer_labels
  ) {
    clear_nodes();
    for (auto& map_value : declared_subtrees_) {
      const Idem& new_idem = map_value.first;
      int new_root = push_back_root(new_idem);
      
      for (auto& vec_value : map_value.second) {
        Gen_type& new_type = vec_value.first;
        int old_root = vec_value.second;
        int new_child = push_back_subtree(
          new_root,
          first_layer_weights[new_type],
          first_layer_labels[new_type],
          old_root,
          old_forest);
        first_layer_nodes_[{new_idem, new_type}] = new_child;
      }
    }
    declared_subtrees_.clear();
  }
  
  /* Lock subtrees, but only the roots */
  void lock_generators() {
    clear_nodes();
    for (auto& map_value : declared_subtrees_) {
      const Idem& new_idem = map_value.first;
      push_back_root(new_idem);
    }
  }
  
  /* Arc creation
   * 
   * Arcs are the most annoying object to store.
   * 
   * Current approach: add arcs in a multi-indexed set that can sort
   * by source and by target.
   * 
   * When checking modulo 2 stuff, take arc, and scan forwards until the source
   * is greater or equal to the after-children index.
   */
  
  /* I don't know if I want this function. There are currently no methods that
   * take an algebra element as argument. Any such candidate methods can
   * also be written for coefficients.
   */
  template< class ...Args >
  Alg_el alg_el(Args&&... args) const {
    return Alg_el(args...);
  }
  
  void add_coef_bundle(
    const Alg_el& new_value,
    const Gen_type back_marking,
    const Gen_type front_marking,
    const Arc& old_arc,
    const Forest& old_forest
  ) {
    if (new_value.is_null()) { return; }
    int source = first_layer_nodes_.at({new_value.source_idem(), back_marking});
    int target = first_layer_nodes_.at({new_value.target_idem(), front_marking});
    source += old_forest.to_root(old_arc.source);
    target += old_forest.to_root(old_arc.target);
    declared_arcs_.emplace_back(source, target, new_value);
  }
  
  void add_coef_bundle(
    const Alg_el& new_value,
    const Gen_type back_marking,
    const Gen_type front_marking,
    const Idem& old_idem
  ) {
    if (new_value.is_null()) { return; }
    int source = first_layer_nodes_.at({new_value.source_idem(), back_marking});
    int target = first_layer_nodes_.at({new_value.target_idem(), front_marking});
    declared_arcs_.emplace_back(source, target, new_value);
  }
  
  void lock_coefficients() {
    insert_arcs(declared_arcs_.begin(), declared_arcs_.end());
    for (auto lower_arc_it = arcs_begin(); lower_arc_it != arcs_end(); ) {
      lower_arc_it = resolve_overlaps_after(lower_arc_it);
    }
  }
  
 public:
  /* Homotopy reduction of the forest to an irreducible one.
   */
  void reduce() {
    bool reduction = true;
    while (reduction) {
      reduction = false;
      for (auto arc_it = arcs_begin(); arc_it != arcs_end(); ) {
        if (arc_it->value.is_invertible()) {
          reduction = true;
          arc_it = contract_(arc_it);
        }
        else {
          ++arc_it;
        }
      }
    }
    const auto offsets = node_offsets();
    prune_nodes(offsets);
    update_arc_endpoints(offsets);
  }
  
 private:
  /* Contract an invertible arc. Return the next arc.
   * 
   * Note for mathematicians: it is not possible to have a back arc equal to
   * front arc, for grading reasons.
   */
  Arc_iterator contract_(const Arc_iterator& reverse_arc_it) {
    std::vector< Arc > zigzag_arcs;
    
    const Arc& reverse_arc = *reverse_arc_it;
    resolve_critical_(reverse_arc);
    
    for (const Arc& back_arc : others_to_target(reverse_arc)) {
      for (const Arc& front_arc : others_from_source(reverse_arc)) {
        std::cout << "[f] Making zig-zag arc from "
                  << back_arc << " "
                  << reverse_arc << " "
                  << front_arc << std::endl;
        zigzag_arcs = add_zigzag_(zigzag_arcs, back_arc, reverse_arc, front_arc);
      }
    }
    
    for (const Arc& zigzag_arc : zigzag_arcs) {
      std::cout << "[f] Adding zig-zag arc " << zigzag_arc << std::endl;
      insert_arc(zigzag_arc);
    }
    
    // get source and target before the iterator is invalidated
    int target = reverse_arc_it->target;
    int source = reverse_arc_it->source;
    
    erase_subtree_(greatest_single_child_ancestor_(target));
    auto next_arc_it = erase_subtree_(greatest_single_child_ancestor_(source));
    
    return next_arc_it;
  }
  
  /* Given a critical arc, raise every arc below it to get critical and non-
   * critical arcs.
   * 
   * A critical arc is one that will be deleted. The input critical arc is
   * generally an invertible arc from the reduction process.
   * 
   * This needs to be done for source and target.
   */
  void resolve_critical_(const Arc& critical_arc) {
    this->template raise_arcs_below_node< Source >(critical_arc.source);
    this->template raise_arcs_below_node< Source >(critical_arc.target);
    this->template raise_arcs_below_node< Target >(critical_arc.source);
    this->template raise_arcs_below_node< Target >(critical_arc.target);
  }
  
/* Check if the a zig-zag concatenation is possible, assuming that the back
   * and front arcs are compatible with the reverse arc, and then concatenate
   * if possible and add to a stream of new arcs.
   * 
   * This will be used in homotopy reduction.
   */
  std::vector< Arc >& add_zigzag_(
    std::vector< Arc >& arc_stream,
    const Arc& back_arc,
    const Arc& reverse_arc,
    const Arc& front_arc
  ) const {
    int source = back_arc.source;
    int target = front_arc.target;
    
    int back_diff = back_arc.target - reverse_arc.target;
    if (back_diff < 0) {
      source -= back_diff;
      back_diff = 0;
    }
    
    int front_diff = front_arc.source - reverse_arc.source;
    if (front_diff < 0) {
      target -= front_diff;
      front_diff = 0;
    }
    
    // At this point, the back and front arcs are at least as high as the
    // reverse arc.
    
    if (back_diff <= front_diff
        and front_diff < back_diff + descendants_size(back_arc.target)) {
      source += front_diff - back_diff;
      arc_stream.emplace_back(source, target, back_arc.value * front_arc.value);
    }
    else if (front_diff <= back_diff
             and back_diff < front_diff + descendants_size(front_arc.source)) {
      target += back_diff - front_diff;
      arc_stream.emplace_back(source, target, back_arc.value * front_arc.value);
    }
    return arc_stream;
  }
  
  /* Return the greatest ancestor of a given node such that the ancestor only
   * has one child.
   */
  int greatest_single_child_ancestor_(int node) const {
    int d_parent = to_parent(node);
    int parent = node - d_parent;
    // Check that node is the first and last child of parent
    while (descendants_size(parent) == descendants_size(node) + to_next(parent)) {
      node = parent;
      d_parent = to_parent(node);
      parent = node - 1;
    }
    return node;
  }
  
  /* Mark all nodes above and including the given node. Erase all arcs to and
   * from above the given node.
   */
  Arc_iterator erase_subtree_(int node) {
    erase_subtree_nodes(node);
    this->template erase_arcs_above_node< Target >(node);
    return this->template erase_arcs_above_node< Source >(node);
  }
  
 public:
  /* I/O interface */
  
  friend std::ostream& operator<<(std::ostream& os, const Forest& forest) {
    os << "Forest nodes: ";
    for (auto& node : forest.nodes_) {
      os << node;
    }
    os << std::endl;
    os << "Forest arcs: ";
    for (auto& arc : forest.coef_bundles()) {
      os << arc << " ";
    }
    os << std::endl;
    return os;
  }

#ifdef DRAW
  /* TeXify */
  void TeXify(std::ofstream& write_file) const {
    write_file << "\\begin{tikzpicture}[suffix forest]" << std::endl;
    
    Node_container::TeXify(write_file, false);
    Arc_container::TeXify(write_file);
    
    write_file << "\\end{tikzpicture}" << std::flush;
  }
#endif  // DRAW
  
 private:
  /* Auxiliary data structures */
  std::map< Idem, std::vector< std::pair< Gen_type, int > > > declared_subtrees_;
  std::vector< Arc > declared_arcs_;
  std::map< std::pair< Idem, Gen_type >, int > first_layer_nodes_;
};

#endif  // TEST_FOREST_H_