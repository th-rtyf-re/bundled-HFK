#ifndef TEST_FOREST_H_
#define TEST_FOREST_H_

#include <algorithm>  // copy_if
#include <fstream>
#include <functional>  // reference_wrapper
#include <iostream>
#include <iterator>  // back_inserter
#include <map>
#include <string>
#include <utility>  // pair
#include <vector>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include "test_forest_options.h"
#include "test_node_container.h"

/* Differential suffix forest.
 * 
 * Arcs are stored in a boost multi-index container, indexed by source and
 * target nodes. Thus, we can iterate over arcs in order of source or target.
 * This is useful for finding neighbors without too much overhead (maybe).
 */
template<
  class Forest_options = Forest_options_default_short,
  template< class > class Node_container_template = Node_container
>
class Forest : private Node_container_template< Forest_options > {
 public:
  using Idem = typename Forest_options::Idem;
  using Gen_type = typename Forest_options::Gen_type;
  using Bordered_algebra = typename Forest_options::Bordered_algebra;
  using Alg_el = typename Forest_options::Alg_el;
  using Weights = typename Forest_options::Weights;
  
  using Node_container = Node_container_template< Forest_options >;
  using typename Node_container::Root_handle;
  using typename Node_container::Root_handle_container;
  
  struct Arc {
    Arc(int s, int t, Alg_el v) :
      source(s),
      target(t),
      value(v)
    { }
    
    int source;
    int target;
    Alg_el value;
    
    bool operator==(const Arc& other) const {
      return (source == other.source
          and target == other.target
          and value == other.value);
    }
    
    bool operator!=(const Arc& other) const {
      return !(*this == other);
    }
    
    friend std::ostream& operator<<(std::ostream& os, const Arc& arc) {
      os << "("
         << arc.source
         << "|"
         << arc.value
         << "|"
         << arc.target
         << ")";
      return os;
    }
  };
  
  using Arc_container = boost::multi_index_container<
    Arc,
    boost::multi_index::indexed_by<
      boost::multi_index::ordered_non_unique<
        boost::multi_index::tag< struct Source >,
        boost::multi_index::member< Arc, int, &Arc::source >
      >,
      boost::multi_index::ordered_non_unique<
        boost::multi_index::tag< struct Target >,
        boost::multi_index::member< Arc, int, &Arc::target >
      >
    >
  >;
  
  using Arc_view = typename Arc_container::template index< Source >::type;
  using Arc_iterator = typename Arc_view::iterator;
  using Arc_reference = std::reference_wrapper< const Arc >;
  
//   using Root_handle_container = std::map< int, Idem >;
//   using Root_handle = typename Root_handle_container::value_type;
  
  /* Aliases for external usage */
  using Gen_bundle_handle_container = Root_handle_container;
  using Gen_bundle_handle = Root_handle;
  using Coef_bundle = Arc;
  using Coef_bundle_container = Arc_view;
  using Coef_bundle_iterator = Arc_iterator;
  using Coef_bundle_reference = Arc_reference;
  
  /* Inherited member functions from node container */
  
  using Node_container::idem;
  using Node_container::poincare_polynomial;
  
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
  using Node_container::get_offsets;
  using Node_container::prune_nodes;
  using Node_container::erase_subtree_nodes;
  
 public:
  const Root_handle_container& gen_bundle_handles() const {
    return this->root_idems_;
  }
  
  const Arc_container& coef_bundles() const {
    return arcs_;
  }
  
  Idem source_idem(const Arc& arc) const {
    return arc.value.source_idem();
  }
  
  Idem target_idem(const Arc& arc) const {
    return arc.value.target_idem();
  }
  
  std::vector< int > U_weights(const Arc& arc) const {
    return arc.value.U_weights();
  }
  
  void set_as_trivial() {
    clear_nodes();
    arcs_.clear();
    add_gen_bundle(Idem("0"));
    lock_generators();
    lock_coefficients();
  }
  
  /* Operations on nodes */
  
  /* Declare subtree to a specified idempotent.
   * Pre-condition: subtrees unlocked.
   */
  void add_gen_bundle(Idem new_idem, Gen_type new_type, Gen_bundle_handle root_handle) {
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
  void lock_generators(const Forest& old_forest,
                     const std::vector< Weights >& first_layer_weights,
                     const std::vector< std::string >& first_layer_labels) {
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
//   Alg_el alg_el(const Idem& source_idem,
//                    const Idem& target_idem,
//                    const std::vector< int >& U_weights) const {
//     return Alg_el(source_idem, target_idem, U_weights);
//   }
  
  void add_coef_bundle(const Idem& source_idem,
                   const Idem& target_idem,
                   const std::vector< int >& U_weights,
                   const Gen_type back_marking,
                   const Gen_type front_marking,
                   const Arc& old_arc,
                   const Forest& old_forest) {
    const Alg_el new_value(source_idem, target_idem, U_weights);
    if (new_value.is_null()) { return; }
    int source = first_layer_nodes_.at({source_idem, back_marking});
    int target = first_layer_nodes_.at({target_idem, front_marking});
    source += old_forest.to_root(old_arc.source);
    target += old_forest.to_root(old_arc.target);
    declared_arcs_.emplace_back(source, target, new_value);
  }
  
  void add_coef_bundle(const Idem& source_idem,
                   const Idem& target_idem,
                   const std::vector< int >& U_weights,
                   const Gen_type back_marking,
                   const Gen_type front_marking,
                   const Idem& old_idem) {
    const Alg_el new_value(source_idem, target_idem, U_weights);
    if (new_value.is_null()) { return; }
    int source = first_layer_nodes_.at({source_idem, back_marking});
    int target = first_layer_nodes_.at({target_idem, front_marking});
    declared_arcs_.emplace_back(source, target, new_value);
  }
  
  void lock_coefficients() {
    arcs_.insert(declared_arcs_.begin(), declared_arcs_.end());
    modulo_2_();
  }
  
 private:
  void modulo_2_() {
    for (auto lower_arc_it = arcs_.begin(); lower_arc_it != arcs_.end(); ) {
      lower_arc_it = resolve_overlaps_after_(lower_arc_it);
    }
  }
  
  /* A bunch of arc-raising methods.
   * 
   * One of the key operations on differential suffix forests is arc-raising.
   * This is not a single method, but rather a family of methods that operate
   * similarly. The idea is to transform a differential suffix forest into
   * another one representing the same D-module, but where some arcs are split
   * into several.
   * 
   * Arc-raising is used to cancel overlapping arcs and to safely homotopy
   * reduce the forest.
   */
  
  /* This is the definition of overlapping arcs. The lower and upper terms
   * are merely suggestions; the two arcs can be arranged arbitrarily.
   */
  bool overlap_(const Arc& lower_arc, const Arc& upper_arc) const {
    return (lower_arc.source + upper_arc.target == lower_arc.target + upper_arc.source
            and lower_arc.value == upper_arc.value);
  }
  
  /* Raise partially overlapping arcs (+ cancel exactly overlapping arcs)
   * 
   * Scan arcs above given arc. If one is found, mark its ancestors and move on
   * to the first arc after children.
   * 
   * The places to add a new arc are the unmarked nodes whose parent is marked.
   * 
   * Return the iterator of the next arc.
   */
  Arc_iterator resolve_overlaps_after_(Arc_iterator lower_arc_it) {
    auto upper_arc_it = std::next(lower_arc_it);
    
    // Scan for arcs that cancel completely
    for (; upper_arc_it != arcs_.end()
         and upper_arc_it->source == lower_arc_it->source; ++upper_arc_it) {
      if (*upper_arc_it == *lower_arc_it) {
        arcs_.erase(upper_arc_it);
        lower_arc_it = arcs_.erase(lower_arc_it);
        return lower_arc_it;
      }
    }
    
    // Scan descendants
    int start_node = lower_arc_it->source;
    int n_nodes = descendants_size(start_node);
    int end_node = start_node + n_nodes;
    std::vector< bool > marked(n_nodes, false);
    
    for (; upper_arc_it != arcs_.end() and upper_arc_it->source < end_node; ) {
      if (overlap_(*lower_arc_it, *upper_arc_it)) {
        // Mark ascendants
        int d_parent = to_parent(upper_arc_it->source);
        int parent = upper_arc_it->source - d_parent;
        while (d_parent > 0 and parent >= start_node) {
          marked[parent - start_node] = true;
          d_parent = to_parent(parent);
          parent = parent - d_parent;
        }
        
        // Skip arcs until after children
        int desc_end = descendants_end(upper_arc_it->source);
        for (; upper_arc_it != arcs_.end()
             and upper_arc_it->source < desc_end; ++upper_arc_it);
      }
      else {
        ++upper_arc_it;
      }
    }
    
    // Place new arcs
    raise_arcs_after_(*lower_arc_it, marked, start_node);
    
    if (marked[0]) {  // arcs were raised
      lower_arc_it = arcs_.erase(lower_arc_it);
      return lower_arc_it;
    }
    else {
      return std::next(lower_arc_it);
    }
  }
  
  /* Raise original arc to avoid all marked nodes. The vector marked is
   * indexed relatively to a start node.
   * 
   * !TBD! Marked nodes are never tbd.
   */
  void raise_arcs_after_(const Arc& old_arc,
                         const std::vector< bool >& marked,
                         const int start_node) {
    for (int rel_node = 1; rel_node != marked.size(); ) {  // relative node index
      if (!marked[rel_node] and marked[rel_node - to_parent(start_node + rel_node)]) {
        arcs_.emplace(start_node + rel_node, old_arc.target + rel_node, old_arc.value);
        rel_node += descendants_size(start_node + rel_node);
      }
      else {
        ++rel_node;
      }
    }
  }
  
  /* Unlike resolving above, where we may cancel more than two arcs at a time,
   * this finds at most one arc below that overlaps. This is used when adding
   * arcs to a container of arcs with no overlaps.
   * 
   * Returns true if a cancelation happened and false otherwise.
   */
  bool resolve_overlaps_before_(Arc_iterator upper_arc_it) {
    int upper_node = upper_arc_it->source;
    std::vector< int > ancestors = { upper_node };
    
    int d_parent = to_parent(upper_node);
    int current_node = upper_node;
    auto lower_arc_it = upper_arc_it;
    while (d_parent != 0 and lower_arc_it != arcs_.begin()) {
      --lower_arc_it;
      int endpoint = lower_arc_it->source;
      if (endpoint == current_node and overlap_(*lower_arc_it, *upper_arc_it)) {
        ancestors.push_back(current_node);
        raise_arcs_before_(*lower_arc_it, ancestors);
        arcs_.erase(lower_arc_it);
        arcs_.erase(upper_arc_it);
        return true;
      }
      else if (endpoint < current_node) {
        while (d_parent != 0 and endpoint < current_node) {
          ancestors.push_back(current_node);
          current_node -= d_parent;
          d_parent = to_parent(current_node);
        }
      }
    }
    
    return false;
  }
  
  /* ancestors holds source node and all of its ancestors up to the lower
   * overlapping arc.
   * 
   * !TBD! Check that children are not tbd.
   */
  void raise_arcs_before_(const Arc& old_arc,
                          const std::vector< int >& ancestors) {
    for (int i = 1; i != ancestors.size(); ++i) {
      int parent = ancestors[i];
      for (int child = descendants_begin(parent);
           child != descendants_end(parent);
           child += descendants_size(child)) {
        if (child != ancestors[i - 1]) {
          arcs_.emplace(child, old_arc.target + child - old_arc.source, old_arc.value);
        }
      }
    }
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
    raise_arcs_below_node_< Source >(critical_arc.source);
    raise_arcs_below_node_< Source >(critical_arc.target);
    raise_arcs_below_node_< Target >(critical_arc.source);
    raise_arcs_below_node_< Target >(critical_arc.target);
  }
  
  /* Auxiliary function for resolve_critical_. Note that this looks a lot like
   * resolve_overlaps_before_.
   */
  template< class Tag >
  void raise_arcs_below_node_(int node) {
    auto& arcs_view = arcs_.template get< Tag >();
    auto get_endpoint = arcs_view.key_extractor();  // how costly is this?
    
    int current_node = node;
    int d_parent = to_parent(node);
    int parent = current_node - d_parent;
    auto arc_it = arcs_view.lower_bound(node);
    
    std::vector< int > new_endpoints = children_(parent);
    
    // invariant at beginning of loop: endpoint of arc_it >= parent, and
    // new_endpoints contains all places to raise an arc at parent.
    while (d_parent != 0 and arc_it != arcs_view.begin()) {
      --arc_it;
      int endpoint = get_endpoint(*arc_it);
      
      if (endpoint == parent) {
        raise_arcs_below_node_aux_(arcs_view, *arc_it, endpoint, new_endpoints);
        arc_it = arcs_view.erase(arc_it);
      }
      else if (endpoint < parent) {
        while (d_parent != 0 and endpoint < parent) {
          current_node = parent;
          d_parent = to_parent(current_node);
          parent -= d_parent;
          add_other_children_(new_endpoints, parent, current_node);
        }
      }
    }
  }
  
  /* !TBD! Check that children are not tbd.
   */
  std::vector< int > children_(int node) const {
    std::vector< int > children;
    for (int child = descendants_begin(node);
         child != descendants_end(node);
         child += descendants_size(child)) {
      children.push_back(child);
    }
    return children;
  }
  
  /* Add the children of a parent to a vector of node indices, except for
   * a child that we avoid.
   * 
   * !TBD! Check that children are not tbd.
   */
  std::vector< int >& add_other_children_(std::vector< int >& node_stream,
                                          int parent,
                                          int avoiding) const {
    for (int child = descendants_begin(parent);
         child != descendants_end(parent);
         child += descendants_size(child)) {
      if (child != avoiding) {
        node_stream.push_back(child);
      }
    }
    return node_stream;
  }
  
  template< class Arc_view >
  void raise_arcs_below_node_aux_(Arc_view& arcs_view,
                                  const Arc& old_arc,
                                  const int old_endpoint,
                                  const std::vector< int >& new_endpoints) {
    for (const int endpoint : new_endpoints) {
      arcs_view.emplace(endpoint + old_arc.source - old_endpoint,
                        endpoint + old_arc.target - old_endpoint,
                        old_arc.value);
    }
  }
  
 public:
  /* Views on arcs */
  
  std::vector< Arc_reference > others_to_source(const Arc& arc) const {
    return arcs_at_node_< Target >(arc.source);
  }
  
  std::vector< Arc_reference > others_from_target(const Arc& arc) const {
    return arcs_at_node_< Source >(arc.target);
  }
  
  std::vector< Arc_reference > others_from_source(const Arc& arc) const {
    return arcs_at_node_< Source >(arc.source, arc);
  }
  
  std::vector< Arc_reference > others_to_target(const Arc& arc) const {
    return arcs_at_node_< Target >(arc.target, arc);
  }
  
 private:
   /* General function to get the arcs such that a specified endpoint is
   * compatible with a given node. The endpoint is specified by Tag.
   * The function returns a vector of references to arcs, probably
   * std::reference_wrapper< const Arc >.
   */
  template< class Tag >
  std::vector< Arc_reference > arcs_at_node_(int node) const {
    const auto& arcs_view = arcs_.template get< Tag >();
    auto arc_begin = arcs_view.lower_bound(node);
    auto arc_end = arcs_view.lower_bound(descendants_end(node));
    
    std::vector< Arc_reference > result(arc_begin, arc_end);
    
    if (arc_begin == arcs_view.begin()) {
      return result;
    }
    
    arcs_at_node_below_(result, arcs_view, arc_begin, node);
    
    return result;
  }
  
  /* Overloaded version where we avoid a specified arc.
   */
  template< class Tag >
  std::vector< Arc_reference > arcs_at_node_(int node, const Arc& avoiding) const {
    const auto& arcs_view = arcs_.template get< Tag >();
    auto arc_begin = arcs_view.lower_bound(node);
    auto arc_end = arcs_view.lower_bound(descendants_end(node));
    
    std::vector< Arc_reference > result;
    
    for (auto arc_it = arc_begin; arc_it != arc_end; ++arc_it) {
      // it sucks that we have to compare data structures instead of pointers...
      if (*arc_it != avoiding) {
        result.emplace_back(*arc_it);
      }
    }
    
    if (arc_begin == arcs_view.begin()) {
      return result;
    }
    
    arcs_at_node_below_(result, arcs_view, arc_begin, node);
    
    return result;
  }
  
  template< class Arc_view, class Arc_iterator = typename Arc_view::iterator >
  std::vector< Arc_reference >& arcs_at_node_below_(std::vector< Arc_reference >& arc_stream,
                                                    const Arc_view& arcs_view,
                                                    Arc_iterator arc_it,  // past-the-end iterator
                                                    int node) const {
    auto get_endpoint = arcs_view.key_extractor();  // how costly is this?
    int d_parent = to_parent(node);
    int current_node = node - d_parent;
    do {
      --arc_it;
      int endpoint = get_endpoint(*arc_it);
      if (endpoint == current_node) {
        arc_stream.emplace_back(*arc_it);
      }
      else if (endpoint < current_node) {
        do {
          d_parent = to_parent(current_node);
          current_node -= d_parent;
        } while (d_parent != 0 and endpoint < current_node);
      }
    } while (d_parent != 0 and arc_it != arcs_view.begin());
    
    return arc_stream;
  }
  
 public:
  
  /* Operations on arcs
   * 
   * These are the elementary operations we need: testing for compatibility,
   * concatenating, and concatenating zig-zags.
   */
  
  bool compatible(const Arc& back_arc, const Arc& front_arc) const {
    return ((back_arc.target <= front_arc.source
             and front_arc.source < descendants_end(back_arc.target))
         or (front_arc.source <= back_arc.target
             and back_arc.target < descendants_end(front_arc.source)));
  }
  
  /* Concatenate two arcs, assuming that they are concatenable */
  Arc concatenate(const Arc& back_arc, const Arc& front_arc) const {
    int difference = front_arc.source - back_arc.target;
    if (difference >= 0) {
      return Arc(back_arc.source + difference, front_arc.target, back_arc.value * front_arc.value);
    }
    else {
      return Arc(back_arc.source, front_arc.target - difference, back_arc.value * front_arc.value);
    }
  }
  
  /* Homotopy reduction of the forest to an irreducible one.
   */
  void reduce() {
    bool reduction = true;
    while (reduction) {
      reduction = false;
      for (auto arc_it = arcs_.begin(); arc_it != arcs_.end();) {
        if (arc_it->value.is_invertible()) {
          reduction = true;
          arc_it = contract_(arc_it);
        }
        else {
          ++arc_it;
        }
      }
    }
    const auto offsets = get_offsets();
    prune_nodes(offsets);
    update_arc_endpoints_(offsets);
  }
  
 private:
  
  
  /* Contract an invertible arc. Return the next
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
      insert_arc_(zigzag_arc);
    }
    
    int target = reverse_arc_it->target;
    int source = reverse_arc_it->source;
    
    delete_subtree_(greatest_single_child_ancestor_(target));
    
    auto next_arc_it = delete_subtree_(greatest_single_child_ancestor_(source));
    
    return next_arc_it;
  }
  
  /* Check if the a zig-zag concatenation is possible, assuming that the back
   * and front arcs are compatible with the reverse arc, and then concatenate
   * if possible and add to a stream of new arcs.
   * 
   * This will be used in homotopy reduction.
   */
  std::vector< Arc >& add_zigzag_(std::vector< Arc >& arc_stream,
                                  const Arc& back_arc,
                                  const Arc& reverse_arc,
                                  const Arc& front_arc) const {
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
  
  /* Insert arc into container of locked arcs. The danger is if this creates
   * overlapping arcs.
   */
  void insert_arc_(const Arc& arc) {
    auto arc_it = arcs_.insert(arc).first;
    bool cancel_below = resolve_overlaps_before_(arc_it);
    if (!cancel_below) {
      resolve_overlaps_after_(arc_it);
    }
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
  Arc_iterator delete_subtree_(int node) {
    erase_subtree_nodes(node);
    delete_arcs_above_< Target >(node);
    return delete_arcs_above_< Source >(node);
  }
    
  /* Delete all arcs whose source or target is above a given node.
   */
  template< class Tag >
  typename Arc_container::template index< Tag >::type::iterator delete_arcs_above_(int node) {
    auto& arcs_view = arcs_.template get< Tag >();
    auto it_begin = arcs_view.lower_bound(node);
    auto it_end = arcs_view.lower_bound(descendants_end(node));
    return arcs_view.erase(it_begin, it_end);
  }
  
  /* Update arc endpoints. No modification should rearrange arcs, so iterating
   * over them naively should work.
   * 
   * Moreover, 
   */
  void update_arc_endpoints_(const std::vector< int >& offsets) {
    for (auto arc_it = arcs_.begin(); arc_it != arcs_.end(); ++arc_it) {
      arcs_.modify(arc_it, [&](Arc& arc) { arc.source -= offsets[arc.source]; });
    }
    
    auto& arcs_target_view = arcs_.template get< Target >();
    for (auto arc_it = arcs_target_view.begin(); arc_it != arcs_target_view.end(); ++arc_it) {
      arcs_target_view.modify(arc_it, [&](Arc& arc) { arc.target -= offsets[arc.target]; });
    }
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
    Grid_ grid;
    
    write_file << "\\begin{tikzpicture}[suffix forest]" << std::endl;
    
    Node_container::TeXify(write_file, false);
    
    std::string extra_options("");
    for (auto& arc : arcs_) {
      write_file << "\\path ("
                 << arc.source
                 << ") edge[differential arc, bend left=10] node[in place]{$"
                 << extra_options << arc.value
                 << "$} ("
                 << arc.target
                 << ");" << std::endl;
    }
    write_file << "\\end{tikzpicture}" << std::flush;
  }
  
 private:
  using Grid_point_ = std::pair< int, float >;  // index of node & x coordinate
  using Grid_layer_ = std::vector< Grid_point_ >;  // left to right?
  using Grid_ = std::vector< Grid_layer_ >;  // lowest to highest
  
  static constexpr float min_x_sep_ = 1.;  // minimal horizontal distance between nodes
  static constexpr float y_sep_ = .8;  // vertical distance between nodes
  //static constexpr float poly_sep_ = .3;  // vertical distance between polynomial and leaf
  
  void add_to_grid_(Grid_& grid, int layer, int node) const {
    float x;
    if (!has_children(node)) {  // no children
      while (grid.size() <= layer) {  // make space for point
        grid.push_back(Grid_layer_());
      }
      if (grid[layer].empty()) {
        x = 0.;
      }
      else {
        x = grid[layer].back().second + min_x_sep_;
      }
      grid[layer].push_back(Grid_point_(node, x));
    }  // no children
    else {  // with children
      while (grid.size() <= layer + 1) {  // make space for children
        grid.push_back(Grid_layer_());
      }
      int first_child_index = grid[layer + 1].size();
      
      for (int child = descendants_begin(node);
           child != descendants_end(node);
           child += descendants_size(child)) {
        add_to_grid_(grid, layer + 1, child);
      }
      x = (grid[layer + 1][first_child_index].second + grid[layer + 1].back().second) / 2.;
      grid[layer].push_back(Grid_point_(node, x));
    }  // with children
  }
#endif  // DRAW
  
  /* Main data structures */
  Arc_container arcs_;
  
  /* Auxiliary data structures */
  std::map< Idem, std::vector< std::pair< Gen_type, int > > > declared_subtrees_;
  std::vector< Arc > declared_arcs_;
  
  std::map< std::pair< Idem, Gen_type >, int > first_layer_nodes_;
};

#endif  // TEST_FOREST_H_