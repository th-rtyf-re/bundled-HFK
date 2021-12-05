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


/* Differential suffix forest.
 * 
 * Arcs are stored in a boost multi-index container, indexed by source and
 * target nodes. Thus, we can iterate over arcs in order of source or target.
 * This is useful for finding neighbors without too much overhead (maybe).
 */
template< class Forest_options = Forest_options_default_short >
class Forest {
 public:
  using Idem = typename Forest_options::Idem;
  using Gen_type = typename Forest_options::Gen_type;
  using Bordered_algebra = typename Forest_options::Bordered_algebra;
  using Alg_el = typename Forest_options::Alg_el;
  using Weights = typename Forest_options::Weights;
    
  struct Node {
    /* Constructor when we know the parent, and when there are no arcs yet.
     * Initializing next_sibling_ is a bit iffy.
     */
    Node(const int d = 0) :
      d_parent(d),
      descendants_offset(1),
      descendants_size(1),
      weights(0, 0),
      to_be_deleted(false)
    { }
    
    int d_parent;
    int descendants_offset;
    int descendants_size;
    Weights weights;
#ifdef DRAW
    std::string label;
#endif  // DRAW
    
    bool to_be_deleted;
    
    friend std::ostream& operator<<(std::ostream& os, const Node& node) {
      os << "<"
         << node.d_parent
         << "|"
#ifdef DRAW
         << node.label
         << "|"
#endif  // DRAW
         << node.descendants_size
         << "> ";
      return os;
    }
  };
  
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
  
  using Root_handle_container = std::map< int, Idem >;
  using Root_handle = typename Root_handle_container::value_type;
  
  /* Aliases for external usage */
  using Gen_bundle_handle_container = Root_handle_container;
  using Gen_bundle_handle = Root_handle;
  using Coef_bundle = Arc;
  using Coef_bundle_container = Arc_view;
  using Coef_bundle_iterator = Arc_iterator;
  using Coef_bundle_reference = Arc_reference;
  
  const Root_handle_container& gen_bundle_handles() const {
    return root_idems_;
  }
  
  const Arc_container& coef_bundles() const {
    return arcs_;
  }
  
  Idem idem(const Root_handle& root_handle) const {
    return root_handle.second;
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
    nodes_.clear();
    arcs_.clear();
    add_gen_bundle(Idem("0"));
    lock_generators();
    lock_coefficients();
  }
  
  /* Operations on nodes */
  
  /* Declare subtree to a specified idempotent. The extra argument at the end
   * is to force the user to acknowledge from which forest the old root comes.
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
    nodes_.clear();
    for (auto& map_value : declared_subtrees_) {
      int new_root = nodes_.size();
      const Idem& new_idem = map_value.first;
      root_idems_[new_root] = new_idem;
      nodes_.push_back(Node());  // add new root node
      
      for (auto& vec_value : map_value.second) {
        Gen_type& new_type = vec_value.first;
        int old_root = vec_value.second;
        int descendants_size = old_forest.descendants_size_(old_root);
        int new_child = nodes_.size();
        // Conditionally copy nodes from old forest
        std::copy_if(old_forest.nodes_.begin() + old_root,
                     old_forest.nodes_.begin() + old_root + descendants_size,
                     std::back_inserter(nodes_),
                     [](const Node& node) { return !node.to_be_deleted; });
        nodes_[new_child].d_parent = new_child - new_root;
        nodes_[new_child].weights = first_layer_weights[new_type];
#ifdef DRAW
        nodes_[new_child].label = first_layer_labels[new_type];
#endif  // DRAW
        nodes_[new_root].descendants_size += descendants_size;
        
        first_layer_nodes_[{new_idem, new_type}] = new_child;
      }
    }
    declared_subtrees_.clear();
  }
  
  /* Lock subtrees, but only the roots - to be removed */
  void lock_generators() {
    nodes_.clear();
    for (auto& map_value : declared_subtrees_) {
      int new_root = nodes_.size();
      const Idem& new_idem = map_value.first;
      root_idems_[new_root] = new_idem;
      nodes_.push_back(Node());
    }
  }
  
 private:
  /* Add old node, taking into account an offset for d_parent. Return the
   * number of newly added nodes and the offset for later siblings.
   */
  std::pair< int, int > lock_generators_rec_(int d_parent_offset, const int old_node, const Forest& old_forest) {
    int new_node = nodes_.size();
    nodes_.push_back(old_forest.nodes_[old_node]);
    nodes_[new_node].d_parent = nodes_[old_node].d_parent - d_parent_offset;
    nodes_[new_node].descendants_offset = 1;
    nodes_[new_node].descendants_size = 1;
    
    for (int child = old_forest.descendants_begin_(old_node);
         child != old_forest.descendants_end_(old_node);
         child += old_forest.descendants_size_(child)) {
      auto value_pair = lock_generators_rec_(d_parent_offset, child, old_forest);
      nodes_[new_node].descendants_size += value_pair.first;
      d_parent_offset += value_pair.second;
    }
    
    if (old_forest.no_children(old_node)) {
      return {1, old_forest.descendants_size_(old_node) - 1};
    }
    else {
      return {nodes_[new_node].descendants_size, d_parent_offset};
    }
  }
  
 public:
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
    source += old_forest.d_root_(old_arc.source);
    target += old_forest.d_root_(old_arc.target);
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
  std::vector< Node > node_container_() const {
    return nodes_;
  }
  
  int d_parent_(int node) const {
    return nodes_[node].d_parent;
  }
  
  int descendants_offset_(int node) const {
    return nodes_[node].descendants_offset;
  }
  
  int descendants_size_(int node) const {
    return nodes_[node].descendants_size;
  }
  
  bool to_be_deleted_(int node) const {
    return nodes_[node].to_be_deleted;
  }
  
  Weights weights_(int node) const {
    return nodes_[node].weights;
  }
  
  int parent_(int node) const {
    return node - d_parent_(node);
  }
  
  int descendants_begin_(int node) const {
    return node + descendants_offset_(node);
  }
  
  int descendants_end_(int node) const {
    return node + descendants_size_(node);
  }
  
  bool no_children_(int node) const {
    return descendants_begin_(node) == descendants_end_(node);
  }
  
  bool has_parent_(int node) const {
    return d_parent_(node) > 0;
  }
  
  void mark_to_be_deleted_(int node) {
    nodes_[node].to_be_deleted = true;
  }
  
  void increase_descendants_offset_(int node, int increase) {
    nodes_[node].descendants_offset += increase;
  }
  
  void increase_descendants_size_(int node, int increase) {
    nodes_[node].descendants_size += increase;
  }
  
  /* Find first non-deleted previous node
   */
  int previous_node_(int node) const {
    int previous_node = node - 1;
    while (previous_node >= 0 and to_be_deleted_(previous_node)) {
      --previous_node;
    }
    return previous_node;
  }
  
  int root_(int node) const {
    return (--root_idems_.upper_bound(node))->first;
  }
  
  /* Distance from a node to its root. Only used internally when declaring
   * arcs.
   */
  int d_root_(int node) const {
    return node - root_(node);
  }
  
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
    int n_nodes = descendants_size_(start_node);
    int end_node = start_node + n_nodes;
    std::vector< bool > marked(n_nodes, false);
    
    for (; upper_arc_it != arcs_.end() and upper_arc_it->source < end_node; ) {
      if (overlap_(*lower_arc_it, *upper_arc_it)) {
        // Mark ascendants
        int d_parent = d_parent_(upper_arc_it->source);
        int parent = upper_arc_it->source - d_parent;
        while (d_parent > 0 and parent >= start_node) {
          marked[parent - start_node] = true;
          d_parent = d_parent_(parent);
          parent = parent - d_parent;
        }
        
        // Skip arcs until after children
        int descendants_end = descendants_end_(upper_arc_it->source);
        for (; upper_arc_it != arcs_.end()
             and upper_arc_it->source < descendants_end; ++upper_arc_it);
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
      if (!marked[rel_node] and marked[rel_node - d_parent_(start_node + rel_node)]) {
        arcs_.emplace(start_node + rel_node, old_arc.target + rel_node, old_arc.value);
        rel_node += descendants_size_(start_node + rel_node);
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
    
    int d_parent = d_parent_(upper_node);
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
          d_parent = d_parent_(current_node);
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
      for (int child = descendants_begin_(parent);
           child != descendants_end_(parent);
           child += descendants_size_(child)) {
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
    int d_parent = d_parent_(node);
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
          d_parent = d_parent_(current_node);
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
    for (int child = descendants_begin_(node);
         child != descendants_end_(node);
         child += descendants_size_(child)) {
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
    for (int child = descendants_begin_(parent);
         child != descendants_end_(parent);
         child += descendants_size_(child)) {
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
    int descendants_end = descendants_end_(node);
    const auto& arcs_view = arcs_.template get< Tag >();
    auto arc_begin = arcs_view.lower_bound(node);
    auto arc_end = arcs_view.lower_bound(descendants_end);
    
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
    int descendants_end = descendants_end_(node);
    const auto& arcs_view = arcs_.template get< Tag >();
    auto arc_begin = arcs_view.lower_bound(node);
    auto arc_end = arcs_view.lower_bound(descendants_end);
    
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
    int d_parent = d_parent_(node);
    int current_node = node - d_parent;
    do {
      --arc_it;
      int endpoint = get_endpoint(*arc_it);
      if (endpoint == current_node) {
        arc_stream.emplace_back(*arc_it);
      }
      else if (endpoint < current_node) {
        do {
          d_parent = d_parent_(current_node);
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
             and front_arc.source < descendants_end_(back_arc.target))
         or (front_arc.source <= back_arc.target
             and back_arc.target < descendants_end_(front_arc.source)));
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
          std::cout << "[f] " << *arc_it << std::endl;
          reduction = true;
          arc_it = contract_(arc_it);
        }
        else {
          ++arc_it;
        }
      }
    }
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
        and front_diff < back_diff + descendants_size_(back_arc.target)) {
      source += front_diff - back_diff;
      arc_stream.emplace_back(source, target, back_arc.value * front_arc.value);
    }
    else if (front_diff <= back_diff
             and back_diff < front_diff + descendants_size_(front_arc.source)) {
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
    int d_parent = d_parent_(node);
    int parent = node - d_parent;
    // Check that node is the first and last child of parent
    while (d_parent == descendants_offset_(parent)
           and descendants_end_(node) == descendants_end_(parent)) {
      node = parent;
      d_parent = d_parent_(node);
      parent = node - 1;
    }
    return node;
  }
  
  /* Mark all nodes above and including the given node. Erase all arcs to and
   * from above the given node.
   */
  Arc_iterator delete_subtree_(int node) {
    update_after_deleting_(node);
    
    for (int i = node; i != descendants_end_(node); ++i) {
      mark_to_be_deleted_(i);
    }
    
    // If we delete a root, remove it from the generator bundle handles
    if (!has_parent_(node)) {
      root_idems_.erase(node);
    }
    
    delete_arcs_above_< Target >(node);
    return delete_arcs_above_< Source >(node);
  }
  
  void update_after_deleting_(int deleted_node) {
    int node = previous_node_(deleted_node);
    
    if (node < 0) { return; }
    
    int d_parent = 1;
    
    // increase descendant size for nodes on the right edge of the previous
    // subtree.
    while (d_parent != 0 and descendants_end_(node) == deleted_node) {
      increase_descendants_size_(node, descendants_size_(deleted_node));
      d_parent = d_parent_(node);
      node -= d_parent;
    }
    
    // Check if deleted node is the first child of its parent.
    d_parent = d_parent_(deleted_node);
    if (d_parent == descendants_offset_(parent_(deleted_node))) {
      increase_descendants_offset_(node, descendants_size_(deleted_node));
    }
  }
  
  /* Delete all arcs whose source or target is above a given node.
   */
  template< class Tag >
  typename Arc_container::template index< Tag >::type::iterator delete_arcs_above_(int node) {
    auto& arcs_view = arcs_.template get< Tag >();
    auto it_begin = arcs_view.lower_bound(node);
    auto it_end = arcs_view.lower_bound(descendants_end_(node));
    return arcs_view.erase(it_begin, it_end);
  }
  
 public:
  template< class Polynomial >
  Polynomial poincare_polynomial(const Idem& idem = Idem("0")) const {
    // find root with given idempotent
    for (const auto& value_pair : root_idems_) {
      if (value_pair.second == idem) {
        return poincare_polynomial_node_< Polynomial >(value_pair.first);
      }
    }
    return Polynomial(0);
  }
  
  
 private:
  template< class Polynomial >
  Polynomial poincare_polynomial_node_(const int node) const {
    if (no_children_(node)) {
      return Polynomial(1);
    }
    else {
      Polynomial poly = 0;
      for (int child = descendants_begin_(node);
           child != descendants_end_(node);
           child += descendants_size_(child)) {
        Polynomial child_poly = poincare_polynomial_node_< Polynomial >(child);
        child_poly *= weights_(node);
        poly += child_poly;
      }
      return poly;
    }
  }
  
 public:
  
  /* I/O interface */
  
  friend std::ostream& operator<<(std::ostream& os, const Forest& forest) {
    os << "Forest nodes: ";
    for (auto& node : forest.node_container_()) {
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
    
    for (auto value_pair : root_idems_) {
      add_to_grid_(grid, 0, value_pair.first);
    }
    
    write_file << "\\begin{tikzpicture}[suffix forest]" << std::endl;
    
    std::string extra_options("");
    
    /* Only put idempotents for the roots */
    int i = 0;
    for (auto value_pair : root_idems_) {
      Grid_point_& grid_point = grid[0][i];
      write_file << "\\node ("
                   << grid_point.first
                   << ") at ("
                   << grid_point.second
                   << ","
                   << 0
                   << ") {" << extra_options << "\\nodeLabel{"
                   << value_pair.second
                   << "}};" << std::endl;
      ++i;
    }
    
    for (int layer = 1; layer < grid.size(); ++layer) {
      for (Grid_point_& grid_point : grid[layer]) {
        write_file << "\\node ("
                   << grid_point.first
                   << ") at ("
                   << grid_point.second
                   << ","
                   << static_cast<float>(layer) * y_sep_
                   << ") {" << extra_options << "\\nodeLabel{0}};\n"
                   << "\\draw[->] ("
                   << grid_point.first
                   << ") -- node[in place]{$"
                   << nodes_[grid_point.first].label
                   << "$} ("
                   << (grid_point.first - nodes_[grid_point.first].d_parent)
                   << ");" << std::endl;
      }
    }
    
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
    if (no_children_(node)) {  // no children
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
      
      for (int child = node + 1;
           child != descendants_end_(node);
           child += descendants_size_(child)) {
        add_to_grid_(grid, layer + 1, child);
      }
      x = (grid[layer + 1][first_child_index].second + grid[layer + 1].back().second) / 2.;
      grid[layer].push_back(Grid_point_(node, x));
    }  // with children
  }
#endif  // DRAW
  
  /* Main data structures */
  std::vector< Node > nodes_;
  Arc_container arcs_;
  
  /* Auxiliary data structures */
  std::map< Idem, std::vector< std::pair< Gen_type, int > > > declared_subtrees_;
  std::vector< Arc > declared_arcs_;
  
  Root_handle_container root_idems_;
  std::map< std::pair< Idem, Gen_type >, int > first_layer_nodes_;
};

#endif  // TEST_FOREST_H_