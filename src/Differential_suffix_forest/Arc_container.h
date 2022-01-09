/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *  Bundled HFK - a knot Floer homology calculator                           *
 *                                                                           *
 *  Copyright (C) 2021-2022  Isaac Ren                                       *
 *  For further details, contact Isaac Ren (gopi3.1415@gmail.com).           *
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

#ifndef ARC_CONTAINER_H_
#define ARC_CONTAINER_H_

#include <iostream>
#include <vector>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include "Differential_suffix_forest_options.h"
#include "Node_container.h"

/* Arc container.
 * 
 * Arcs correspond to coefficient bundles (see explanation at
 * Differential_suffix_forest).
 * 
 * Arcs are stored in a boost multi-index container, indexed by source and
 * target nodes. Thus, we can iterate over arcs in order of source or target.
 * This is useful for finding neighbors without too much overhead (maybe).
 * 
 * Mathematically speaking, this is a D-module without homotopy reduction
 */
template< class Forest_options = Forest_options_default_short >
class Arc_container : protected Node_container< Forest_options > {
 public:
  using Idem = typename Forest_options::Idem;
  using Alg_el = typename Forest_options::Alg_el;
  
  using Node_container = Node_container< Forest_options >;
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
  
  struct Source { };
  struct Target { };
  
  using Arc_multi_index_container = boost::multi_index_container<
    Arc,
    boost::multi_index::indexed_by<
      boost::multi_index::ordered_non_unique<
        boost::multi_index::tag< Source >,
        boost::multi_index::member< Arc, int, &Arc::source >
      >,
      boost::multi_index::ordered_non_unique<
        boost::multi_index::tag< Target >,
        boost::multi_index::member< Arc, int, &Arc::target >
      >
    >
  >;
  
  using Arc_view = typename Arc_multi_index_container::template index< Source >::type;
  using Arc_iterator = typename Arc_view::iterator;
  using Arc_reference = std::reference_wrapper< const Arc >;
  
 private:
  // Relative distances
  using Node_container::to_parent;
  
  // Absolute node positions
  using Node_container::descendants_begin;
  using Node_container::descendants_end;
  using Node_container::descendants_size;
  
 public:
  
  /* OBSERVERS */
    
  Idem source_idem(const Arc& arc) const {
    return arc.value.source_idem();
  }
  
  Idem target_idem(const Arc& arc) const {
    return arc.value.target_idem();
  }
  
  std::vector< int > U_weights(const Arc& arc) const {
    return arc.value.U_weights();
  }
  
  int U_weight(const Arc& arc, int position) const {
    return arc.value.U_weight(position);
  }
  
  /* Views on arcs
   * 
   * There are two contexts in which we need access to arcs: while constructing
   * the set of arcs, and after constructing the set of arcs. In the first case,
   * we use the get_[arc type] functions, which calculate arcs at the moment
   * of the call. In the second case, we use the [arc type] functions, which
   * access pre-calculated vectors of arc references.
   * 
   * Our usage of these functions is as follows:
   * - The get_[arc type] functions are used when inserting zig-zag arcs during
   * reduction.
   * - The [arc type] functions are used by Morse events.
   */
  
  Arc_iterator arcs_begin() const {
    return arcs_.begin();
  }
  
  Arc_iterator arcs_end() const {
    return arcs_.end();
  }
  
  std::vector< Arc_reference > get_others_to_source(const Arc& arc) {
    return get_arcs_at_node_< Target >(arc.source);
  }
  
  std::vector< Arc_reference > get_others_from_target(const Arc& arc) {
    return get_arcs_at_node_< Source >(arc.target);
  }
  
  std::vector< Arc_reference > get_others_from_source(const Arc& arc) {
    return get_arcs_at_node_< Source >(arc.source, arc);
  }
  
  std::vector< Arc_reference > get_others_to_target(const Arc& arc) {
    return get_arcs_at_node_< Target >(arc.target, arc);
  }
  
  const std::vector< Arc_reference > others_to_source(const Arc& arc) const {
    return arcs_to_node_.at(arc.source);
  }
  
  const std::vector< Arc_reference > others_from_target(const Arc& arc) const {
    return arcs_from_node_.at(arc.target);
  }
  
  const std::vector< Arc_reference > others_from_source(const Arc& arc) const {
    return arcs_from_node_.at(arc.source); // need to avoid arc
  }
  
  const std::vector< Arc_reference > others_to_target(const Arc& arc) const {
    return arcs_to_node_.at(arc.target); // need to avoid arc
  }
  
 private:
  /* General function to get the arcs such that a specified endpoint is
   * compatible with a given node. The endpoint is specified by Tag.
   * The function returns a vector of references to arcs.
   */
  template< class Tag >
  std::vector< Arc_reference > get_arcs_at_node_(int node) const {
    const auto& arcs_view = arcs_.template get< Tag >();
    auto arc_begin = arcs_view.lower_bound(node);
    auto arc_end = arcs_view.lower_bound(descendants_end(node));
    
    std::vector< Arc_reference > result(arc_begin, arc_end);
    
    if (arc_begin == arcs_view.begin()) {
      return result;
    }
    
    arcs_below_node_(result, arcs_view, arc_begin, node);
    
    return result;
  }
  
  /* Overloaded version where we avoid a specified arc.
   */
  template< class Tag >
  std::vector< Arc_reference > get_arcs_at_node_(int node, const Arc& avoiding) const {
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
    
    arcs_below_node_(result, arcs_view, arc_begin, node);
    
    return result;
  }
  
  /* Arcs strictly below a node */
  template< class Arc_view, class Arc_iterator >
  std::vector< Arc_reference >& arcs_below_node_(
    std::vector< Arc_reference >& arc_stream,
    const Arc_view& arcs_view,
    Arc_iterator arc_it,  // past-the-end iterator
    int node
  ) const {
    auto get_endpoint = arcs_view.key_extractor();  // how costly is this?
    
    auto parent_it = ++this->ascender(node);
    int endpoint;
    while (parent_it.valid() and arc_it != arcs_view.begin()) {
      --arc_it;
      endpoint = get_endpoint(*arc_it);
      
      if (endpoint < *parent_it) {
        while (parent_it.valid() and endpoint < *parent_it) {
          ++parent_it;
        }
      }
      
      if (endpoint == *parent_it) {
        arc_stream.emplace_back(*arc_it);
      }
    }
    
    return arc_stream;
  }
  
 public:
  /* MODIFIERS */
  
  void clear_arcs() {
    arcs_.clear();
  }
  
  template< class Iterator >
  void insert_arcs(Iterator first, Iterator last) {
    arcs_.insert(first, last);
  }
  
  /* Insert arc into container of locked arcs. The danger is if this creates
   * overlapping arcs.
   */
  void insert_arc(const Arc& arc) {
    auto arc_it = arcs_.insert(arc).first;
    bool cancel_below = resolve_overlaps_before_(arc_it);
    if (!cancel_below) {
      resolve_overlaps_after_(arc_it);
    }
  }
  
  void basic_insert_arc(const Arc& arc) {
    arcs_.insert(arc);
  }
  
  /* Delete all arcs whose source or target is above a given node.
   */
  template<
    class Tag,
    class Iterator =
      typename Arc_multi_index_container::template index< Tag >::type::iterator
  >
  Iterator erase_arcs_above_node(int node) {
    auto& arcs_view = arcs_.template get< Tag >();
    auto it_begin = arcs_view.lower_bound(node);
    auto it_end = arcs_view.lower_bound(descendants_end(node));
    return arcs_view.erase(it_begin, it_end);
  }
  
  /* Calculate which arcs are adjacent to which nodes. This is used in the box
   * tensor product, where we no longer modify the forest.
   */
  void compute_arcs_at_nodes() {
    arcs_from_node_.clear();
    arcs_to_node_.clear();
    for (auto& arc : arcs_) {
      arcs_from_node_[arc.source];
      arcs_from_node_[arc.target];
      arcs_to_node_[arc.source];
      arcs_to_node_[arc.target];
    }
    compute_arcs_at_node_< Source >(arcs_from_node_);
    compute_arcs_at_node_< Target >(arcs_to_node_);
  }
  
 private:
  template< class Tag >
  void compute_arcs_at_node_(std::map< int, std::vector< Arc_reference > >& arcs_at_node) const {
    const auto& arcs_view = arcs_.template get< Tag >();
    auto get_endpoint = arcs_view.key_extractor();  // how costly is this?
    
    auto arc_it = arcs_view.begin();
    
    for (auto lower_node_it = arcs_at_node.begin(); lower_node_it != arcs_at_node.end(); ++lower_node_it) {
      const int lower_node = lower_node_it->first;
      
      for (; arc_it != arcs_view.end() and get_endpoint(*arc_it) < lower_node; ++arc_it) { }
      auto lower_arc_begin = arc_it;
      
      for (; arc_it != arcs_view.end() and get_endpoint(*arc_it) == lower_node; ++arc_it) {
        arcs_at_node[lower_node].emplace_back(*arc_it);
      }
      auto lower_arc_end = arc_it;
      
      auto upper_arc_it = arc_it;
      for (
        auto upper_node_it = std::next(lower_node_it);
        upper_node_it != arcs_at_node.end()
          and upper_node_it->first < this->descendants_end(lower_node);
        ++upper_node_it
      ) {
        int upper_node = upper_node_it->first;
        for (; upper_arc_it != arcs_view.end() and get_endpoint(*upper_arc_it) == upper_node; ++upper_arc_it) {
          lower_node_it->second.emplace_back(*upper_arc_it);
        }
        for (auto lower_arc_it = lower_arc_begin; lower_arc_it != lower_arc_end; ++lower_arc_it) {
          arcs_at_node[upper_node].emplace_back(*lower_arc_it);
        }
      }
    }
  }
  
 public:
  
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
  
  void modulo_2() {
    for (
      auto lower_arc_it = arcs_.begin();
      lower_arc_it != arcs_.end();
    ) {
      lower_arc_it = resolve_overlaps_after_(lower_arc_it);
    }
  }
  
  /* Raise all arcs strictly below the given node. As a result, there are no
   * arcs below the node.
   */
  template< class Tag >
  void raise_arcs_below_node(int node) {
    auto& arcs_view = arcs_.template get< Tag >();
    auto get_endpoint = arcs_view.key_extractor();  // how costly is this?
    
    auto parent_it = ++this->ascender(node);
    auto arc_it = arcs_view.lower_bound(node);
    
    std::vector< int > new_endpoints = children_(*parent_it);
    
    // new_endpoints contains all places to raise an arc at *parent_it.
    while (parent_it.valid() and arc_it != arcs_view.begin()) {
      --arc_it;
      int endpoint = get_endpoint(*arc_it);
      
      if (endpoint < *parent_it) {
        do {
          node = *parent_it;
          ++parent_it;
          add_other_children_(new_endpoints, *parent_it, node);
        } while (parent_it.valid() and endpoint < *parent_it);
      }
      
      if (endpoint == *parent_it) {
        // raise the arc pointed to by arc_it
        for (const int new_endpoint : new_endpoints) {
          arcs_view.emplace(
            new_endpoint + arc_it->source - endpoint,
            new_endpoint + arc_it->target - endpoint,
            arc_it->value
          );
        }
        arc_it = arcs_view.erase(arc_it);
      }
    }
  }
  
 private:
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
    std::vector< bool > marked(n_nodes, false);  // parent of new endpoints
    std::vector< bool > except(n_nodes, false);  // except these enpoints
    
    while (upper_arc_it != arcs_.end() and upper_arc_it->source < end_node) {
      if (overlap_(*lower_arc_it, *upper_arc_it)) {
        // Mark ascendants
        except[upper_arc_it->source - start_node] = true;
        auto parent_it = ++this->ascender(upper_arc_it->source);
        while (
          parent_it.valid()
          and *parent_it >= start_node
          and !marked[*parent_it - start_node]
        ) {
          marked[*parent_it - start_node] = true;
          ++parent_it;
        }
        
        // Skip arcs until after children
        int desc_end = descendants_end(upper_arc_it->source);
        upper_arc_it = arcs_.erase(upper_arc_it);
        for (; upper_arc_it != arcs_.end()
             and upper_arc_it->source < desc_end; ++upper_arc_it);
      }
      else {
        ++upper_arc_it;
      }
    }
    
    // Place new arcs
    raise_arcs_after_(*lower_arc_it, marked, except, start_node);
    
    if (marked[0]) {  // arcs were raised
      lower_arc_it = arcs_.erase(lower_arc_it);
      return lower_arc_it;
    }
    else {
      return std::next(lower_arc_it);
    }
  }
  
  /* Unlike resolving above, where we may cancel more than two arcs at a time,
   * this finds at most one arc below that overlaps. This is used when adding
   * arcs to a container of arcs with no overlaps.
   * 
   * Returns true if a cancelation happened and false otherwise.
   */
  bool resolve_overlaps_before_(Arc_iterator upper_arc_it) {
    auto current_it = this->ascender(upper_arc_it->source);
    std::vector< int > ancestors = { *current_it };
    
    auto lower_arc_it = upper_arc_it;
    while (current_it.valid() and lower_arc_it != arcs_.begin()) {
      --lower_arc_it;
      int endpoint = lower_arc_it->source;
      
      while (current_it.valid() and endpoint < *current_it) {
        ++current_it;
        ancestors.push_back(*current_it);
      }
      
      if (endpoint == *current_it and overlap_(*lower_arc_it, *upper_arc_it)) {
        raise_arcs_before_(*lower_arc_it, ancestors);
        arcs_.erase(lower_arc_it);
        arcs_.erase(upper_arc_it);
        return true;
      }
    }
    
    return false;
  }
  
  /* This is the definition of overlapping arcs. The lower and upper terms
   * are merely suggestions; the two arcs can be arranged arbitrarily.
   */
  bool overlap_(const Arc& lower_arc, const Arc& upper_arc) const {
    return (lower_arc.source + upper_arc.target == lower_arc.target + upper_arc.source
            and lower_arc.value == upper_arc.value);
  }
  
  /* Raise original arc to avoid all marked nodes. The vector marked is
   * indexed relatively to a start node.
   * 
   * Raising arcs means placing an arc at each unmarked node whose parent is
   * marked.
   */
  void raise_arcs_after_(
    const Arc& old_arc,
    const std::vector< bool >& marked,
    const std::vector< bool >& except,
    const int start_node
  ) {
    for (int rel_node = 1; rel_node != marked.size(); ) {  // relative node index
      if (
        !marked[rel_node]
        and !except[rel_node]
        and marked[rel_node - to_parent(start_node + rel_node)]
      ) {
        arcs_.emplace(start_node + rel_node, old_arc.target + rel_node, old_arc.value);
        rel_node += descendants_size(start_node + rel_node);
      }
      else {
        rel_node += this->to_next(start_node + rel_node);
      }
    }
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
  
  std::vector< int > children_(int node) const {
    std::vector< int > children;
    for (
      int child = descendants_begin(node);
      child != descendants_end(node);
      child += descendants_size(child)
    ) {
      children.push_back(child);
    }
    return children;
  }
  
  /* Add the children of a parent to a vector of node indices, except for
   * a child that we avoid.
   */
  std::vector< int >& add_other_children_(
    std::vector< int >& node_stream,
    int parent,
    int avoiding
  ) const {
    for (
      int child = descendants_begin(parent);
      child != descendants_end(parent);
      child += descendants_size(child)
    ) {
      if (child != avoiding) {
        node_stream.push_back(child);
      }
    }
    return node_stream;
  }
  
 public:
  /* Miscellaneous other functions */
  
  bool compatible(const Arc& back_arc, const Arc& front_arc) const {
    return (
      (
        back_arc.target <= front_arc.source
        and front_arc.source < descendants_end(back_arc.target)
      )
      or (
        front_arc.source <= back_arc.target
        and back_arc.target < descendants_end(front_arc.source)
      )
    );
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
  
  /* Update arc endpoints. No modification should rearrange arcs, so iterating
   * over them naively should work, with linear cost.
   * 
   * Moreover, we cannot factor these two loops using a key extractor, since
   * we are modifying the key...
   */
  void update_arc_endpoints(const std::vector< int >& offsets) {
    for (auto arc_it = arcs_.begin(); arc_it != arcs_.end(); ++arc_it) {
      arcs_.modify(arc_it, [&](Arc& arc) { arc.source -= offsets[arc.source]; });
    }
    
    auto& arcs_target_view = arcs_.template get< Target >();
    for (auto arc_it = arcs_target_view.begin(); arc_it != arcs_target_view.end(); ++arc_it) {
      arcs_target_view.modify(arc_it, [&](Arc& arc) { arc.target -= offsets[arc.target]; });
    }
  }
  
  /* I/O interface */
  
  friend std::ostream& operator<<(std::ostream& os, const Arc_container& ac) {
    for (auto& arc : ac.arcs_) {
      os << arc << " ";
    }
    return os;
  }
  
#ifdef BUNDLED_HFK_DRAW_
  /* TeXify */
  void TeXify(std::ofstream& write_file) const {
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
  }
#endif  // BUNDLED_HFK_DRAW_
  
 protected:
  Arc_multi_index_container arcs_;
  
  std::map< int, std::vector< Arc_reference > > arcs_from_node_;
  std::map< int, std::vector< Arc_reference > > arcs_to_node_;
};

#endif  // ARC_CONTAINER_H_