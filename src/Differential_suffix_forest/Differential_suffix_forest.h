/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 *                                                                           *
 *  Bundled HFK - a knot Floer homology calculator                           *
 *                                                                           *
 *  Copyright (C) 2021  Isaac Ren                                            *
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

#ifndef DIFFERENTIAL_SUFFIX_FOREST_H_
#define DIFFERENTIAL_SUFFIX_FOREST_H_

#include <algorithm>  // max_element
#include <fstream>
#include <initializer_list>
#include <iterator>  // next
#include <list>
#include <map>
#include <vector>

#include "Utility/Iterators.h"

/* Differential suffix forest.
 * 
 * This is our implementation of bundled D-modules. This class satisfies the
 * concept BundledDModule.
 * 
 * This is also where bundled discrete Morse theory is implemented.
 * 
 * TO DO: Completely rework this class, as well as arcs and nodes, so that arcs
 * no longer need to store their neighbors explicitly. This involves:
 * - putting all nodes in one vector, sorted by depth-first traversal
 * - having two lists of pointers to arcs, sorting them by depth-first traversal
 *   based on source and target
 * - defining a custom iterator that will find all neighboring arcs on the fly
 * - etc.
 */
template< class Forest_options >
class Differential_suffix_forest {
 public:
  using Node = typename Forest_options::Node;
  using Arc = typename Forest_options::Arc;
  using Idem = typename Forest_options::Idem;
  using Bordered_algebra = typename Forest_options::Bordered_algebra;
  using Alg_el = typename Forest_options::Alg_el;
  using Polynomial = typename Forest_options::Polynomial;
  
  using Gen_type = typename Forest_options::Gen_type;
  using Node_container = typename Forest_options::Node_container;
  using Node_iterator = typename Node_container::iterator;
  using Node_const_iterator = typename Node_container::const_iterator;
  using Arc_container = typename Forest_options::Arc_container;
  using Arc_iterator = typename Arc_container::iterator;
  using Arc_const_iterator = typename Arc_container::const_iterator;
  
  /* Handles.
   * 
   * Nodes have a handle, a const_iterator.
   * 
   * Arcs do not have a handle. Users can access them directly, and call any
   * member function of the AgglomeratedCoefficient concept.
   */
  using Gen_bundle_handle = Node_const_iterator;
  using Gen_bundle_handle_iterator = Self_indirection_const_iterator< Node_container >;
  using Gen_bundle_handle_range = Custom_const_range< Gen_bundle_handle_iterator >;
  
  //using Neighbor_range = typename Arc::Neighbor_range;
  //using Neighbor_const_range = typename Arc::Neighbor_const_range;
  using Coef_bundle = Arc;
  
  const Gen_bundle_handle_range gen_bundle_handles() const {
    return Gen_bundle_handle_range(root_nodes_);
  }
  
  const Arc_container coef_bundles() const {
    return differential_arcs_;
  }
  
  /* SECTION: Tentative for a new handle, for lower generators added during
   * tensor product. Only if further optimization is needed
   */
//   std::vector< Node* > lower_gen_handles(const Idem& idem) const {
//     return layer_one_nodes_[idem];
//   }
//   
//   Gen_type lower_gen_marking(const Node* lower_gen_handle) const {
//     return lower_gen_handle->type();
//   }
  /* end section */
  
  /* SECTION: If we shared arcs using a handle, these would be the operations:
   * - algebra_element(const Aggl_coef_handle&);
   * - concatenate(const Aggl_coef_handle&, const Aggl_coef_handle&);
   * - light_concatenate(const Aggl_coef_handle&, const Aggl_coef_handle&);
   */
  
  /* Constructors, destructors
   */
  Differential_suffix_forest() = default;
  Differential_suffix_forest(Differential_suffix_forest&) = delete;
  Differential_suffix_forest& operator=(Differential_suffix_forest&) = delete;
  Differential_suffix_forest(Differential_suffix_forest&&) = default;
  Differential_suffix_forest& operator=(Differential_suffix_forest&&) = default;
  ~Differential_suffix_forest() = default;
  
  /* Access and creation using handles */
  const Idem idem(const Gen_bundle_handle& gen_bundle_handle) const {
    return gen_bundle_handle->first;
  }
  
  /* Emplace node and return the emplaced node. Currently returns a pointer to
   * a node on layer 1
   * 
   * TO DO: replace gen_handle with a handle to an agglomerated generator
   */
  void add_gen_bundle(const Idem& new_idem,
                      const Gen_type marking,
                      const Gen_bundle_handle& gen_handle,
                      std::string new_chain,
                      std::pair< int, int > weights = {0, 0}  // temporary?
                      ) {
    auto result = root_nodes_.emplace(new_idem, new_idem);
    // result.first->second is the emplaced node:
    Node* new_child = result.first->second.add_child(marking, node_(gen_handle), new_chain, weights.first, weights.second);
    layer_one_nodes_[node_(gen_handle).idem()].push_back(new_child);
  }
  
  /* Version when there is no old agglomerated generator to copy over */
  void add_gen_bundle(const Idem& new_idem) {
    root_nodes_.emplace(new_idem, new_idem);
  }
  
  /* Set the D-module as trivial, i.e. representing a knot diagram with no
   * strands.
   */
  void set_as_trivial() {
    root_nodes_.clear();
    differential_arcs_.clear();
    add_gen_bundle(Idem("0"));
  }
  
  /* The construction of the upper arc is left to the DA-bimodules. The D-module
   * takes this upper arc, and composes it with lower markings and the new
   * algebra element to make the new arc.
   */
  void add_coef_bundle(const Alg_el& alg_el,
                       const Gen_type back_marking,
                       const Gen_type front_marking,
                       const typename Arc::Light& upper_arc) {
     Node* new_source = find_layer_one_node_(upper_arc.algebra_element.source_idem(), back_marking)
      ->descendant_like(upper_arc.source, upper_arc.ht).first;
    Node* new_target = find_layer_one_node_(upper_arc.algebra_element.target_idem(), front_marking)
      ->descendant_like(upper_arc.target, upper_arc.ht).first;
    differential_arcs_.emplace_back(new_source, new_target, alg_el, upper_arc.ht + 1);
  }
  
  void add_coef_bundle(const Alg_el& alg_el,
                       const Gen_type back_marking,
                       const Gen_type front_marking,
                       const Arc& upper_arc) {
     Node* new_source = find_layer_one_node_(upper_arc.algebra_element().source_idem(), back_marking)
      ->descendant_like(upper_arc.source(), upper_arc.ht()).first;
    Node* new_target = find_layer_one_node_(upper_arc.algebra_element().target_idem(), front_marking)
      ->descendant_like(upper_arc.target(), upper_arc.ht()).first;
    differential_arcs_.emplace_back(new_source, new_target, alg_el, upper_arc.ht() + 1);
  }
  
  /* If there are no upper arcs, use this version */
  void add_coef_bundle(const Alg_el& alg_el,
                       const Gen_type back_marking,
                       const Gen_type front_marking,
                       const Idem& upper_idem) {
    Node* new_source = find_layer_one_node_(upper_idem, back_marking);
    Node* new_target = find_layer_one_node_(upper_idem, front_marking);
    differential_arcs_.emplace_back(new_source, new_target, alg_el, 1);
  }
  
 private:
  const Node& node_(const Gen_bundle_handle& handle) const {
    return handle->second;
  }
  
  /* Assume that the node exists */
  Node* find_layer_one_node_(const Idem& idem, const Gen_type marking) {
    auto& node_list = layer_one_nodes_[idem];
    return *std::find_if(node_list.begin(), node_list.end(),
                         [marking](Node* node) { return node->type() == marking; });
  }
  
 public:
  /* General functions.
   */
  
  /* Returns the polynomial at the given root node. In practice, this is called
   * on the root node of a tree without arcs, i.e., a fully reduced D-module.
   */
  template< class DA_bimodule_it >
  Polynomial poincare_polynomial(DA_bimodule_it da_bimodule_it, const Idem& idem = Idem("0")) const {
    if (root_nodes_.find(idem) == root_nodes_.end()) {
      return Polynomial();
    }
    return root_nodes_.at(idem).poincare_polynomial(da_bimodule_it);
  }
  
  void dualize() {
    for (auto& arc : differential_arcs_) {
      arc.dualize();
    }
  }
  
  /* Pre-condition: nodes have no incoming or outgoing arcs */
  void post_construction_processing() {
    modulo_2_();
    connect_nodes_();
    connect_arcs_();
  }
  
 private:
  /* Modulo 2 needs to be called when adding an arc to differential_arcs, as
   * well as when adding a radient path to gradient paths/partial gradient path
   * to the incoming partial gradient paths of a tbd node...
   * We need to check if there exists an older arc with same alg el, and source
   * and target that are either more specific or more general.
   * 
   * The issue is that calculating neighbors live raises the complexity from
   * linear in #layers to linear in tree size.
   */
  
  /* Post-condition: node and arc neighbors are invalidated.
   * 
   * Cancel pairs of arcs with the same source and target.
   * This is probably the main reason that the container for arcs is a list,
   * and not a vector.
   * 
   * TO DO: Complete this function with partial cancelations
   */
  void modulo_2_() {
    raise_partially_overlapping_arcs_();  // partial cancellations
    
    bool canceled_pair;
    for (Arc_const_iterator arc_it1 = differential_arcs_.cbegin();
                            arc_it1 != differential_arcs_.cend(); ) {
      canceled_pair = false;
      for (Arc_const_iterator arc_it2 = std::next(arc_it1);
                              arc_it2 != differential_arcs_.cend(); ++arc_it2) {
        if (*arc_it1 == *arc_it2) {
          arc_it1 = differential_arcs_.erase(arc_it1);  // erase first arc
          // erase second arc
          if (arc_it1 == arc_it2) {
            arc_it1 = differential_arcs_.erase(arc_it2);  // update arc_it1
          }
          else {
            differential_arcs_.erase(arc_it2);
          }
          canceled_pair = true;
          break;  // exit arc_it2 loop
        }
      }  // arc_it2 loop
      if (!canceled_pair) {
        ++arc_it1;
      }
    }
  }  // modulo_2_
  
  /* Pre-condition: nodes have no neighbors */
  void connect_nodes_() {
    for (Arc_iterator arc_it = differential_arcs_.begin();
                      arc_it != differential_arcs_.end(); ++arc_it) {
      arc_it->source()->add_outgoing_arc_it(arc_it);
      arc_it->target()->add_incoming_arc_it(arc_it);
    }
  }
  
  /* Pre-condition: nodes are connected and arcs have no neighbors */
  void connect_arcs_() {
    for (Arc_iterator arc_it = differential_arcs_.begin();
                      arc_it != differential_arcs_.end(); ++arc_it) {
      arc_it->source()->calculate_arcs_neighboring_source(arc_it);
      arc_it->target()->calculate_arcs_neighboring_target(arc_it);
    }
  }
  
  void clear_neighbors_() {
    for (auto& value : root_nodes_) {
      value.second.clear_neighbors();
    }
    for (auto& arc : differential_arcs_) {
      arc.clear_neighbors();
    }
  }
  
  void recalculate_neighbors_() {  // clear and reconnect
    clear_neighbors_();  // O(#nodes)
    connect_nodes_();    // O(#arcs)
    connect_arcs_();     // O(#arcs#layers)
  }
  
  template< class DA_bimodule_it >
  void trim_leaves_(DA_bimodule_it da_bimodule_it) {
    int max_arc_ht = 0;
    for (auto& arc : differential_arcs_) {
      if (arc.ht() > max_arc_ht) {
        max_arc_ht = arc.ht();
      }
    }
    for (auto& value : root_nodes_) {
      value.second.trim_leaves(max_arc_ht, da_bimodule_it);
    }
  }
  
  /* TO DO: Redo this with something less than quadratic complexity... */  
  void raise_partially_overlapping_arcs_() {
    for (auto arc_it = differential_arcs_.begin();
              arc_it != differential_arcs_.end(); ++arc_it) {
      Node* current_source = arc_it->source()->parent();
      Node* current_target = arc_it->target()->parent();
      for (int rel_ht = 1; current_source != nullptr and current_target != nullptr; ++rel_ht) {
        for (auto lower_arc_it = differential_arcs_.begin();
                  lower_arc_it != differential_arcs_.end(); ) {
          if (lower_arc_it != arc_it and lower_arc_it->source() == current_source
              and lower_arc_it->target() == current_target
              and arc_it->algebra_element() == lower_arc_it->algebra_element()) {
            auto raised_arcs = lower_arc_it->raised_arcs(rel_ht);
            differential_arcs_.insert(differential_arcs_.end(), raised_arcs.begin(), raised_arcs.end());
            lower_arc_it = differential_arcs_.erase(lower_arc_it);
          }
          else {
            ++lower_arc_it;
          }
        }
        current_source = current_source->parent();
        current_target = current_target->parent();
      }
    }
  }
  
 public:  
  /* Discrete Morse theory */
  
  /* reduce.
   * Apply discrete Morse theory.
   * Pre-condition: nodes and arcs are connected.
   */
  template< class DA_bimodule_it >
  bool reduce(DA_bimodule_it da_bimodule_it, std::ofstream& write_file) {
#ifdef DEBUG
    int n_gen = 0;
    for (auto& value : root_nodes_) {
      n_gen += value.second.n_gen();
    }
    std::cout << "(old n_gen = " << n_gen << ")" << std::flush;
#endif  // DEBUG
    
    std::vector< Arc* > tbr_arcs = mark_reversible_arcs_();
    mark_critical_nodes_();
    
#ifdef DRAW
#ifdef DEBUG
    std::cout << "Draw..." << std::flush;
#endif  // DEBUG
    write_file << "Before reducing:" << std::endl;
    TeXify(write_file, da_bimodule_it);
    write_file << std::endl;
#endif  // DRAW
    
    extend_undecided_arcs_();
    modulo_2_();
    recalculate_neighbors_();  // because things are invalidated
    
    /* Discrete Morse theory */
    Arc_container gradient_paths = initialize_gradient_paths_();
    extend_gradient_paths_(gradient_paths, tbr_arcs);
    
    differential_arcs_.swap(gradient_paths);
    delete_nodes_();
    //trim_leaves_(da_bimodule_it);
    modulo_2_();
    
    //differential_arcs_.remove_if([](const Arc& arc) { return !arc.same_paths_exist(); });
    
    recalculate_neighbors_();
    
#ifdef DRAW   
#ifdef DEBUG
    std::cout << "ing... " << std::flush;
#endif  // DEBUG
    write_file << "After reducing:" << std::endl;
    TeXify(write_file, da_bimodule_it);
    write_file << "\n" << std::endl;
#endif  // DRAW
    
    bool done = true;
    for (auto& arc : differential_arcs_) {
      arc.reset_flags();
      if (arc.is_closed()) {
        done = false;
      }
    }
    return done;
  }
  
#ifdef DEBUG
  void debug() {
    int n_gen = 0;
    for (auto& value : root_nodes_) {
      n_gen += value.second.n_gen();
    }
    std::cout << "(n_gen = " << n_gen << ")" << std::flush;
    
    for (auto& value : root_nodes_) {
      value.second.leaves();
    }
    for (auto& arc : differential_arcs_) {
      arc.debug_out("[dsf] reduced");
    }
  }
#endif  // DEBUG
  
 private:  
  /* Consider the following arrows in the Hasse diagram associated to the
   * D-module, where each arrow points downwards, and the arrow drawn with # is
   * the one that interests us:
   * 
   *                 *
   *                /
   *               /
   *        *     *
   *       / \   # \
   *      /   \ #   \
   *     *     *     *
   *          /
   *         /
   *        *
   * 
   * If the # arrow is closed, and has not yet been marked as not-reversed,
   * then we mark it as to-be-reversed, and mark all of the adjacent arrows,
   * as drawn in the diagram, as not-reversed.
   * 
   * In reality, we work with differential arcs, but the idea is the same.
   */
  std::vector< Arc* > mark_reversible_arcs_() {
    std::vector< Arc* > to_be_reversed_arcs;
    for (auto& arc : differential_arcs_) {
      if (arc.is_closed() and !arc.not_reversed()) {
        arc.set_to_be_reversed();
        to_be_reversed_arcs.push_back(&arc);
        /* Incidentally, mark relevant leaves as to-be-deleted */
        //arc.mark_tbd_leaves();
        arc.source()->set_leaves_tbd();
        arc.target()->set_leaves_tbd();
        
        for (auto& arc : arc.others_to_source()) {
          arc.set_not_reversed();
        }
        for (auto& arc : arc.others_from_source()) {
          arc.set_not_reversed();
        }
        for (auto& arc : arc.others_from_target()) {
          arc.set_not_reversed();
        }
        for (auto& arc_zig : arc.others_to_target()) {
          arc_zig.set_not_reversed();
          for (auto& arc_zag : arc_zig.others_from_source()) {
            if (!arc_zag.to_be_reversed()) {
              arc_zag.set_not_reversed();
            }
          }
        }
        break;  // only reverse one arc
      }  // closed, not not-reversed arc
    }
    return to_be_reversed_arcs;
  }  // mark_reversible_arcs_
  
  /* Pre-condition: arcs have been marked as to be reversed, leaves tbd */
  void mark_critical_nodes_() {
    for (auto& value : root_nodes_) {
      value.second.set_criticalness();
    }
  }
  
  /* Post-condition: node and arc neighbors are invalidated
   */
  void extend_undecided_arcs_() {
    for (Arc_iterator arc_it = differential_arcs_.begin();
                      arc_it != differential_arcs_.end(); ) {
      // with this line, we basically get all coefficients
      //if (!arc_it->source()->children().empty() and !arc_it->target()->children().empty()) {
      if (arc_it->source()->criticalness() == 1 or arc_it->target()->criticalness() == 1) {
        auto ext_arcs = arc_it->raised_arcs(1);
        differential_arcs_.insert(differential_arcs_.end(), ext_arcs.begin(), ext_arcs.end());
        arc_it = differential_arcs_.erase(arc_it);
      }
      else {  // move to next arc
        ++arc_it;
      }
    }
  }
  
  /* Initialize gradient paths. Returns list of length 1 gradient paths, and
   * adds partial gradient paths to tbr arcs.
   */
  Arc_container initialize_gradient_paths_() {
    Arc_container gradient_paths;
    for (auto& arc : differential_arcs_) {
      if (arc.source()->criticalness() == 0) {
        if (arc.target()->criticalness() == 0) {
          gradient_paths.push_back(arc);
        }
        else {  // partial gradient path, target is tbd
          for (auto& arc2 : arc.others_to_target()) {
            if (arc2.to_be_reversed()) {  // add partial gradient path
              arc2.add_incoming_partial_gradient_path(arc);
            }
          }
        }  // partial gradient path
      }
    }
    return gradient_paths;
  }  // initialize_gradient_paths_
  
  /* Extend partial gradient paths.
   * We do this in the order in which the tbr arcs were added.
   * Note that a gradient path may pass by the same arc several times
   * consecutively, which is possible because we iterate over a list that may
   * call push_back.
   */
  void extend_gradient_paths_(Arc_container& gradient_paths, std::vector< Arc* >& to_be_reversed_arcs) {
    for (auto rev_arc : to_be_reversed_arcs) {
      rev_arc->reverse();
    }
    for (auto rev_arc : to_be_reversed_arcs) {
      //rev_arc->debug_out("[dmt] rev");
      for (auto& gradient_path : rev_arc->gradient_paths()) {  // incoming partial gp
        //gradient_path.debug_out("[dmt] pre");
        if (!gradient_path.compatible_above(*rev_arc)) { continue; }
        Arc arc_zig = gradient_path.concatenate(*rev_arc);
        //arc_zig.debug_out("[dmt] zig");
        for (auto& post_arc : rev_arc->others_from_target()) {  // rev_arc has been reversed
          //post_arc.debug_out("[dmt] post");
          if (!arc_zig.compatible_above(post_arc, rev_arc->ht())) { continue; }
          if (gradient_path.source_idem().too_far_from(post_arc.target_idem())) { continue; }
          Arc arc_zigzag = arc_zig.concatenate(post_arc);
          if (arc_zigzag.algebra_element().is_null()) { continue; }
          //arc_zigzag.debug_out("[dmt] zigzag");
          
          if (arc_zigzag.target()->criticalness() == 0) {
            arc_zigzag.gradient = true;
            gradient_paths.push_back(arc_zigzag);
          }
          else {
            for (auto& next_arc : post_arc.others_to_target()) {  // next_arc may == rev_arc
              if (next_arc.to_be_reversed() and arc_zigzag.compatible_above(next_arc, post_arc.ht())) {
                //next_arc.debug_out("[dmt] next");
                next_arc.add_incoming_partial_gradient_path(arc_zigzag);
              }
            }
          }
        }  // post_arc
      }
    }
  }  // extend_gradient_paths_
  
  /* Delete nodes marked as to-be-deleted */
  void delete_nodes_() {
    for (auto map_it = root_nodes_.begin(); map_it != root_nodes_.end(); ) {
      if (map_it->second.criticalness() == 2) {
        map_it = root_nodes_.erase(map_it);
      }
      else {
        map_it->second.erase_children();
        ++map_it;
      }
    }
  }
    
 public:
  /* TeXify */
  
  /* Note: the DA-bimodules are ordered top-down, but we build the grid
   * bottom-up.
   */
  template< class DA_bimodule_it >
  void TeXify(std::ofstream& write_file, DA_bimodule_it da_bimodule_it) {
    Grid_ grid;
    
    for (auto& value : root_nodes_) {
      add_to_grid_(grid, 0, &value.second);
    }
    
    write_file << "\\begin{tikzpicture}[suffix forest]" << std::endl;
    
    std::string extra_options;
    
    ++da_bimodule_it;  // we never use the iterator if it goes past the end
    for (int layer = 0; layer < grid.size(); ++layer) {
      for (Grid_point_& grid_point : grid[layer]) {
        if (grid_point.first->criticalness() == 0) {
          extra_options = "\\color{blue}";
        }
        else if (grid_point.first->criticalness() == 1) {
          extra_options = "\\color{teal}";
        }
        else {
          extra_options = "\\color{red}";
        }
        write_file << "\\node ("
                   << static_cast<void*>(grid_point.first)
                   << ") at (-"
                   << grid_point.second
                   << ","
                   << static_cast<float>(layer) * y_sep_
                   << ") {" << extra_options << "\\nodeLabel{"
                   << grid_point.first->idem().to_string()
                   << "}};" << std::endl;
        if (grid_point.first->has_parent()) {
          write_file << "\\draw[->] ("
                     << static_cast<void*>(grid_point.first)
                     << ") -- node[in place]{$"
                     << (*da_bimodule_it)->get_label(grid_point.first->type())
                     << "$} ("
                     << static_cast<void*>(grid_point.first->parent())
                     << ");" << std::endl;
        }
      }
      --da_bimodule_it;
    }
    if (!grid.empty()) {
      float y_poly = static_cast<float>(grid.size() - 1) * y_sep_ + poly_sep_;
      for (Grid_point_& grid_point : grid.back()) {
        write_file << "\\node at (-"
                   << grid_point.second
                   << ","
                   << y_poly
                   << ") {\\footnotesize$"
                   << grid_point.first->poincare_polynomial().to_string()
                   << ", "
                   << grid_point.first->weights().first
                   << ", "
                   << grid_point.first->weights().second
                   << "$};" << std::endl;
      }
    }
    
    for (auto& arc : differential_arcs_) {
      //if (!arc.is_closed()) { continue; }
      if (arc.to_be_reversed()) {
        extra_options = "\\color{red}";
      }
      else if (arc.not_reversed()) {
        extra_options = "\\color{blue}";
      }
      if (arc.gradient) {
        extra_options = "\\color{green}";
      }
      write_file << "\\path ("
                 << static_cast<void*>(arc.source())
                 << ") edge[differential arc, bend left=10] node[in place]{$"
                 << extra_options << arc.algebra_string()
                 << "$} ("
                 << static_cast<void*>(arc.target())
                 << ");" << std::endl;
    }
    write_file << "\\end{tikzpicture}" << std::flush;
  }
  
 private:
  typedef std::pair< Node*, float > Grid_point_;
  typedef std::vector< Grid_point_ > Grid_layer_;
  typedef std::vector< Grid_layer_ > Grid_;
  
  static constexpr float min_x_sep_ = 1.;  // minimal horizontal distance between nodes
  static constexpr float y_sep_ = .8;  // vertical distance between nodes
  static constexpr float poly_sep_ = .3;  // vertical distance between polynomial and leaf
  
  void add_to_grid_(Grid_& grid, int layer, Node* node_ptr) {
    float x;
    if (node_ptr->children_begin() == node_ptr->children_end()) {  // no children
      while (grid.size() <= layer) {  // make space for point
        grid.push_back(Grid_layer_());
      }
      if (grid[layer].size() == 0) {
        x = 0.;
      }
      else {
        x = grid[layer].back().second + min_x_sep_;
      }
      grid[layer].push_back(Grid_point_(node_ptr, x));
    }  // no children
    else {  // with children
      while (grid.size() <= layer + 1) {  // make space for children
        grid.push_back(Grid_layer_());
      }
      int first_child_index = grid[layer + 1].size();
      
      for (auto child_it = node_ptr->children_begin();
                child_it != node_ptr->children_end(); ++child_it) {
        add_to_grid_(grid, layer + 1, &*child_it);
      }
      x = (grid[layer + 1][first_child_index].second + grid[layer + 1].back().second) / 2.;
      grid[layer].push_back(Grid_point_(node_ptr, x));
    }  // with children
  }
  
  /* Private members */
  Node_container root_nodes_;
  Arc_container differential_arcs_;
  
  std::map< const Idem, std::vector< Node* > > layer_one_nodes_;  // for \gamma_1
  
  /* Discrete Morse theory */
  std::vector< Arc* > to_be_reversed_arcs_;
  Arc_container gradient_paths_;
};

#endif  // DIFFERENTIAL_SUFFIX_FOREST_H_