#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>  // shared_ptr
#include <set>
#include <string>
#include <vector>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include "Differential_suffix_forest/Differential_suffix_forest_options.h"

#include "test_arc.h"
#include "test_node.h"
#include "test_utils.h"

/* Differential suffix forest.
 * 
 * Arcs are stored in a boost multi-index container, indexed by source and
 * target nodes. Thus, we can iterate over arcs in order of source or target.
 * This is useful for finding neighbors without too much overhead (maybe).
 */
class Forest {
 public:
  using Idem = typename Forest_options_default_short::Idem;
  using Node = Node< Forest_options_default_short >;
  using Arc = Arc< Forest_options_default_short >;
  
  using Gen_type = typename Node::Gen_type;
  using Alg_el = typename Arc::Coefficient_value;
  
  struct Source { };
  struct Target { };
  
  struct Extract_source {
    using result_type = int;
    
    result_type operator()(Arc arc) const {
      return arc.source();
    }
  };
  
  struct Extract_target {
    using result_type = int;
    
    result_type operator()(Arc arc) const {
      return arc.target();
    }
  };
  
  using Arc_container = boost::multi_index_container<
    Arc,
    boost::multi_index::indexed_by<
      boost::multi_index::ordered_non_unique< boost::multi_index::tag< Source >, Extract_source >,
      boost::multi_index::ordered_non_unique< boost::multi_index::tag< Target >, Extract_target >
    >
  >;
  using Arc_source_container = typename Arc_container::index<Source>::type;
  using Arc_target_container = typename Arc_container::index<Target>::type;
  using Arc_source_iterator = typename Arc_source_container::iterator;
  using Arc_target_iterator = typename Arc_target_container::iterator;
  
  std::vector< Node > nodes() const {
    return nodes_;
  }
  
  Arc_container arcs() const {
    return arcs_;
  }
  
  /* Operations on nodes */
  
  /* Declare subtree to a specified idempotent. The extra argument at the end
   * is to force the user to acknowledge from which forest the old root comes.
   * Pre-condition: subtrees unlocked.
   */
  void declare_subtree(Idem new_idem, Gen_type new_type, int old_root, Forest&) {
    declared_subtrees_[new_idem].emplace_back(new_type, old_root);
  }
  
  void declare_subtree(Idem new_idem) {
    declared_subtrees_[new_idem];
  }
  
  /* Lock subtrees using another forest.
   * Post-condition: subtrees locked.
   */
  void lock_subtrees(const Forest& old_forest) {
    nodes_.clear();
    for (auto& map_value : declared_subtrees_) {
      int new_root = nodes_.size();
      const Idem& new_idem = map_value.first;
      roots_[new_idem] = new_root;
      nodes_.push_back(Node());  // add new root node
      
      for (auto& vec_value : map_value.second) {
        Gen_type& new_type = vec_value.first;
        int old_root = vec_value.second;
        int d_after_children = old_forest.nodes_[old_root].d_after_children();
        int new_child = nodes_.size();
        nodes_.insert(nodes_.end(),
                      old_forest.nodes_.begin() + old_root,
                      old_forest.nodes_.begin() + old_root + d_after_children);
        nodes_[new_child].set_d_parent(new_child - new_root);
        nodes_[new_child].set_type(new_type);
        nodes_[new_root].add_d_after_children(d_after_children);
      }
    }
    declared_subtrees_.clear();
  }
  
  /* Lock subtrees, but only the roots */
  void lock_subtrees() {
    nodes_.clear();
    for (auto& map_value : declared_subtrees_) {
      int new_root = nodes_.size();
      const Idem& new_idem = map_value.first;
      roots_[new_idem] = new_root;
      nodes_.push_back(Node());
    }
  }
  
  /* Find and return child node with a certain type to parent. Only used
   * internally when declaring arcs.
   */
  int child_of_type_(int node, Gen_type type) const {
    for (int child = node + 1; child < node + nodes_[node].d_after_children(); child += nodes_[child].d_after_children()) {
      if (nodes_[child].type() == type) {
        return child;
      }
    }
    return -1;
  }
  
  /* Distance from a node to its root. Only used internally when declaring
   * arcs.
   */
  int d_root_(int node) const {
    int distance = 0;
    while (nodes_[node].d_parent() > 0) {
      distance += nodes_[node].d_parent();
      node -= nodes_[node].d_parent();
    }
    return distance;
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
  
  void declare_arc(Alg_el new_value, Gen_type back_marking, Gen_type front_marking, const Arc& old_arc, const Forest& old_forest) {
    int source = child_of_type_(roots_[new_value.source_idem()], back_marking);
    int target = child_of_type_(roots_[new_value.target_idem()], front_marking);
    source += old_forest.d_root_(old_arc.source());
    target += old_forest.d_root_(old_arc.target());
    declared_arcs_.emplace_back(source, target, new_value);
  }
  
  void declare_arc(Alg_el new_value, Gen_type back_marking, Gen_type front_marking, const Idem& old_idem) {
    int source = child_of_type_(roots_[new_value.source_idem()], back_marking);
    int target = child_of_type_(roots_[new_value.target_idem()], front_marking);
    declared_arcs_.emplace_back(source, target, new_value);
  }
  
  void declare_arc(int source, int target, Alg_el value) {
    declared_arcs_.emplace_back(source, target, value);
  }
  
  void lock_arcs() {
    arcs_.insert(declared_arcs_.begin(), declared_arcs_.end());
    modulo_2();
  }
  
 private:
  void modulo_2() {
    Arc_source_container& arcs_source_view = arcs_.get<Source>();
    for (auto lower_arc_it = arcs_.begin(); lower_arc_it != arcs_.end(); ) {
      lower_arc_it = resolve_overlapping_arcs_(lower_arc_it, arcs_source_view);
    }
  }
  
  /* Raise partially overlapping arcs (+ cancel exactly overlapping arcs)
   * 
   * Scan arcs above given arc. If one is found, mark its ancestors and move on
   * to the first arc after children.
   * 
   * The places to add a new arc are the unmarked nodes whose parent is marked.
   * 
   * Return the iterator of the next arc.
   * 
   * Modifies arcs_.
   */
  Arc_source_iterator resolve_overlapping_arcs_(Arc_source_iterator lower_arc_it, Arc_source_container& arcs_source_view) {
    auto upper_arc_it = std::next(lower_arc_it);
    
    // Scan for arcs that cancel completely
    for (; upper_arc_it != arcs_source_view.end() and upper_arc_it->source() == lower_arc_it->source(); ++upper_arc_it) {
      if (*upper_arc_it == *lower_arc_it) {
        arcs_source_view.erase(upper_arc_it);
        lower_arc_it = arcs_source_view.erase(lower_arc_it);
        return lower_arc_it;
      }
    }
    
    // Scan descendants
    int start_node = lower_arc_it->source();
    int end_node = lower_arc_it->source() + nodes_[lower_arc_it->source()].d_after_children();
    std::vector< bool > marked(nodes_[lower_arc_it->source()].d_after_children(), false);
    
    for (; upper_arc_it != arcs_source_view.end() and upper_arc_it->source() < end_node; ) {
      if (upper_arc_it->value() == lower_arc_it->value()) {
        // Mark ascendants
        int d_parent = nodes_[upper_arc_it->source()].d_parent();
        int parent = upper_arc_it->source() - d_parent;
        while (d_parent > 0 and parent >= lower_arc_it->source()) {
          marked[parent - start_node] = true;
          d_parent = nodes_[parent].d_parent();
          parent = parent - d_parent;
        }
        
        // Skip arcs until after children
        int after_children = upper_arc_it->source() + nodes_[upper_arc_it->source()].d_after_children();
        for (; upper_arc_it != arcs_source_view.end() and upper_arc_it->source() < after_children; ++upper_arc_it);
      }
      else {
        ++upper_arc_it;
      }
    }
    
    // Place new arcs; upper_arc_it is now a hint for placement
    for (int node = start_node + 1; node != end_node; ++node) {
      if (!marked[node] and marked[node - nodes_[node].d_parent()]) {
        arcs_source_view.emplace_hint(upper_arc_it, node, lower_arc_it->target() + node - lower_arc_it->source(), lower_arc_it->value());
      }
    }
    
    if (marked[0]) {
      lower_arc_it = arcs_source_view.erase(lower_arc_it);
      return lower_arc_it;
    }
    else {
      return std::next(lower_arc_it);
    }
  }
  
  
 public:
  /* Views on arcs */
  
  /* Produces a composite range.
   * 
   * This might be costly, because composite ranges are iffy.
   */
  Composite_range< Arc_target_container > arcs_to_source(Arc arc) {
    Composite_range< Arc_target_container > comp_range;
    int source = arc.source();
    int after_children = source + nodes_[source].d_after_children();
    
    const Arc_target_container& arcs_target_view = arcs_.get<Target>();
    Arc_target_iterator arc_target_begin = arcs_target_view.lower_bound(source);
    Arc_target_iterator arc_target_end = arcs_target_view.lower_bound(after_children);
    comp_range.add_range(arc_target_begin, arc_target_end);
    
    for (int d_parent = nodes_[source].d_parent(); d_parent > 0; ) {
      source = source - d_parent;
      d_parent = nodes_[source].d_parent();
      
      auto it_pair = arcs_target_view.equal_range(source);
      comp_range.add_range(it_pair);
    }
    return comp_range;
  }
  
  /* Operations on arcs
   * 
   * These are the elementary operations we need: testing for compatibility,
   * concatenating, and concatenating zig-zags.
   */
  
  bool compatible(const Arc& back_arc, const Arc& front_arc) const {
    return (back_arc.target() <= front_arc.source()
            and front_arc.source() < back_arc.target() + nodes_[back_arc.target()].d_after_children());
  }
  
  Arc concatenate(const Arc& back_arc, const Arc& front_arc) const {
    int difference = front_arc.source() - back_arc.target();
    if (difference >= 0) {
      return Arc(back_arc.source() + difference, front_arc.target(), back_arc.value() * front_arc.value());
    }
    else {
      return Arc(back_arc.source(), front_arc.target() - difference, back_arc.value() * front_arc.value());
    }
  }
  
  /* TeXify */
  
  void TeXify(std::ofstream& write_file) {
    Grid_ grid;
    
    for (auto value_pair : roots_) {
      add_to_grid_(grid, 0, value_pair.second);
    }
    
    write_file << "\\begin{tikzpicture}[suffix forest]" << std::endl;
    
    std::string extra_options("");
    
    /* Only put idempotents for the roots */
    int i = 0;
    for (auto value_pair : roots_) {
      Grid_point_& grid_point = grid[0][i];
      write_file << "\\node ("
                   << grid_point.first
                   << ") at ("
                   << grid_point.second
                   << ","
                   << 0
                   << ") {" << extra_options << "\\nodeLabel{"
                   << value_pair.first
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
                   << static_cast<int>(nodes_[grid_point.first].type())
                   << "$} ("
                   << (grid_point.first - nodes_[grid_point.first].d_parent())
                   << ");" << std::endl;
      }
    }
    
    for (auto& arc : arcs_) {
      write_file << "\\path ("
                 << arc.source()
                 << ") edge[differential arc, bend left=10] node[in place]{$"
                 << extra_options << arc.value()
                 << "$} ("
                 << arc.target()
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
  static constexpr float poly_sep_ = .3;  // vertical distance between polynomial and leaf
  
  void add_to_grid_(Grid_& grid, int layer, int node) {
    float x;
    if (nodes_[node].no_children()) {  // no children
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
      
      for (int child = node + 1; child < node + nodes_[node].d_after_children(); child += nodes_[child].d_after_children()) {
        add_to_grid_(grid, layer + 1, child);
      }
      x = (grid[layer + 1][first_child_index].second + grid[layer + 1].back().second) / 2.;
      grid[layer].push_back(Grid_point_(node, x));
    }  // with children
  }
  
  std::vector< Node > nodes_;
  
  std::vector< Arc > declared_arcs_;
  Arc_container arcs_;
    
  std::map< const Idem, std::vector< std::pair< Gen_type, int > > > declared_subtrees_;
  
  std::map< const Idem, int > roots_;
};

std::ostream& operator<<(std::ostream& os, const Forest& forest) {
  os << "Forest nodes: " << std::flush;
  for (auto& node : forest.nodes()) {
    os << node;
  }
  os << std::endl;
  os << "Forest arcs: " << std::flush;
  for (auto& arc : forest.arcs()) {
    os << arc;
  }
  os << std::endl;
  return os;
}