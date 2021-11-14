#include <fstream>
#include <functional>  // reference_wrapper
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>

#include "Differential_suffix_forest/Differential_suffix_forest_options.h"

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
  using Alg_el = typename Forest_options::Alg_el;
    
  struct Node {
    /* Constructor for a root node */
    Node() :
      d_parent(0),
      d_after_children(1),
      type(-1)
    { }
    
    /* Constructor when we know the parent, and when there are no arcs yet.
     * Initializing next_sibling_ is a bit iffy.
     */
    Node(const int d) :
      d_parent(d),
      d_after_children(1),
      type(-1)
    { }
    
    int d_parent;
    int d_after_children;
    Gen_type type;
    
    /* Evaluation functions */
    bool has_parent() const {
      return d_parent > 0;
    }
    
    bool no_children() const {
      return d_after_children == 1;
    }
    
    friend std::ostream& operator<<(std::ostream& os, const Node& node) {
      os << "<" << node.d_parent << "|" << static_cast<int>(node.type) << "|" << node.d_after_children << "> " << std::flush;
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
      return (source == other.source and target == other.target and value == other.value);
    }
    
    bool operator!=(const Arc& other) const {
      return !(*this == other);
    }
    
    friend std::ostream& operator<<(std::ostream& os, const Arc& arc) {
      os << "(" << arc.source << "|" << arc.value << "|" << arc.target << ")" << std::flush;
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
  
  using Arc_source_container = typename Arc_container::template index< Source >::type;
  using Arc_target_container = typename Arc_container::template index< Target >::type;
  using Arc_source_iterator = typename Arc_source_container::iterator;
  using Arc_target_iterator = typename Arc_target_container::iterator;
  
  using Arc_view = Arc_source_container;
  using Arc_iterator = Arc_source_iterator;
  
  using Arc_reference = std::reference_wrapper< const Arc >;
  
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
      add_root_(new_idem, new_root);
      nodes_.push_back(Node());  // add new root node
      
      for (auto& vec_value : map_value.second) {
        Gen_type& new_type = vec_value.first;
        int old_root = vec_value.second;
        int d_after_children = old_forest.nodes_[old_root].d_after_children;
        int new_child = nodes_.size();
        nodes_.insert(nodes_.end(),
                      old_forest.nodes_.begin() + old_root,
                      old_forest.nodes_.begin() + old_root + d_after_children);
        nodes_[new_child].d_parent = new_child - new_root;
        nodes_[new_child].type = new_type;
        nodes_[new_root].d_after_children += d_after_children;
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
      add_root_(new_idem, new_root);
      nodes_.push_back(Node());
    }
  }
  
 private:
  void add_root_(const Idem& idem, int root) {
    roots_[idem] = root;
    root_indices_.insert(root);
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
    source += old_forest.d_root_(old_arc.source);
    target += old_forest.d_root_(old_arc.target);
    declared_arcs_.emplace_back(source, target, new_value);
  }
  
  void declare_arc(Alg_el new_value, Gen_type back_marking, Gen_type front_marking, const Idem& old_idem) {
    int source = child_of_type_(roots_[new_value.source_idem()], back_marking);
    int target = child_of_type_(roots_[new_value.target_idem()], front_marking);
    declared_arcs_.emplace_back(source, target, new_value);
  }
  
  /* Only for debugging purposes */
  void declare_arc(int source, int target, Alg_el value) {
    declared_arcs_.emplace_back(source, target, value);
  }
  
  void lock_arcs() {
    arcs_.insert(declared_arcs_.begin(), declared_arcs_.end());
    modulo_2_();
  }
  
 private:
  /* Find and return child node with a certain type to parent. Only used
   * internally when declaring arcs.
   */
  int child_of_type_(int node, Gen_type type) const {
    for (int child = node + 1; child != node + nodes_[node].d_after_children; child += nodes_[child].d_after_children) {
      if (nodes_[child].type == type) {
        return child;
      }
    }
    return -1;
  }
  
  int root_of_(int node) const {
    return *--root_indices_.upper_bound(node);
  }
  
  /* Distance from a node to its root. Only used internally when declaring
   * arcs.
   */
  int d_root_(int node) const {
    return node - root_of_(node);
  }
  
  
  void modulo_2_() {
    for (auto lower_arc_it = arcs_.begin(); lower_arc_it != arcs_.end(); ) {
      lower_arc_it = resolve_overlaps_after_(arcs_, lower_arc_it);
    }
  }
  
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
   * 
   * Modifies the container viewed by arcs_view.
   * 
   * Templated by the view of the arcs. In practice, this can be an index of
   * the multi-indexed container of arcs.
   * 
   * WARNING: currently hard-coded for source-view... probably just keep this
   */
  Arc_iterator resolve_overlaps_after_(Arc_view& arcs_view, Arc_iterator lower_arc_it) {
    auto upper_arc_it = std::next(lower_arc_it);
    
    // Scan for arcs that cancel completely
    for (; upper_arc_it != arcs_view.end() and upper_arc_it->source == lower_arc_it->source; ++upper_arc_it) {
      if (*upper_arc_it == *lower_arc_it) {
        arcs_view.erase(upper_arc_it);
        lower_arc_it = arcs_view.erase(lower_arc_it);
        return lower_arc_it;
      }
    }
    
    // Scan descendants
    int start_node = lower_arc_it->source;
    int n_nodes = nodes_[start_node].d_after_children;
    int end_node = start_node + n_nodes;
    std::vector< bool > marked(n_nodes, false);
    
    for (; upper_arc_it != arcs_view.end() and upper_arc_it->source < end_node; ) {
      if (overlap_(*lower_arc_it, *upper_arc_it)) {
        // Mark ascendants
        int d_parent = nodes_[upper_arc_it->source].d_parent;
        int parent = upper_arc_it->source - d_parent;
        while (d_parent > 0 and parent >= start_node) {
          marked[parent - start_node] = true;
          d_parent = nodes_[parent].d_parent;
          parent = parent - d_parent;
        }
        
        // Skip arcs until after children
        int after_children = upper_arc_it->source + nodes_[upper_arc_it->source].d_after_children;
        for (; upper_arc_it != arcs_view.end() and upper_arc_it->source < after_children; ++upper_arc_it);
      }
      else {
        ++upper_arc_it;
      }
    }
    
    // Place new arcs; upper_arc_it is now a hint for placement
    raise_arcs_above_hint_(arcs_view, upper_arc_it, *lower_arc_it, marked, start_node, 1, n_nodes);
    
    if (marked[0]) {  // arcs were raised
      lower_arc_it = arcs_view.erase(lower_arc_it);
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
  bool resolve_overlaps_before_(Arc_view& arcs_view, Arc_iterator upper_arc_it) {
    int upper_node = upper_arc_it->source;
    std::vector< int > ancestors = { upper_node };
    
    int d_parent = nodes_[upper_node].d_parent;
    int current_node = upper_node;
    auto lower_arc_it = upper_arc_it;
    do {
      --lower_arc_it;
      int endpoint = lower_arc_it->source;
      if (endpoint == current_node and overlap_(*lower_arc_it, *upper_arc_it)) {
        ancestors.push_back(current_node);
        raise_arcs_below_hint_(arcs_view, upper_arc_it, *lower_arc_it, ancestors);
        arcs_view.erase(lower_arc_it);
        arcs_view.erase(upper_arc_it);
        return true;
      }
      else if (endpoint < current_node) {
        do {
          ancestors.push_back(current_node);
          current_node -= d_parent;
          d_parent = nodes_[current_node].d_parent;
        } while (d_parent != 0 and endpoint < current_node);
      }
    } while (d_parent != 0 and lower_arc_it != arcs_view.begin());
    
    return false;
  }
  
  /* Given a critical arc, raise every arc below it to get critical and non-
   * critical arcs.
   * 
   * A critical arc is one that will be deleted. The input critical arc is
   * generally an invertible arc from the reduction process.
   * 
   * This needs to be done for source and target...
   */
  void resolve_criticalness_(Arc_view& arcs_view, Arc_iterator critical_arc_it) {
  }
  
  /* Raise original arc to avoid all marked nodes. The vector marked is
   * indexed relatively to a start node.
   * 
   * This function requires a hint for emplacement of new arcs. The new arcs
   * will be placed before the hint.
   */
  void raise_arcs_above_hint_(Arc_view& arcs_view,
                              const Arc_iterator& hint,
                              const Arc& original_arc,
                              const std::vector< bool >& marked,
                              const int start_node,
                              const int rel_begin,
                              const int rel_end) {
    for (int rel_node = rel_begin; rel_node != rel_end; ) {  // relative node index
      if (!marked[rel_node] and marked[rel_node - nodes_[start_node + rel_node].d_parent]) {
        arcs_view.emplace_hint(hint, start_node + rel_node, original_arc.target + rel_node, original_arc.value);
        rel_node += nodes_[start_node + rel_node].d_after_children;
      }
      else {
        ++rel_node;
      }
    }
  }
  
  /* ancestors holds source node and all of its ancestors up to the lower
   * overlapping arc.
   */
  void raise_arcs_below_hint_(Arc_view& arcs_view,
                              const Arc_iterator hint,
                              const Arc& original_arc,
                              const std::vector< int >& ancestors) {
    for (int i = 1; i != ancestors.size(); ++i) {
      int parent = ancestors[i];
      for (int child = parent + 1; child != parent + nodes_[parent].d_after_children; child += nodes_[child].d_after_children) {
        if (child != ancestors[i - 1]) {
          arcs_view.emplace_hint(hint, child, original_arc.target + child - original_arc.source, original_arc.value);
        }
      }
    }
  }
  
 public:
  /* Views on arcs */
  
  std::vector< Arc_reference > arcs_to_source(const Arc& arc) const {
    return arcs_at_node_< Target >(arc.source);
  }
  
  std::vector< Arc_reference > arcs_from_target(const Arc& arc) const {
    return arcs_at_node_< Source >(arc.target);
  }
  
  std::vector< Arc_reference > arcs_from_source(const Arc& arc) const {
    return arcs_at_node_< Source >(arc.source, arc);
  }
  
  std::vector< Arc_reference > arcs_to_target(const Arc& arc) const {
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
    int after_children = node + nodes_[node].d_after_children;
    const auto& arcs_view = arcs_.template get< Tag >();
    auto arc_begin = arcs_view.lower_bound(node);
    auto arc_end = arcs_view.lower_bound(after_children);
    
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
    int after_children = node + nodes_[node].d_after_children;
    const auto& arcs_view = arcs_.template get< Tag >();
    auto arc_begin = arcs_view.lower_bound(node);
    auto arc_end = arcs_view.lower_bound(after_children);
    
    std::vector< Arc_reference > result;
    
    for (auto arc_it = arc_begin; arc_it != arc_end; ++arc_it) {
      if (*arc_it != avoiding) {  // it sucks that we have to compare data structures instead of pointers...
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
    int d_parent = nodes_[node].d_parent;
    int current_node = node - d_parent;
    do {
      --arc_it;
      int endpoint = get_endpoint(*arc_it);
      if (endpoint == current_node) {
        arc_stream.emplace_back(*arc_it);
      }
      else if (endpoint < current_node) {
        do {
          d_parent = nodes_[current_node].d_parent;
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
             and front_arc.source < back_arc.target + nodes_[back_arc.target].d_after_children)
         or (front_arc.source <= back_arc.target
             and back_arc.target < front_arc.source + nodes_[front_arc.source].d_after_children));
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
  
  void reduce() {
    // find invertible arc
    // determine criticalness of all arcs whose source or target are below
    //   this arc
    
    for (const auto& arc : arcs_) {
      if (arc.value.is_invertible()) {
        contract_(arc);
      }
    }
  }
  
 private:
  /* Contract an invertible arc.
   */
  void contract_(const Arc& reverse_arc) {
    std::vector< Arc > zigzag_arcs;
    
    for (const Arc& back_arc : arcs_to_target(reverse_arc)) {
      for (const Arc& front_arc : arcs_from_source(reverse_arc)) {
        std::cout << "Making zig-zag arc from " << back_arc << " " << reverse_arc << " " << front_arc << std::endl;
        zigzag_arcs = add_zigzag_(zigzag_arcs, back_arc, reverse_arc, front_arc);
      }
    }
    
    for (const Arc& zigzag_arc : zigzag_arcs) {
      std::cout << "Adding zig-zag arc " << zigzag_arc << std::endl;
      insert_arc_(zigzag_arc);
    }
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
    
    if (back_diff <= front_diff and front_diff < back_diff + nodes_[back_arc.target].d_after_children) {
      source += front_diff - back_diff;
      arc_stream.emplace_back(source, target, back_arc.value * front_arc.value);
    }
    else if (front_diff <= back_diff and back_diff < front_diff + nodes_[front_arc.source].d_after_children) {
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
    bool cancel_below = resolve_overlaps_before_(arcs_, arc_it);
    if (!cancel_below) {
      std::cout << "checking above" << std::endl;
      resolve_overlaps_after_(arcs_, arc_it);
    }
  }
  
 public:
  /* I/O interface */
  
  friend std::ostream& operator<<(std::ostream& os, const Forest& forest) {
    os << "Forest nodes: " << std::flush;
    for (auto& node : forest.nodes()) {
      os << node;
    }
    os << std::endl;
    os << "Forest arcs: " << std::flush;
    for (auto& arc : forest.arcs()) {
      os << arc << " ";
    }
    os << std::endl;
    return os;
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
                   << static_cast<int>(nodes_[grid_point.first].type)
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
      
      for (int child = node + 1; child < node + nodes_[node].d_after_children; child += nodes_[child].d_after_children) {
        add_to_grid_(grid, layer + 1, child);
      }
      x = (grid[layer + 1][first_child_index].second + grid[layer + 1].back().second) / 2.;
      grid[layer].push_back(Grid_point_(node, x));
    }  // with children
  }
  
  std::vector< Node > nodes_;
  Arc_container arcs_;
  
  std::map< const Idem, std::vector< std::pair< Gen_type, int > > > declared_subtrees_;
  std::vector< Arc > declared_arcs_;
  
  std::map< const Idem, int > roots_;
  std::set< int > root_indices_;
};