#ifndef TEST_NODE_CONTAINER_H_
#define TEST_NODE_CONTAINER_H_

#include <map>
#include <string>
#include <vector>

#include "test_forest_options.h"

/* Nodes and roots.
 * 
 * Mathematically, this is a (bundled) left module over a bordered algebra.
 */
template< class Forest_options = Forest_options_default_short >
class Node_container {
 public:
  using Idem = typename Forest_options::Idem;
  using Weights = typename Forest_options::Weights;
  
  using Root_handle_container = std::map< int, Idem >;
  using Root_handle = typename Root_handle_container::value_type;
  
  struct Node {
    Node() :
      to_parent(0),
      to_next(1),
      descendants_size(1),
      weights(0, 0)
    { }
    
    Node(const int to_parent,
         const int to_next,
         const int descendants_size,
         const Weights weights
#ifdef DRAW
         , std::string label
#endif  // DRAW
         ) :
      to_parent(to_parent),
      to_next(to_next),
      descendants_size(descendants_size),
      weights(weights)
#ifdef DRAW
      , label(label)
#endif  // DRAW
    { }
    
    int to_parent;
    int to_next;
    int descendants_size;
    Weights weights;
#ifdef DRAW
    std::string label;
#endif  // DRAW
        
    friend std::ostream& operator<<(std::ostream& os, const Node& node) {
      os << "<"
         << node.to_parent
         << "|"
#ifdef DRAW
         << node.label
         << "|"
#endif  // DRAW
         << node.descendants_size
         << ">";
      return os;
    }
  };
  
  /* A sort-of iterator for parents.
   * 
   * This should satisfy LegacyIterator requirement.
   */
  class Ascender {
   public:
    Ascender(int node, const Node_container& nc) :
      node_(node),
      to_node_(1),
      nc_(nc)
    { }
    
    /* "indirection" */
    int operator*() const {
      return node_;
    }
    
    /* "increment" */
    Ascender& operator++() {
      to_node_ = nc_.to_parent(node_);
      node_ -= to_node_;
      return *this;
    }
    
    /* " != end " */
    bool valid() const {
      return to_node_ != 0;
    }
    
   private:
    int node_;
    int to_node_;
    const Node_container& nc_;
  };
  
  Node_container() { }
  
  /* OBSERVERS */
  
  /* Getters */
  
  Node operator[](const int i) const {
    return nodes_[i];
  }
  
  int to_next(int i) const {
    return nodes_[i].to_next;
  }
  
  int to_parent(int i) const {
    return nodes_[i].to_parent;
  }
  
  int descendants_size(int i) const {
    return nodes_[i].descendants_size;
  }
  
  Weights weights(int i) const {
    return nodes_[i].weights;
  }
  
#ifdef DRAW
  std::string label(int i) const {
    return nodes_[i].label;
  }
#endif  // DRAW
  
  int size() const {
    return nodes_.size();
  }
  
  Idem idem(const Root_handle& root_handle) const {
    return root_handle.second;
  }
  
  /* Other observers (constant cost) */
  
  int next(int i) const {
    return i + to_next(i);
  }
  
  int parent(int i) const {
    return i - to_parent(i);
  }
  
  bool is_root(int i) const {
    return to_parent(i) == 0;
  }
  
  bool is_first_child(int i) const {
    return to_parent(i) == to_next(parent(i));
  }
  
  // alias for next
  int descendants_begin(int i) const {
    return next(i);
  }
  
  int descendants_end(int i) const {
    return i + descendants_size(i);
  }
  
  bool has_children(int i) const {
    return to_next(i) != descendants_size(i);
  }
  
  Ascender ascender(int node) const {
    return Ascender(node, *this);
  }
  
  /* Other observers, any cost */
  
  int root(int i) const {
    return (--root_idems_.upper_bound(i))->first;
  }
  
  /* Distance from a node to its root. Only used internally when declaring
   * arcs.
   */
  int to_root(int i) const {
    return i - root(i);
  }
  
  int last_child(int i) const {
    int result = 0;
    for (int child = descendants_begin(i);
         child != descendants_end(i);
         child += descendants_size(child)) {
      result = child;
    }
    return result;
  }
  
  int n_leaves() const {
    int count = 0;
    for (auto& node : nodes_) {
      if (node.to_next == node.descendants_size) {
        count += 1;
      }
    }
    return count;
  }
  
  /* MODIFIERS */
  
  int push_back_root(const Idem idem) {
    int new_root = nodes_.size();
    root_idems_[new_root] = idem;
    nodes_.push_back(Node());
    return new_root;
  }
  
  /* Add subtree to a late node: make old subroot a child of the new subroot.
   * We do not update the new subroot (in particular, descendants_size).
   */
  int push_back_subtree(const int new_subroot,
                         const Weights new_weights,
                         const std::string new_label,
                         const int old_subroot,
                         const Node_container& old_nodes) {
    int new_child = nodes_.size();
    int subtree_size = old_nodes.descendants_size(old_subroot);
    nodes_.emplace_back(new_child - new_subroot, 1, subtree_size, new_weights
#ifdef DRAW
    , new_label
#endif  // DRAW
    );
    nodes_.insert(
      nodes_.end(),
      old_nodes.nodes_.begin() + old_subroot + 1,
      old_nodes.nodes_.begin() + old_subroot + subtree_size
    );
    nodes_[new_subroot].descendants_size += subtree_size;
    return new_child;
  }
  
  void erase_subtree_nodes(const int subroot) {
    if (subroot == 0) {
      root_idems_.erase(0);
    }
    else if (is_root(subroot)) {
      root_idems_.erase(subroot);
      increase_right_edge_(root(subroot), descendants_size(subroot));
    }
    else if (is_first_child(subroot)) {
      nodes_[parent(subroot)].to_next += descendants_size(subroot);
    }
    else {  // other child
      int child = descendants_begin(parent(subroot));
      int previous_child = 0;
      
      while (child != subroot) {
        previous_child = child;
        child += descendants_size(child);
      }
      
      increase_right_edge_(previous_child, descendants_size(subroot));
    }
  }
  
  void clear_nodes() {
    nodes_.clear();
    root_idems_.clear();
  }
  
 private:
  void increase_right_edge_(const int node, const int offset) {
    if (has_children(node)) {
      int child = last_child(node);
      nodes_[node].descendants_size += offset;
      increase_right_edge_(child, offset);
    }
    else {
      nodes_[node].descendants_size += offset;
      nodes_[node].to_next += offset;
    }
  }
  
 public:
  /* Pruning stuff */
  
  std::vector< int > node_offsets() const {
    std::vector< int > offsets(nodes_.size() + 1, -1);
    
    std::vector< int > something;
    
    int node = root_idems_.begin()->first;
    int offset = node;
    
    for (; node != nodes_.size(); node += to_next(node)) {
      offsets[node] = offset;
      offset += to_next(node) - 1;
    }
    
    offsets.back() = offset;
    return offsets;
  }
  
  void prune_nodes(const std::vector< int >& offsets) {
    std::vector< Node > new_nodes;
    Root_handle_container new_root_idems;
    
    for (int i = 0; i != nodes_.size(); ++i) {
      if (offsets[i] >= 0) {
        if (is_root(i)) {
          new_root_idems[new_nodes.size()] = root_idems_[i];
          new_nodes.emplace_back(0, 1,
            descendants_size(i) + offsets[i] - offsets[descendants_end(i)],
            weights(i)
#ifdef DRAW
            , label(i)
#endif  // DRAW
          );
        }
        else {
          new_nodes.emplace_back(
            to_parent(i) - offsets[i] + offsets[parent(i)],
            1,
            descendants_size(i) + offsets[i] - offsets[descendants_end(i)],
            weights(i)
#ifdef DRAW
            , label(i)
#endif  // DRAW
          );
        }
      }
    }
    
    nodes_.swap(new_nodes);
    root_idems_.swap(new_root_idems);
  }
  
  /* Poincaré polynomial stuff */
  
  template< class Polynomial >
  Polynomial poincare_polynomial(const Idem& idem = Idem("0")) const {
    // find root with given idempotent
    for (const auto& value_pair : root_idems_) {
      if (value_pair.second == idem) {
        return poincare_polynomial_at_< Polynomial >(value_pair.first);
      }
    }
    return Polynomial(0);
  }
  
  Weights generator_weights(int leaf) const {
    Weights weights = {0, 0};
    for (auto node_it = ascender(leaf); node_it.valid(); ++node_it) {
      Weights new_weights = this->weights(*node_it);
      weights.first += new_weights.first;
      weights.second += new_weights.second;
    }
    return weights;
  }
  
 private:
  template< class Polynomial >
  Polynomial poincare_polynomial_at_(const int node) const {
    if (!has_children(node)) {
      return Polynomial(1);
    }
    else {
      Polynomial poly = 0;
      for (int child = descendants_begin(node);
           child != descendants_end(node);
           child += descendants_size(child)) {
        Polynomial child_poly = poincare_polynomial_at_< Polynomial >(child);
        child_poly *= weights(node);
        poly += child_poly;
      }
      return poly;
    }
  }
  
 public:
  /* I/O interface */
  
  friend std::ostream& operator<<(std::ostream& os, const Node_container& nc) {
    for (auto& node : nc.nodes_) {
      os << node << " ";
    }
    return os;
  }
  
#ifdef DRAW
  public:
  /* TeXify */
  
  void TeXify(std::ofstream& write_file, const bool independent = true) const {
    Grid_ grid;
    
    if (independent) {
      write_file << "\\begin{tikzpicture}[suffix forest]" << std::endl;
    }
    
    for (auto value_pair : root_idems_) {
      add_to_grid_(grid, 0, value_pair.first);
    }
    
    /* Put idempotents for the roots */
    int i = 0;
    for (auto value_pair : root_idems_) {
      Grid_point_& grid_point = grid[0][i];
      write_file << "\\node ("
                   << grid_point.first
                   << ") at ("
                   << grid_point.second
                   << ","
                   << 0
                   << ") {" << value_pair.first << "\\nodeLabel{"
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
                   << ") {" << grid_point.first << "};\n"
                   << "\\draw[->] ("
                   << grid_point.first
                   << ") -- node[in place]{$"
                   << nodes_[grid_point.first].label
                   << "$} ("
                   << (grid_point.first - nodes_[grid_point.first].to_parent)
                   << ");" << std::endl;
      }
    }
    
    if (!grid.empty()) {
      float y_poly = static_cast<float>(grid.size() - 1) * y_sep_ + poly_sep_;
      for (Grid_point_& grid_point : grid.back()) {
        const Weights weights = generator_weights(grid_point.first);
        write_file << "\\node at ("
                   << grid_point.second
                   << ","
                   << y_poly
                   << ") {\\footnotesize$"
                   << weights.first
                   << ", "
                   << weights.second
                   << "$};" << std::endl;
      }
    }
    
    if (independent) {
      write_file << "\\end{tikzpicture}" << std::flush;
    }
  }
  
 private:
  using Grid_point_ = std::pair< int, float >;  // index of node, x coordinate
  using Grid_layer_ = std::vector< Grid_point_ >;  // left to right?
  using Grid_ = std::vector< Grid_layer_ >;  // lowest to highest
  
  // minimal horizontal distance between nodes
  static constexpr float min_x_sep_ = 1.;
  // vertical distance between nodes
  static constexpr float y_sep_ = .8;
  // vertical distance between polynomial and leaf
  static constexpr float poly_sep_ = .3;
  
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
      int first_child = grid[layer + 1].size();
      
      for (int child = descendants_begin(node);
           child != descendants_end(node);
           child += descendants_size(child)) {
        add_to_grid_(grid, layer + 1, child);
      }
      x = (grid[layer + 1][first_child].second + grid[layer + 1].back().second)
        / 2.;
      grid[layer].push_back(Grid_point_(node, x));
    }  // with children
  }
#endif  // DRAW
  
 protected:
  std::vector< Node > nodes_;
  
  Root_handle_container root_idems_;
};

#endif  // TEST_NODE_CONTAINER_H_