template< class Forest_options >
class Node {
 public:
  using Gen_type = typename Forest_options::Gen_type;
  
  /* Constructor for a root node */
  Node() :
    d_parent_(0),
    d_after_children_(1),
    type_(-1)
  { }
  
  /* Constructor when we know the parent, and when there are no arcs yet.
   * Initializing next_sibling_ is a bit iffy.
   */
  Node(const int d_parent) :
    d_parent_(d_parent),
    d_after_children_(1),
    type_(-1)
  { }
  
  /* Access functions */
  int d_parent() const {
    return d_parent_;
  }
  
  int d_after_children() const {
    return d_after_children_;
  }
  
  Gen_type type() const {
    return type_;
  }
  
  /* Update functions */
  void add_d_after_children(int extra) {
    d_after_children_ += extra;
  }
  
  void set_d_parent(const int d_parent) {
    d_parent_ = d_parent;
  }
  
  void set_type(Gen_type type) {
    type_ = type;
  }
  
  /* Evaluation functions */
  bool has_parent() const {
    return d_parent_ > 0;
  }
  
  bool no_children() const {
    return d_after_children_ == 1;
  }
  
 private:
  int d_after_children_;  // relative position
  int d_parent_;  // relative position, positive integer
  
  Gen_type type_;
};

template< class Forest_options >
std::ostream& operator<<(std::ostream& os, const Node< Forest_options >& node) {
  os << "<" << node.d_parent() << "|" << static_cast<int>(node.type()) << "|" << node.d_after_children() << "> " << std::flush;
  return os;
}