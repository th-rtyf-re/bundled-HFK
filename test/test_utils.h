#include <iostream>
#include <utility>
#include <vector>

template< class Container >
class Range {
 public:
  using Iterator = typename Container::iterator;
  
  Range(Container container) :
    begin_(container.begin()),
    end_(container.end())
  { }
  
  Iterator begin() {
    return begin_;
  }
  
  Iterator end() {
    return end_;
  }
 private:
  Iterator begin_;
  Iterator end_;
};

template< class Container >
class Composite_range;

template< class Container >
class Composite_iterator {
 public:
  using Iterator = typename Container::iterator;
  using Value = typename Container::value_type;
  
  Composite_iterator() { }
  
  Composite_iterator(Composite_range< Container >* comp_range, Iterator iterator, int current_range) :
    iterator_(iterator),
    comp_range_(comp_range),
    current_range_(current_range)
  { }
  
  Value operator*() const {
    return *iterator_;
  }
  
  Composite_iterator& operator++() {
    ++iterator_;
    if (iterator_ == comp_range_->ends_[current_range_]) {
      ++current_range_;
      if (current_range_ < comp_range_->ends_.size()) {
        iterator_ = comp_range_->begins_[current_range_];
      }
    }
    return *this;
  }
  
  bool operator!=(const Composite_iterator& other) {
    return current_range_ != other.current_range_ or iterator_ != other.iterator_;
  }
  
 private:
  Iterator iterator_;
  
  Composite_range< Container >* comp_range_;
  
  int current_range_;
};


template< class Container >
class Composite_range {
 public:
  using Iterator = typename Container::iterator;
  using Comp_iterator = Composite_iterator< Container >;
  
  friend Comp_iterator;
  
  Composite_range() { }
  
  Comp_iterator begin() noexcept {
    if (begins_.empty()) {
      return empty_composite_iterator_;
    }
    return Comp_iterator(this, begins_.front(), 0);
  }
  
  Comp_iterator end() noexcept {
    if (ends_.empty()) {
      return empty_composite_iterator_;
    }
    return Comp_iterator(this, ends_.back(), ends_.size());
  }
  
  void add_range(Iterator begin, Iterator end) {
    if (begin != end) {
      begins_.push_back(begin);
      ends_.push_back(end);
    }
  }
  
  void add_range(std::pair< Iterator, Iterator > it_pair) {
    if (it_pair.first != it_pair.second) {
      begins_.push_back(it_pair.first);
      ends_.push_back(it_pair.second);
    }
  }
  
 private:
  std::vector< Iterator > begins_;
  std::vector< Iterator > ends_;
  
  Comp_iterator empty_composite_iterator_;
};