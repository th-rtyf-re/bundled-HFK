#include <iostream>
#include <fstream>

#include "test_forest.h"

int main() {
  using Idem = typename Forest::Idem;
  using Alg_el = typename Forest::Arc::Coefficient_value;
  
  Forest f0;
  
  f0.declare_subtree(Idem("0"));
  f0.declare_subtree(Idem("1"));
  f0.lock_subtrees();
  
  std::cout << f0;
  
  Forest f1;
  f1.declare_subtree(Idem("0"), 1, 0, f0);
  f1.declare_subtree(Idem("0"), 2, 1, f0);
  f1.declare_subtree(Idem("1"), 3, 0, f0);
  f1.declare_subtree(Idem("1"), 4, 1, f0);
  f1.lock_subtrees(f0);
  
  f1.declare_arc(Alg_el(Idem("0"), Idem("1")), 1, 3, Idem("0"));
  f1.lock_arcs();
  
  std::cout << f1;
  
  Forest f2;
  f2.declare_subtree(Idem("0"), 5, 0, f1);
  f2.declare_subtree(Idem("0"), 6, 3, f1);
  f2.declare_subtree(Idem("1"), 7, 0, f1);
  f2.declare_subtree(Idem("1"), 8, 3, f1);
  f2.lock_subtrees(f1);
  
  Alg_el a1(Idem("0"), Idem("0"));
  Alg_el a2(Idem("1"), Idem("1"));
  
  f2.declare_arc(0, 7, a1);
  f2.declare_arc(3, 10, a1);
  f2.declare_arc(3, 2, a2);
  f2.declare_arc(11, 1, a2);
  
  f2.declare_arc(Alg_el(Idem("0"), Idem("1")), 5, 8, *f1.arcs().begin(), f1);
  f2.lock_arcs();
  
  std::cout << f2;
  
  for (const auto& arc : f2.arcs()) {
    for (const auto& arc2 : f2.arcs_to_source(arc)) {
      std::cout << "Concatenating " << arc2 << " and " << arc << " gives " << std::flush;
      std::cout << f2.concatenate(arc2, arc) << std::endl;
    }
  }
  
  f2.concatenate_zigzag(*(std::prev(std::prev(f2.arcs().end()))), *f2.arcs().begin(), *std::next(f2.arcs().begin()));
  
  std::ofstream suffix_forest("differential_suffix_forest.tex");
  f1.TeXify(suffix_forest);
  f2.TeXify(suffix_forest);
  suffix_forest.close();
  
  return 0;
}