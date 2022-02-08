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

#ifndef PLANAR_DIAGRAM_H_
#define PLANAR_DIAGRAM_H_

#include <algorithm>  // min_element
#include <array>
#include <sstream>  // stringstream
#include <string>
#include <utility>  // pair
#include <vector>

/* Planar diagrams. Adapted from ComputeHFKv2/Diagrams.h and .cpp
 * 
 * A planar diagram is stored as a vector of integers.
 */
class Planar_diagram {
 public:
  Planar_diagram() { }
  
  Planar_diagram(std::vector< int > data) :
    data_(data)
  { }
  
  Planar_diagram(std::string pd_string) :
    data_()
  {
    int start = pd_string.find_first_of("0123456789", 0);
    while (start < std::string::npos) {
      int stop = pd_string.find_first_not_of("0123456789", start);
      int edge = 0;
      for (int i = start; i < stop; ++i) {
        edge = 10 * edge + (pd_string[i] - '0');
      }
      data_.push_back(edge);
      start = pd_string.find_first_of("0123456789", stop);
    }
    
    // Planar diagram number edges starting from 1
    if (!data_.empty()) {
      int min_edge = *std::min_element(data_.begin(), data_.end());
      for (int& edge : data_) {
        edge += 1 - min_edge;
      }
    }
  }
  
  std::vector< int > data() const {
    return data_;
  }
  
  std::string to_string() const {
    std::stringstream ss;
    ss << "PD[";
    for (int crossing = 0; crossing < data_.size() / 4; ++crossing) {
      ss << "X["
        << data_[4 * crossing] << ","
        << data_[4 * crossing + 1] << ","
        << data_[4 * crossing + 2] << ","
        << data_[4 * crossing + 3] << "],";
    }
    ss.seekp(-1, ss.cur);
    ss << "]" << std::flush;
    
    return ss.str();
  }
  
  /* Search for a small-girth Morse list. Adapted from ComputeHFKv2/Diagrams.cpp
   * 
   * We construct an upper knot diagram, adding crossings one by one. We choose
   * a crossing as the first crossing, and then we successively choose
   * crossings that are strongly connected to the upper knot diagram. In
   * certain cases, we will need to add local extrema.
   * 
   * The number of possible Morse lists may be exponentially large, so we limit
   * to a maximal number of tries and test choices randomly.
   * 
   * To evaluate the size of the Morse list, we also use a cost function. This
   * cost increases when we need to add local minima, and is greater the more
   * crossings we have added.
   */
  std::vector< int > get_legacy_morse_code(int max_attempts = 10000) const {
    int n_crossings = data_.size() / 4;
    std::vector< int > smallest_morse_list;
    int min_girth = 100;  // girth can't go over integer bit size, anyways
    long long min_cost = 1000000000;  // = 1e9
    max_attempts = std::min(100 + n_crossings * n_crossings, max_attempts);
    
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
      /* STEP 1: Initialize the Morse list with the first crossing */
      
      int girth = 4;  // girth of current Morse list
      long long cost = 0;
      
      // Choose first crossing and orientation, in order
      auto first_choice = div(attempt % data_.size(), 4);
      int first_crossing = first_choice.quot;
      int shift = first_choice.rem;
      
      // Edges at the bottom of the upper knot diagram
      std::vector< int > edges = {
        data_[4 * first_crossing + shift],
        data_[4 * first_crossing + (shift + 1) % 4],
        data_[4 * first_crossing + (shift + 2) % 4],
        data_[4 * first_crossing + (shift + 3) % 4]
      };
     
      /* Add two maxima and the first crossing to Morse list. We also calculate
       * the orientation of the maxima. The orientation of the knot follows the
       * edge numbering in increasing order.
       * There are twice as many edges as crossings.
       */
      std::vector< int > morse_list = {1001, 1, 1001, 3, -2};
      if ((edges[2] - edges[0]) % (2 * n_crossings) == 1) {
        morse_list[0] = 1000;
      }
      if ((edges[3] - edges[1]) % (2 * n_crossings) == 1) {
        morse_list[2] = 1000;
      }
      if (shift % 2 == 0) {
        morse_list[4] = 2;
      }
      
      std::vector< bool > added(n_crossings, false);  // added crossings
      added[first_crossing] = true;
      //std::cout << "\n[Da] first crossing " << first_crossing << ": " << edges << std::endl;
      
      /* STEP 2: Iteratively add a maximally connected crossing */
      for (int n_added = 1; n_added < n_crossings; ++n_added) {
        // Get maximally connected crossings
        auto result = get_max_connections_(data_, added, edges);
        auto max_connections = result.first;
        auto max_connected_crossings = result.second;
        
        if (max_connected_crossings.empty()) {
          //std::cout << "[Da] no good crossings found" << std::endl;
         goto attempt_end;  // give up on this attempt
        }
        
        // Nearly uniform random choice of next crossing
        int next_crossing = max_connected_crossings[rand() % max_connected_crossings.size()];
        added[next_crossing] = true;
        
        // Update morse_list, edges, and cost
        extend_Morse_list_(morse_list, edges, cost, data_, next_crossing, max_connections, n_added);
        
        if (edges.size() > girth) {
          girth = edges.size();
        }
        if (girth > min_girth) {
          //std::cout << "[Da] girth too big" << std::endl;
          goto attempt_end;  // give up on this attempt
        }
      }
      
      // Update smallest Morse list in the case of a successful attempt
      if (girth < min_girth or (girth == min_girth and cost < min_cost)) {
        min_girth = girth;
        min_cost = cost;
        smallest_morse_list = morse_list;
      }
      
      attempt_end: ;
    }
    
    return smallest_morse_list;
  }
    
 private:
  /* Find the crossings, not yet added to the Morse list, that share a maximal
   * number of edges with the vector edges.
   */
  inline std::pair< int, std::vector< int > > get_max_connections_(
    const std::vector< int >& pd,
    const std::vector< bool >& added,
    const std::vector< int >& edges
  ) const {
    std::vector< int > max_connected_crossings;
    int max_connections = 1;
  
    //std::cout << "[Da] current edges " << edges << std::endl;
    // Find maximally connected crossings
    for (int crossing = 0; crossing < (pd.size() / 4); ++crossing) {
    if (added[crossing]) {
      continue;
    }
  
    // static variables for memory optimization
    static std::vector< int > pos_in_edges;
    static std::vector< int > pos_in_crossing;
    pos_in_edges.clear();
    pos_in_crossing.clear();
  
    for (int i = 0; i < edges.size(); ++i) {
      int j = 0;
      while (j < 4 and pd[4 * crossing + j] != edges[i]) {
        ++j;
      }
    
      if (j != 4) {
        pos_in_edges.push_back(i);
        pos_in_crossing.push_back(j);
      }
    }
    int connections = pos_in_edges.size();
  
    if (connections == 0) continue;
    else if (pos_in_edges.back() - pos_in_edges.front() > connections) {
      continue;  // crossing attaches to temp in disjoint intervals
    }
  
    //std::std::cout << "pos in temp " << pos_in_edges << " pos in crossing " << pos_in_crossing << std::std::endl;
  
    for (int i = 0; i + 1 < connections; ++i) {
      if ((4 - pos_in_crossing[i + 1] + pos_in_crossing[i]) % 4 != 1) {
      // crossing does not attach counter-clockwise
      //std::std::cout << "[Da] data " << pos_in_crossing[i] << " " << pos_in_crossing[i + 1] << std::std::endl;
      //std::std::cout << "[Da] not counterclockwise " << crossing_edges << std::std::endl;
      goto crossing_loop_end;  // exit loop + jump using goto
      }
    }
  
    if (connections == max_connections) {
      max_connected_crossings.push_back(crossing);
    }
    else if (connections > max_connections) {
      max_connections = connections;
      max_connected_crossings.clear();
      max_connected_crossings.push_back(crossing);
    }
    crossing_loop_end: ;
    }
  
    return std::make_pair(max_connections, max_connected_crossings);
  }

  /* Extend the Morse list based on the new crossing and its connectivity.
   * Update the cost of the current Morse list.
   * 
   * Depending on the connectivity, we may need to add local extrema.
   */
  inline void extend_Morse_list_(
    std::vector< int >& morse_list,
    std::vector< int >& edges,
    long long& cost,
    const std::vector< int >& pd,
    const int next_crossing,
    const int connectivity,
    const int n_added
  ) const {
    const std::vector< int > crossing_edges(pd.begin() + 4 * next_crossing, pd.begin() + 4 * (next_crossing + 1));
    //std::cout << "[Da] next crossing " << next_crossing << ": " << crossing_edges << std::endl;
    int first_pos = 0;
    while (
      first_pos < edges.size()
      and edges[first_pos] != crossing_edges[0]
      and edges[first_pos] != crossing_edges[1]
      and edges[first_pos] != crossing_edges[2]
      and edges[first_pos] != crossing_edges[3]
    ) {
      ++first_pos;
    }
    int crossing_first_pos = 0;
    while (crossing_edges[crossing_first_pos] != edges[first_pos]) {
      ++crossing_first_pos;
    }
  
    // Compute cost of adding Morse event
    long long added_cost = n_added * edges.size() * edges.size();
    if (connectivity == 3) {
      cost += (first_pos + 1) * added_cost;
    }
    else if (connectivity == 4) {
      cost += (2 * first_pos + 1) * added_cost;
    }
  
    // Add Morse events according to the new crossing
    /* If connectivity == 1, add a maximum and a crossing as follows:
     * 
     *     | |   _   |
     *     |  \ / \  |
     *     |   X   | |
     *     |  | |  | |
     * 
     */
    if (connectivity == 1) {
    int left_pos = (crossing_first_pos + 1) % 4;
    int right_pos = (crossing_first_pos + 3) % 4;
    if ((crossing_edges[right_pos] - crossing_edges[left_pos]) % (pd.size() / 2) == 1) {
      morse_list.push_back(1000);  // maximum oriented left to right
    }
    else {
      morse_list.push_back(1001);  // maximum oriented right to left
    }
    morse_list.push_back(first_pos + 2);  // position of maximum
    if (crossing_first_pos % 2 == 0) {
      morse_list.push_back(first_pos + 1);  // positive crossing
    }
    else {
      morse_list.push_back(-first_pos - 1);  // negative crossing
    }
    }
  
    // If connectivity == 2, nothing complicated happens.
    else if (connectivity == 2 and crossing_first_pos % 2 == 0) {
      morse_list.push_back(first_pos + 1);
    }
    else if (connectivity == 2) {
      morse_list.push_back(-first_pos - 1);
    }
  
    /* If connectivity >= 3, then add a crossing and a minimum, and move
     * the minimum to the left:
     * 
     *     | |  | | | |
     *     | | /  _X  |
     *     |  X_ /  | |
     *     | /  X   | |
     *      X_ / |  | |
     *     /  X  |  | |
     *     \_/ | |  | |
     * 
     * If connectivity == 4, then we need to do this process twice.
     */
    else if (connectivity >= 3) {
      if (crossing_first_pos % 2 == 0) {
        morse_list.push_back(-first_pos - 2);  // negative crossing
      }
      else {
        morse_list.push_back(first_pos + 2);  // positive crossing
      }
      for (int i = first_pos; i > 0; --i) {
        morse_list.push_back(i);
        morse_list.push_back(i + 1);
      }
      morse_list.push_back(-1000);  // minimum
    
      // repeat if connectivity == 4
      if (connectivity == 4) {
        for (int i = first_pos; i > 0; --i) {
          morse_list.push_back(i);
          morse_list.push_back(i + 1);
        }
        morse_list.push_back(-1000);
      }
    }
  
  
    // Update edge list
    for (int i = 0; i < connectivity; ++i) {
      edges.erase(edges.begin() + first_pos);  // erase crossing_edges[(4 - i crossing_first_pos) % 4]
    }
    for (int i = connectivity; i < 4; ++i) {
      edges.insert(edges.begin() + first_pos, crossing_edges[(4 - i + crossing_first_pos) % 4]);
    }
  }
  
  std::vector< int > data_;
};

#endif  // PLANAR_DIAGRAM_H_