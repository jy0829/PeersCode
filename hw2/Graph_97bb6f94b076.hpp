#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP


/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <cmath>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E> class Graph {
 private:

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;

  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** HW1: Synonym for node value type. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
 
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Synonym for edge value type. */
  using edge_value_type = E;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
 
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
 
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
 
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes_(), adjacency_list() {

  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered <Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return ((*graph_).nodes_[node_id].first);
    }

    /** nonconst version of position. */
    Point& position() {
      return ((*graph_).nodes_[node_id].first);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_id;
    }

    incident_iterator edge_begin() const {
        return IncidentIterator(graph_, node_id, 0);
    }

    incident_iterator edge_end() const {
        return IncidentIterator(graph_, node_id, degree());
    }

    node_value_type& value(){
        return (*graph_).nodes_[node_id].second;
    }

    const node_value_type& value() const {
        return (*graph_).nodes_[node_id].second;
    }

    size_type degree() const {
      return graph_->adjacency_list[node_id].size();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (node_id == n.node_id and graph_ == n.graph_) {
          return true;
      } else {
          return false;
      }
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if (graph_ == n.graph_ and node_id < n.node_id) {
          return true;
      } else if (graph_ < n.graph_) {
          return true;
      } else {
          return false;
      }
    }

   private:
    friend class Graph;

    Graph* graph_;
    size_type node_id;

    //Construct a valid node:
    Node(const Graph* graph, size_type node_id) : graph_(const_cast<Graph*>(graph)),
        node_id(node_id) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& myval = node_value_type()) {
    nodes_.push_back(std::make_pair(position, myval));
    adjacency_list.push_back({});
    return Node(this, size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this and n.node_id < size()) {
      return true;
    } else {
      return false;
    }

  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered <Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, index_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //given how I am storing the second node I need to use the fetch func
      return Node(graph_, fetch_node2_id());
    }
 
    /** Return the value of the edge. */
    edge_value_type& value() {
      return graph_->adjacency_list[index_1][aux_index_2].second;
    }

    /** Const version of value. */
    const edge_value_type& value() const {
      return graph_->adjacency_list[index_1][aux_index_2].second;
    }
    
    /** Return the Euclidean length of the straight edge. */ 
    double length() const {
      return sqrt(pow(node1().position().x - node2().position().x, 2) +
                  pow(node1().position().y - node2().position().y, 2) +
                  pow(node1().position().z - node2().position().z, 2)); 
   
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ == e.graph_ and
          ((index_1 == e.index_1 and fetch_node2_id() == e.fetch_node2_id()) or
          (fetch_node2_id() == e.index_1 and index_1 == e.fetch_node2_id()))) {
         return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type auxmax = std::max(index_1, fetch_node2_id());
      size_type auxmin = std::min(index_1, fetch_node2_id());
      size_type newmax = std::max(e.index_1, e.fetch_node2_id());
      size_type newmin = std::min(e.index_1, e.fetch_node2_id());
      if ((graph_ < e.graph_) or ((graph_ == e.graph_) and 
         ((auxmin < newmin) or ((auxmin == newmin) and (auxmax < newmax))))) {
  	  return true;
  	} else {
          return false;
        }
    }

   private:
    friend class Graph;

    Graph* graph_;
    size_type index_1;
    size_type aux_index_2;

    //fetch node 2 id
    size_type fetch_node2_id() const{
        return graph_->adjacency_list[index_1][aux_index_2].first;
    }

    //valid edge constructor
    Edge (const Graph* graph, size_type index_1, size_type aux_index_2) :
      graph_(const_cast<Graph*>(graph)), index_1(index_1), aux_index_2(aux_index_2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return std::distance(edge_begin(),edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(), i);
  }
 

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type a_id = a.node_id;
    size_type b_id = b.node_id;

    for (size_type i = 0; i < adjacency_list[a_id].size(); ++i) {
      if (adjacency_list[a_id][i].first == b_id) { 
        return true; }
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    size_type a_id = a.node_id;
    size_type b_id = b.node_id;

    for(size_type i = 0; i < adjacency_list[a_id].size(); ++i){
      if(adjacency_list[a_id][i].first == b_id) { 
        return Edge(this, a_id, i); 
      }
    }
    adjacency_list[a_id].push_back(std::make_pair(b_id, edge_value_type()));
    adjacency_list[b_id].push_back(std::make_pair(a_id, edge_value_type()));
    return Edge(this, a_id, adjacency_list[a_id].size() - 1);
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adjacency_list.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable <NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** NodeIterator dereference operator */
    Node operator*() const  {
        return Node(graph_,current_);
    }

    /** NodeIterator increment operator */
    NodeIterator& operator++() {
      ++current_; 
      return *this;
    }
    
    /** NodeIterator equality operator */
    bool operator==(const NodeIterator& otheriter) const {
        return (graph_ == otheriter.graph_ and current_ == otheriter.current_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type current_;

    /** NodeIterator private valid constructor*/
    NodeIterator(const Graph* graph, size_type current)
     : graph_(const_cast<Graph*>(graph)), current_(current){}
  };

  /** instantiate a node_iterator */
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
      //return nodeiter;
  }

  /** instantiate a node_iterator */
  node_iterator node_end() const {
      return NodeIterator(this, size());
      //return nodeiter;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private equality_comparable <IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** IncidentIterator dereference operator */
    Edge operator*() const {
      return Edge(graph_, main_node_, other_aux_);
    }

    /** IncidentIterator increment operator */
    IncidentIterator& operator++() {
      ++other_aux_;
      return *this;
    }

    /** IncidentIterator equality operator */
    bool operator==(const IncidentIterator& incidentiter) const {
      return (graph_ == incidentiter.graph_ 
              and main_node_ == incidentiter.main_node_
              and other_aux_ == incidentiter.other_aux_);
    } 

   private:
    friend class Graph;

    Graph* graph_;
    size_type main_node_;
    size_type other_aux_;

    /** IncidentIterator private valid constructor */
    IncidentIterator(const Graph* graph, size_type main_node, size_type other_aux) :
    graph_(const_cast<Graph*>(graph)), main_node_(main_node), other_aux_(other_aux){};

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private equality_comparable <EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** EdgeIterator dereference operator */
    Edge operator*() const {
      return Edge(graph, main_node_, other_aux_);
    }

    /** EdgeIterator increment operator */
    EdgeIterator& operator++() {
      ++other_aux_;
      fix();
      return *this;
    }

    /** EdgeIterator equality operator */
    bool operator==(const EdgeIterator& edgeiter) const {
      return(graph == edgeiter.graph and main_node_ == edgeiter.main_node_
             and other_aux_ == edgeiter.other_aux_);
    }

   private:
    friend class Graph;
    
    Graph* graph;
    size_type main_node_;
    size_type other_aux_;

    /** this function is the nicest part of the whole homework. It
     ** finds the next valid edge and allowed me to get rid of the edge vector. */
    void fix() {
      while (true) {
        if (main_node_ == graph->num_nodes()) {
            break;
        } else if (other_aux_ == graph->adjacency_list[main_node_].size()){
          other_aux_ = 0;
          main_node_++;
        } else if ( main_node_ > graph->adjacency_list[main_node_][other_aux_].first){
          other_aux_++;
        } else {
          break;
        }
      }
    }

    /** EdgeIterator private valid constructor*/
    EdgeIterator(const Graph* g, size_type main_node, size_type other_aux) :
                 graph(const_cast<Graph*>(g)), main_node_(main_node), other_aux_(other_aux) {
                 fix();
    }

  };

  /** instantiate an EdgeIterator */
  edge_iterator edge_begin() const {
      return EdgeIterator(this,0,0);
  }
 
  /** instantiate an the first invalid EdgeItertor, the args required some thinking */
  edge_iterator edge_end() const {
      return EdgeIterator(this, num_nodes(), 0);
  }

  
  /** Remove the node @a n from the graph
   * @return the number of nodes deleted by this function 
   * @post has_node(@a n) == false
   * @post If old !has_node(@a n),   new num_nodes() == old num_nodes().
   *       Else,                     new num_nodes() == old num_nodes() - 1.
   *
   * @post If has_edge(edge.node1(), edge.node2()) 
   *       && (edge.node1() == n || edge.node2() == n), edge is removed from graph. 
   *
   * Can invalidate node indices -- old node(@a i) might not equal new node(@a i). 
   * Can invalidate node iterators
   * Can invalidate edge iterators 
   * Can invalidate incident iterators
   * Complexity: O(num_nodes())
   */
  size_type remove_node(const Node& n) {
  
    if (has_node(n)) { 
      for (size_type n_id  = 0; n_id < num_nodes(); ++n_id) {
        if (n_id != n.index()) {
          for (auto ii = adjacency_list[n_id].begin(); ii != adjacency_list[n_id].end(); ){ 
            if ( (*ii).first == n.index()) {
              ii = adjacency_list[n_id].erase(ii);
            } else if ((*ii).first > n.index()) {
              --((*ii).first);
              ++ii;
            } else {
              ++ii;
            }
          }  
        }
      }
      adjacency_list.erase(adjacency_list.begin() + n.index());
      nodes_.erase(nodes_.begin() + n.index());
      return 1;
    } else {
      return 0;
    }
  }

  /** Remove the node pointed by @a n_it from the graph
   * @return the number of nodes deleted by this function 
   * @post has_node(*(@a n_it)) == false
   * @post If old !has_node(@a n),   new num_nodes() == old num_nodes().
   *       Else,                     new num_nodes() == old num_nodes() - 1.
   *
   * @post If has_edge(edge.node1(), edge.node2()) 
   *       && (edge.node1() == n || edge.node2() == n), edge is removed from graph. 
   *
   * Can invalidate node indices -- old node(@a i) might not equal new node(@a i). 
   * Can invalidate node iterators
   * Can invalidate edge iterators 
   * Can invalidate incident iterators
   * Complexity: O(num_nodes()) 
   */
  node_iterator remove_node(node_iterator n_it) {

    node_type n = *n_it;
    for (size_type n_id  = 0; n_id < num_nodes(); ++n_id) {
      if (n_id != n.index()) {
        for (auto ii = adjacency_list[n_id].begin(); ii != adjacency_list[n_id].end(); ){ 
          if ( (*ii).first == n.index()) {
            ii = adjacency_list[n_id].erase(ii);
          } else if ((*ii).first > n.index()) {
            --((*ii).first);
            ++ii;
          }
        }
      }
    }
    adjacency_list.erase(adjacency_list.begin() + n.index());
    auto next = nodes_.erase(nodes_.begin() + n.index());   
    return (node_begin() + std::distance(nodes_.begin, next));
  } 
  

  /** Remove edge with node @a a and @a b from the graph.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return the number of edges deleted from the graph.
   * @post new has_edge(@a a, @a b) == false
   * @post If old !has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                         new num_edges() == old num_edges() - 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   * Can invalidate edge iterators. 
   * Can invalidate incident iterators.i
   *
   * Complexity: O(num_edges())
   */
  size_type remove_edge(const Node& a, const Node& b) {
   
    if (has_edge(a, b)) {
      for (size_type i = 0; i < adjacency_list[a.index()].size(); i++) {
        if (adjacency_list[a.index()][i].first == b.index()) {
          adjacency_list[a.index()].erase(adjacency_list[a.index()].begin() + i);
          break;
        }
      } 
      for (size_type i = 0; i < adjacency_list[b.index()].size(); i++) {
        if (adjacency_list[b.index()][i].first == a.index()) {
          adjacency_list[b.index()].erase(adjacency_list[b.index()].begin() + i);
          break;
        } 
      }
      return 1;
    } else {
      return 0;
    }   
  }

  /** Remove edge @a e from the graph.
   * @pre  e is a valid edge.
   * @return the number of edges deleted from the graph.
   * @post a = e.node1() and b = e.node2(), then new has_edge(@a a, @a b) == false
   * @post If old !has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                         new num_edges() == old num_edges() - 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   * Can invalidate edge iterators. 
   * Can invalidate incident iterators.i
   *
   * Complexity: O(num_edges())
   */
  size_type remove_edge(const Edge& e) {

    node_type a = e.node1();
    node_type b = e.node2();
    return remove_edge(a, b);
  }

  /** Remove edge pointed by @a e_it from the graph.
   * @pre  *(@a e_it) is a valid edge.
   * @return the edge iterator pointing to the next edge before deletion.
   * @post a = (*e_it).node1() and b = (*e_it).node2(), then new has_edge(@a a, @a b) == false
   * @post If old !has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                         new num_edges() == old num_edges() - 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i).
   * Can invalidate edge iterators. 
   * Can invalidate incident iterators.i
   *
   * Complexity: O(num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    edge_type next_edge = *(e_it + 1);
    remove_edge(*e_it);
    for (auto it = edge_begin(); it != edge_end(); ++it) {
      if (next_edge == *it) {
        return it;
      }
    }  
  }

 private:
  
  std::vector<std::pair<Point, node_value_type>> nodes_;
  std::vector<std::vector<std::pair<size_type, edge_value_type>>> adjacency_list;
 
};

#endif // CME212_GRAPH_HPP
