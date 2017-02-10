#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph : private totally_ordered<Graph<V>> {
 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

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

  /** Predeclaration of Node value type from template. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

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
  Graph() {
    size_ = 0;
    next_id_ = 0;
    size_edges_ = 0;
    next_e_id_ = 0;
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
  class Node : private totally_ordered<Node> {
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
      return graph_->nodes_[id_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

    /** return a reference to the value of this (non-constant) node.  
        The value can be reassigned through this reference.
    */
    node_value_type& value() {
      return graph_->nodes_[id_].second;
    }

    /** return a reference to the value of this (constant) node.
        The value cannot be reassigned through this reference. 
    */
    const node_value_type& value() const {
      return graph_-> nodes_[id_].second;
    }

    /** return the number of edges adjacent to this node */
    size_type degree() const {
      return graph_->edge_adj_[id_].size();
    }

    /** return incident iterator iin such that if this node has adjacent edges,
     * *iin is one of these edges; if not, iin is the null pointer. Also, any
     * other edge adjacent to this node can be given by incrementing iin
     * some finite number of times and dereferencing the result */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, id_, 0);
    }

    /** return incident iterator iin_end such that iin_end is the null pointer,
     * and if we have an incident iterator iin of this node such that *iin = e
     * then incrementing iin some number of times will yield iin_end. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, id_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.index() != id_) {
        return false;
      }
      if (n.graph_ != graph_) {
        return false;
      }
      return true;
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
      if (graph_ < n.graph_) {
        return true;
      }
      if (n.graph_ < graph_) {
        return false;
      }
      if (id_ < n.index()) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type id_;
    /** Construct a valid node.
     */
    Node(const graph_type* graph, size_type id) {
      //assert(graph != nullptr);
      //assert(id < graph.size());
      graph_ = const_cast<graph_type*>(graph);
      id_ = id; 
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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
  Node add_node(const Point& position, const node_value_type& value = 
   node_value_type()) {
    Node n(this, next_id_);
    ++next_id_;
    ++size_;
    std::pair<Point, node_value_type> p;
    p = std::make_pair(position, value);
    nodes_.push_back(p);
    std::vector<std::pair<size_type, size_type>> v;
    edge_adj_.push_back(v);
    return n;     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this) {
      if (size_ > n.id_) {
        return true;
      }
    }
    return false;
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, end_points_.first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, end_points_.second);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ != e.graph_) {
        return false;
      }
      if ((end_points_.first == e.end_points_.first) and 
        (end_points_.second == e.end_points_.second)) {
        return true;
      }
      if ((end_points_.first == e.end_points_.second) and 
        (end_points_.second == e.end_points_.first)) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ < e.graph_) {
        return true;
      }
      if (e.graph_ < graph_) {
        return false;
      }
      if (id_ < e.id_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type id_;
    std::pair<size_type, size_type> end_points_;
    /** Create valid edge
     */
    Edge(const graph_type* graph, size_type id, Node a, Node b) {
      //assert(graph != nullptr);
      //assert(id < graph->edges_.size());
      //assert(a != b);
      id_ = id;
      graph_ = const_cast<graph_type*>(graph);
      end_points_.first = a.id_;
      end_points_.second = b.id_;
    }
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return size_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    Node n1(this, edges_[i].first);
    Node n2(this, edges_[i].second);
    Edge e(this, i, n1, n2);
    return e;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned i = 0; i < edge_adj_[a.id_].size(); i++) {
      if (edge_adj_[a.id_][i].first == b.id_) {
        return true;
      }
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
    if (has_edge(a, b)) {
      size_type index = find_edge_index(a, b); 
      return Edge(this, index, a, b);
    }
    else {
      Edge e(this, next_e_id_, a, b);
      std::pair<size_type, size_type> p1;
      std::pair<size_type, size_type> p2;
      std::pair<size_type, size_type> p3;
      p1 = std::make_pair(b.id_, next_e_id_);
      p2 = std::make_pair(a.id_, next_e_id_);
      p3 = std::make_pair(a.id_, b.id_);
      edge_adj_[a.id_].push_back(p1);
      edge_adj_[b.id_].push_back(p2);
      edges_.push_back(p3);
      ++next_e_id_;
      ++size_edges_; 
      return e;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    size_ = 0;
    next_id_ = 0;
    next_e_id_ = 0;
    size_edges_ = 0;
    nodes_.clear();
    edges_.clear(); 
    edge_adj_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Dereference NodeIterator. 
     * @pre id_ < graph_.size()
     * @return node
     * @post node is a node of graph_. */
    Node operator*() const {
      return Node(graph_, id_);
    }

    /** Increment node Iterator. */
    NodeIterator& operator++() {
      ++id_;
      return *this;
    }

    /** if the private elements id_ and graph_ are
     the same, return true. Otherwise return false*/
    bool operator==(const NodeIterator& ni) const {
      return ((ni.graph_ == graph_) && (ni.id_ == id_));
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type id_;
    /** Construct a valid Node Iterator. */
    NodeIterator(const graph_type* graph, size_type id) {
      //assert(graph_ != nullptr);
      id_ = id;
      graph_ = graph;
    }
  };

    /** returns node iterator ni such that if graph_.size() > 0 then *ni is a
     * node; otherwise, ni is the nullpointer. Also, for every node in the graph     * incrementing ni some finite number of times and dereferencing yields 
     * this node. */
    node_iterator node_begin() const {
      return NodeIterator(this, 0);
    }

    /** returns node iterator ni_end such that ni_end is the null pointer, and 
      * if ni is a node iterator such that *ni = n, for some node n, then 
      * incrementing ni some finite number of times yields ni_end. */
    node_iterator node_end() const {
      return NodeIterator(this, size_);
    }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    /** @pre Node n1(graph_, node_id_).degree() < neigh_id_. 
     * @return Edge e. 
     * @post e is an edge incident on n1.  */
    Edge operator*() const {
      Node n1(graph_, node_id_);
      Node n2(graph_, graph_->edge_adj_[node_id_][neigh_id_].first);
      return Edge(graph_, graph_->edge_adj_[node_id_][neigh_id_].second, n1, 
        n2); 
    }

    /** Increment incident iterator */
    IncidentIterator& operator++() {
      ++neigh_id_;
      return *this;
    }
 
    /** return true if graph and node id are the same, false otherwise*/
    bool operator==(const IncidentIterator& iit) const {
      return ((iit.graph_ == graph_) && (iit.node_id_ == node_id_) && 
        (iit.neigh_id_ == neigh_id_)); 
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type node_id_;
    size_type neigh_id_;
    /** Construct valid Incident Iterator */
    IncidentIterator(const graph_type* graph, size_type node_id,
      size_type neigh_id) {
      graph_ = graph;
      node_id_ = node_id;
      neigh_id_ = neigh_id;
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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

    /** @pre @a e_id_ < @a graph_->edges_.size(). 
     * @return Edge. 
     * @post return edge in the graph*/
    Edge operator*() const {
      Node n1(graph_, graph_->edges_[e_id_].first);
      Node n2(graph_, graph_->edges_[e_id_].second);
      return Edge(graph_, e_id_, n1, n2);
    }

    /** Increment edge iterator */
    EdgeIterator& operator++() {
      ++e_id_;
      return *this;
    }

    /** return true if graph and id are equal, false otherwise*/
    bool operator==(const EdgeIterator& ei) const {
      return ((ei.graph_ == graph_) && (ei.e_id_ == e_id_));
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type e_id_;
    /** Construct valid Edge iterator*/
    EdgeIterator(const graph_type* graph, size_type e_id) {
      graph_ = graph;
      e_id_ = e_id;
    }
  };

  /** @return edge_iterator ei.
   * @post *ei == e for some edge e such that for each edge e' in the edges of
   * the graph there exists a nonnegative integer n such that incrementing ei 
   * n times will yield *ei == e'
   * */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** @return edge_iterator ei.
   * @post *ei is not a reference to an edge, but if ei' is the edge iterator
   * that has been incremented the most while still remaining an edge when 
   * dereferenced, then ++ei' == ei*/
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.size());
  }

 private:
  size_type size_;
  size_type next_id_;
  size_type next_e_id_;
  size_type size_edges_;
  std::vector<std::pair<Point, node_value_type>> nodes_;
  std::vector<std::vector<std::pair<size_type, size_type>>> edge_adj_;
  std::vector<std::pair<size_type, size_type>> edges_;
  /** Find the index of an edge in a graph.
  * @pre @a a and @a b are distinct valid nodes of this graph for which 
  * has_edge(@a a, @a b) == true 
  * @return @a i, where edges_[i] is a pair containing @a a.id_ and @a b.id_ 
  */
   size_type find_edge_index(const Node& a, const Node& b) {
     for (unsigned i = 0; i < edge_adj_[a.id_].size(); i++) {
       if (edge_adj_[a.id_][i].first == b.id_) {
         return edge_adj_[a.id_][i].second;
       }
     }
     return size_edges_ + 1;
   }
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
