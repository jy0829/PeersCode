#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
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
class Graph {
 private:
  // all private members are at the bottom

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;
  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

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
    next_node_id_ = 0;
    next_edge_id_ = 0;
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
      // do nothing
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Returns the value associated with this node. */
    node_value_type& value() {
      return graph_->node_values_[uid_];
    }

    /** Returns the value associated with this node, as a const reference. */ 
    const node_value_type& value() const {
      return graph_->node_values_[uid_];
    }

    /** Returns the number of other nodes connected to this node by an edge. */
    size_type degree() const {
      return graph_->adjacencyList_[uid_].size();
    }
    
    /** Returns an iterator pointing to the first edge adjacent to
     *  this node.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);  
    }

    /** Returns an iterator that points to the "past-the-end" edge
     *  adjacent to this node.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return n.uid_ == uid_ && n.graph_ == graph_;
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
      if (graph_ < n.graph_) // pointer comparison
        return true;
      if (graph_ > n.graph_)
        return false;
      if (uid_ < n.uid_)
        return true;
      return false;
    }

   private:
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
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
  Node add_node(const Point& position,
                const node_value_type& val = node_value_type()) {
    Node n = Node(this, next_node_id_++);
    nodes_.push_back(position);
    node_values_.push_back(val);
    adjacencySet_.push_back(std::unordered_map<size_type,size_type>());
    adjacencyList_.push_back(std::vector<size_type>());
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this && n.uid_ < nodes_.size();
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
      return graph_->node(n1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(n2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1() == e.node1() && node2() == e.node2()) ||
             (node1() == e.node2() && node2() == e.node1());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ < e.graph_)
        return true;
      if (graph_ > e.graph_)
        return false;
      Node n1 = node1(), n2 = node2(), temp;
      if (n2 < n1) { // assign n1 to be the "lesser" node, to ease comparison
        temp = n1;
        n1 = n2;
        n2 = temp;
      }
      Node n3 = e.node1(), n4 = e.node2();
      if (n4 < n3) { // similarly, assign n3 to be the "lesser" node
        temp = n3;
        n3 = n4; 
        n4 = temp;
      }
      return n1 < n3 || (n1 == n3 && n2 < n4);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type eid_;
    size_type n1_, n2_;
    Edge(const Graph* graph, size_type eid, size_type n1, size_type n2)
        : graph_(const_cast<Graph*>(graph)), eid_(eid), n1_(n1), n2_(n2) {
      assert(graph_ != nullptr);
      assert(n1_ != n2_);
      assert(n1_ < graph_->num_nodes());
      assert(n2_ < graph_->num_nodes());
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i, edges_[i].first, edges_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return adjacencySet_[a.uid_].count(b.uid_) ||
           adjacencySet_[b.uid_].count(a.uid_);
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
    if (adjacencySet_[a.uid_].count(b.uid_))
      return edge(adjacencySet_[a.uid_][b.uid_]);
    if (adjacencySet_[b.uid_].count(a.uid_)) { // edge (b,a) is in graph
      size_type i = adjacencySet_[b.uid_][a.uid_];
      // swap the order of the nodes in the vector storing edge data
      std::pair<size_type,size_type> p(a.uid_, b.uid_);
      edges_[i] = p;
      return edge(i);
    }
    // edge not found, so add it
    std::pair<size_type,size_type> p(a.uid_, b.uid_);
    edges_.push_back(p);
    adjacencySet_[a.uid_][b.uid_] = next_edge_id_;
    adjacencyList_[a.uid_].push_back(next_edge_id_);
    adjacencyList_[b.uid_].push_back(next_edge_id_);
    return edge(next_edge_id_++);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adjacencySet_.clear();
    adjacencyList_.clear();
    edges_.clear();
    next_node_id_ = 0;
    next_edge_id_ = 0;
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

    /** Return the Node that this NodeIterator is pointing to.
     *  Behavior is undefined if this NodeIterator is at or past the end. 
     */
    Node operator*() const {
      return graph_->node(index_);
    }
    
    /** Increment this NodeIterator, so it points to the next node.
     *  Returns a reference to this NodeIterator.
     */
    NodeIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Determine whether this NodeIterator is equal to the given one {ni}. */
    bool operator==(const NodeIterator& ni) const {
      return graph_ == ni.graph_ && index_ == ni.index_;
    }

   private:
    friend class Graph;
    const Graph* graph_;
    size_type index_;
    NodeIterator(const Graph* graph, size_type index)
                : graph_(graph), index_(index) {
    }
  };

  /** Returns an iterator pointing to the first node in this graph. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Returns an iterator that points to the "past-the-end" node
   *  in this graph.
   */
  node_iterator node_end() const {
    return NodeIterator(this, nodes_.size());
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
    
    /** Return the Edge that this IncidentIterator is pointing to.
     *  Behavior is undefined if this IncidentIterator is at or past the end. 
     */
    Edge operator*() const {
      size_type edge_index = graph_->adjacencyList_[node_index_][index_];
      Edge e = graph_->edge(edge_index);
      if (e.node1() == graph_->node(node_index_))
        return e;
      return graph_->reverseEdge(edge_index);
    }
    
    /** Increment this IncidentIterator, so it points to the next edge.
     *  Returns a reference to this IncidentIterator.
     */
    IncidentIterator& operator++() {
      ++index_;
      return *this;
    }
    
    /** Determine whether this IncidentIterator is equal to 
     *  the given one {iter}.
     */
    bool operator==(const IncidentIterator& iter) const {
      return graph_ == iter.graph_ && node_index_ == node_index_ &&
             index_ == iter.index_;
    }

   private:
    friend class Graph;
    const Graph* graph_;
    int node_index_;
    int index_;
    IncidentIterator(const Graph* graph, size_type node_index, size_type index)
                    : graph_(graph), node_index_(node_index), index_(index) {
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

    /** Return the Edge that this EdgeIterator is pointing to.
     *  Behavior is undefined if this EdgeIterator is at or past the end. 
     */
    Edge operator*() const {
      return graph_->edge(index_);
    }

    /** Increment this EdgeIterator, so it points to the next edge.
     *  Returns a reference to this EdgeIterator.
     */ 
    EdgeIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Determine whether this EdgeIterator is equal to the given one {ei}. */
    bool operator==(const EdgeIterator& ei) const {
      return graph_ == ei.graph_ && index_ == ei.index_;
    }

   private:
    friend class Graph;
    const Graph* graph_;
    size_type index_;
    EdgeIterator(const Graph* graph, size_type index)
                : graph_(graph), index_(index) {
    }
  };

  /** Returns an iterator pointing to the first edge in this graph. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
 
  /** Returns an iterator that points to the "past-the-end" edge
   *  in this graph.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.size());
  }

 private:
  size_type next_node_id_, next_edge_id_;
  std::vector<Point> nodes_;
  std::vector<V> node_values_;
  std::vector<std::pair<size_type,size_type>> edges_;
  std::vector<std::unordered_map<size_type,size_type>> adjacencySet_;
  std::vector<std::vector<size_type>> adjacencyList_;
  Edge reverseEdge(int i) const {
    return Edge(this, i, edges_[i].second, edges_[i].first);
  }
};

#endif // CME212_GRAPH_HPP
