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
template <typename V, typename E>
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
  using edge_value_type = E;
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

  /** Struct that stores information about an edge:
   *  the uids of its nodes, and the index of its value in the
   *  data structure that stores edge values. 
   */
  struct edge_info {
    edge_info(size_type n1, size_type n2, size_type val_index)
                   : n1_(n1), n2_(n2), val_index_(val_index) {
    }

    size_type n1() const {
      return n1_;
    }

    size_type n2() const {
      return n2_;
    }

    private:
      size_type n1_;
      size_type n2_;
    public:
      size_type val_index_;
  };

  /** Struct that stores information about an edge:
   *  the uids of its nodes, and its value.
   *  Used so that we do not store duplicate values for both
   *  directions of an edge.
   */
  struct edge_value_info {
    edge_value_info(size_type n1, size_type n2, edge_value_type val)
                   : n1_(n1), n2_(n2), val_(val) {
    }

    size_type n1() const {
      return n1_;
    }

    size_type n2() const {
      return n2_;
    }

    private:
      size_type n1_;
      size_type n2_;
    public:
      edge_value_type val_;
  };

  /** Struct that stores information about a node:
   *  its position, its value, its index, and its neighbors.
   */
  struct node_info {
    node_info() {
      // Construct an invalid node_info instance
    }
    node_info(Point pt, V val, size_type index)
            : pt_(pt), val_(val), index_(index) {
    }

    Point pt_;
    V val_;
    size_type index_;
    std::vector<edge_info> adj_;
  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    next_node_id_ = 0;
    num_edges_ = 0;
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
      return graph_->nodes_[uid_].pt_;
    }

    Point& position() {
      return graph_->nodes_[uid_].pt_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].index_;
    }

    /** Returns the value associated with this node. */
    node_value_type& value() {
      return graph_->nodes_[uid_].val_;
    }

    /** Returns the value associated with this node, as a const reference. */ 
    const node_value_type& value() const {
      return graph_->nodes_[uid_].val_;
    }

    /** Returns the number of other nodes connected to this node by an edge. */
    size_type degree() const {
      return graph_->nodes_[uid_].adj_.size();
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
      return n.index() == index() && n.graph_ == graph_;
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
      if (index() < n.index())
        return true;
      return false;
    }

   private:
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
      assert(graph_ != nullptr);
    }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
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
    Node n = Node(this, next_node_id_);
    nodes_.push_back(node_info(position, val, i2u_.size()));
    i2u_.push_back(next_node_id_++);
    return n;
  }

  /** Removes the given node from the graph.
   *  Also removes all of the node's incident edges from the graph.
   * @param[in] n   Node to be removed.
   * @return        The index of the node that was removed.
   * @pre           @a n does belong to this graph.
   * @post          @a n does not belong to this graph, i.e. !this.has_node(n).
   * @post          All edges @a e such that e.node1() == n or e.node2() == n
   *                no longer belong to this graph, i.e.
   *                !this.has_edge(e.node1(), e.node2()).
   * Some node indices may be changed. However, for all nodes @a n1 that belong
   * to the graph, the invariant n1.index() < num_nodes() is maintained.
   * Any outstanding node iterators @a ni with ni.index_ >= n.index()
   * are invalidated.
   * All edge iterators may be invalidated. All incident iterators of nodes
   * adjacent to @a n (as well as those of @a n itself) are invalidated.
   * Complexity: If we consider the maximum degree of the graph to be a
   *             fixed constant, the runtime is O(1).
   */
  size_type remove_node(const Node& n) {
    auto it = n.edge_begin();
    while (it != n.edge_end()) {
      remove_edge(*it);
      it = n.edge_begin();
    }
    size_type index = n.index(), uid = i2u_[index];
    if (index != size() - 1) {
      size_type u2 = i2u_.back();
      i2u_[index] = u2;
      nodes_[u2].index_ = index;
    }
    i2u_.pop_back();
    nodes_[uid] = {}; // clear the data, but don't modify the nodes_ vector
    return index;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this && n.uid_ < nodes_.size() &&
           0 <= nodes_[n.uid_].index_ && nodes_[n.uid_].index_ < i2u_.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u_[i]);
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
    
    /** Returns the value associated with this edge. */
    edge_value_type& value() {
      for (unsigned i = 0; i < graph_->nodes_[n1_].adj_.size(); i++) {
        if (graph_->nodes_[n1_].adj_[i].n2() == n2_)
          return graph_->edge_values_[
                         graph_->nodes_[n1_].adj_[i].val_index_].val_;
      }
      return graph_->edge_values_[0].val_; // silence compiler warnings, but
                                           // should never reach this point
    }

    /** Returns the value associated with this edge, as a const reference. */
    const edge_value_type& value() const {
      for (unsigned i = 0; i < graph_->nodes_[n1_].adj_.size(); i++) {
        if (graph_->nodes_[n1_].adj_[i].n2() == n2_)
          return graph_->edge_values_[
                         graph_->nodes_[n1_].adj_[i].val_index_].val_;
      }
      return graph_->edge_values[0].val_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, n1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, n2_);
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
    size_type n1_, n2_;
    Edge(const Graph* graph, size_type n1, size_type n2)
        : graph_(const_cast<Graph*>(graph)), n1_(n1), n2_(n2) {
      assert(graph_ != nullptr);
      assert(n1_ != n2_);
      assert(n1_ < graph_->nodes_.size());
      assert(n2_ < graph_->nodes_.size());
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return num_edges_;
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
   * Complexity: O(d), where d is the maximum degree of any node in the graph.
   */
  bool has_edge(const Node& a, const Node& b) const {
    auto v = nodes_[a.uid_].adj_;
    for (unsigned i = 0; i < v.size(); i++) {
      if (v[i].n2() == b.uid_)
        return true;
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
   * Complexity: O(1)
   */
  Edge add_edge(const Node& a, const Node& b,
                edge_value_type val = edge_value_type()) {
    if (!has_edge(a, b)) {
      edge_values_.push_back(edge_value_info(a.uid_, b.uid_, val));
      edge_info p1(a.uid_, b.uid_, edge_values_.size() - 1);
      edge_info p2(b.uid_, a.uid_, edge_values_.size() - 1);
      nodes_[a.uid_].adj_.push_back(p1);
      nodes_[b.uid_].adj_.push_back(p2);
      num_edges_++;
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Removes the edge with the given endpoints from the graph, if it exists.
   *  @param[in] n1   One endpoint of the edge to remove.
   *  @param[in] n2   The other endpoint of the edge to remove.
   *  @return The number of edges removed (because test_edges.cpp expects this)
   *  @post If there is an edge between @a n1 and @a n2 in the graph,
   *        it will be removed, i.e. !this.has_edge(n1, n2).
   *  If the edge was in the graph, all edge iterators are invalidated.
   *  The incident iterators of n1 and n2 are invalidated.
   *  Complexity: O(d), where d is the maximum degree of any node in the graph.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!has_node(n1) || !has_node(n2))
      return 0;
    bool found = false;
    std::vector<edge_info>& v1 = nodes_[n1.uid_].adj_;
    for (unsigned i = 0; i < v1.size(); i++) {
      if (v1[i].n2() == n2.uid_) {
        found = true;
        v1.erase(v1.begin() + i);
        break;
      }
    }
    if (!found)
      return 0;

    size_type val_index = 0;
    std::vector<edge_info>& v2 = nodes_[n2.uid_].adj_;
    for (unsigned i = 0; i < v2.size(); i++) {
      if (v2[i].n2() == n1.uid_) {
        val_index = v2[i].val_index_;
        v2.erase(v2.begin() + i);
        break;
      }
    }

    // Remove value of the deleted edge from edge_values_
    edge_values_[val_index] = edge_values_.back();
    edge_values_.pop_back();
    size_type aid = edge_values_[val_index].n1(),
              bid = edge_values_[val_index].n2();

    // Since we have moved the last entry in edge_values_,
    // update the adjacency list data to reflect this.
    std::vector<edge_info>& av = nodes_[aid].adj_;
    for (unsigned i = 0; i < av.size(); i++) {
      if (av[i].n2() == bid) {
        av[i].val_index_ = val_index;
        break;
      }
    }
    std::vector<edge_info>& bv = nodes_[bid].adj_;
    for (unsigned i = 0; i < bv.size(); i++) {
      if (bv[i].n2() == aid) {
        bv[i].val_index_ = val_index;
        break;
      }
    }

    num_edges_--;
    return 1;
  }
  
  /** Removes the given edge from the graph, if it is in the graph.
   *  @param[in] e    The edge to remove.
   *  @return         The number of edges removed.
   *  @post           @a e will be removed from the graph if it was in
   *                  the graph, i.e. !this.has_edge(e.node1(), e.node2()).
   *  If the edge was in the graph, all edge iterators are invalidated.
   *  The incident iterators of e.node1() and e.node2() are invalidated.
   *  Complexity: O(d), where d is the maximum degree of any node in the graph.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    i2u_.clear();
    next_node_id_ = 0;
    num_edges_ = 0;
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
    return NodeIterator(this, size());
  }

  /**
   * Removes the node currently pointed to by the given node iterator
   * (as well as all its incident edges) from the graph, if the iterator
   * is not at the end.
   * @param[in] n_it   Node iterator pointing the node to remove.
   * @return           A node iterator pointing to the next node (the node
   *                   that @a n_it would have pointed to after ++n_it), or to
   *                   the end, if @a n_it was one before the end.
   * @pre              @a n_it is a valid node iterator of this graph.
   * @post             If @a n_it != node_end(), the node pointed to by @a n_it
   *                   before the removal no longer belongs to this graph,
   *                   i.e. !this.has_node(*n_it).
   * @post             If @a n_it != node_end(), all edges @a e such that
   *                   e.node1() == *n_it or e.node2() == *n_it
   *                   no longer belong to this graph, i.e.
   *                   !this.has_edge(e.node1(), e.node2()).
   * Some node indices may be changed. However, for all nodes @a n that belong
   * to the graph, the invariant n.index() < num_nodes() is maintained.
   * Any outstanding node iterators @a ni with ni.index_ >= n_it.index_
   * are invalidated.
   * All edge iterators may be invalidated. All incident iterators of nodes
   * adjacent to @a *n_it (as well as those of @a *n_it itself) are invalidated.
   * Complexity: If we consider the maximum degree of the graph to be a
   *             fixed constant, the runtime is O(1).
   */
  node_iterator remove_node(node_iterator n_it) {
    if (n_it == node_end()) {
      return n_it;
    }
    size_type removed_index = remove_node(*n_it);
    return NodeIterator(this, removed_index);
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
      size_type n2 = graph_->nodes_[node_index_].adj_[index_].n2();
      return Edge(graph_, node_index_, n2);
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
      size_type node_uid = graph_->i2u_[node_index_];
      edge_info ei = graph_->nodes_[node_uid].adj_[index_];
      return Edge(graph_, node_uid, ei.n2());
    }

    /** Increment this EdgeIterator, so it points to the next edge.
     *  Returns a reference to this EdgeIterator.
     */ 
    EdgeIterator& operator++() {
      if (node_index_ != graph_->size()) {
        ++index_;
        fix();
      }
      return *this;
    }

    /** Determine whether this EdgeIterator is equal to the given one {ei}. */
    bool operator==(const EdgeIterator& ei) const {
      return graph_ == ei.graph_ && node_index_ == ei.node_index_ && index_ == ei.index_;
    }

   private:
    friend class Graph;
    const Graph* graph_;
    size_type node_index_;
    size_type index_;
    EdgeIterator(const Graph* graph, size_type node_index, int index)
                : graph_(graph), node_index_(node_index), index_(index) {
      fix();
    }

    /** Ensure that this iterator always points to a valid edge or to the end.
     *  Skips over edges such that node1() < node2(), so that each direction of
     *  an edge is only seen once.
     */
    void fix() {
      while (node_index_ < graph_->size()) {
        const std::vector<edge_info>& neighbors =
          graph_->nodes_[graph_->i2u_[node_index_]].adj_;
        while (index_ < neighbors.size()) {
          if (neighbors[index_].n2()
              > graph_->i2u_[node_index_]) {
            return;
          }
          index_++;
        }
        // if we reach this point, index must equal neighbors.size()
        index_ = 0;
        node_index_++;
      }
    }
  };

  /** Returns an iterator pointing to the first edge in this graph. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }
 
  /** Returns an iterator that points to the "past-the-end" edge
   *  in this graph.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, size(), 0);
  }

  /** Removes the edge currently pointed to by the given edge iterator,
   *  if the iterator is not at the end.
   *  @param[in] e_it   Edge iterator pointing the node to remove.
   *  @return           An edge iterator pointing to the next edge (the edge
   *                    that @a e_it would have pointed to after ++e_it), or to
   *                    the end, if @a e_it was one before the end.
   *  @pre              @a e_it is a valid edge iterator of this graph.
   *  @post             If @a e_it != edge_end(), *e_it will be removed
   *                    from the graph if it was in the graph,
                        i.e. !this.has_edge((*e_it).node1(), (*e_it).node2()).
   *  If @a e_it != edge_end(), all edge iterators are invalidated.
   *  The incident iterators of (*e_it).node1() and (*e_it).node2()
   *  are invalidated.
   *  Complexity: O(d), where d is the maximum degree of any node in the graph.
   */
   edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it == edge_end()) {
      return e_it;
    }
    remove_edge(*e_it);
    e_it.fix();
    return e_it;
  }

 private:
  size_type next_node_id_, num_edges_;
  std::vector<node_info> nodes_;
  std::vector<size_type> i2u_;
  std::vector<edge_value_info> edge_values_;
};

#endif // CME212_GRAPH_HPP
