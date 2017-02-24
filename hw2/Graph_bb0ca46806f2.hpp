#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <tuple>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

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
  using graph_type = Graph<V, E>;

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

  /** Type of node values */
  using node_value_type = V;

  /** Type of edge values */
  using edge_value_type = E;

  /** Type of indexes and sizes. */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : points(), i2u(), no_nodes(0), no_edges(0) {}

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
    Node() {}

    /** Return this node's position. */
    const Point &position() const { return get(); }

    /** Return this node's position. */
    Point &position() { return this->graph_->points[uid_].p; }

    /** Check if the node it valid. */
    size_type is_valid() const {
      return uid_ >= 0 && uid_ < graph_->points.size() &&
             graph_->points[uid_].idx < graph_->i2u.size() &&
             graph_->i2u[graph_->points[uid_].idx] == uid_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return graph_->points[uid_].idx; }

    /** Returns the reference to @a val of node, a value of
     *  type node_value_type&
     * @return a reference to @a val of this node
     *
     * @pre the node calling should be valid
     */
    node_value_type &value() { return this->graph_->points[uid_].val; }

    /** Returns the reference to @a val of node, a value of
     *  type node_value_type&
     * @return a const reference to @a val of this node
     *
     * @pre the node calling should be valid
     */
    const node_value_type &value() const {
      return (const node_value_type &)this->graph_->points[uid_].val;
    }

    /** Returns the number of neighbors of node
     * @return the number of edges incident to current node.
     *
     * @pre the node calling should be valid
     */
    size_type degree() const { return graph_->points[uid_].adj.size(); }

    /** Returns an incident iterator pointing to the
     * "first" neigbor of this node
     * @return iterator of type incident_iterator
     *
     * @pre the node calling should be valid
     */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_->points[this->uid_].idx, this->graph_, 0);
    }

    /** Returns an incident iterator pointing to one past the
     * "last" neigbor of @a this node
     * @return iterator of type incident_iterator
     *
     * @pre the node calling should be valid
     */
    incident_iterator edge_end() const {
      return incident_iterator(graph_->points[this->uid_].idx, this->graph_,
                               this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {

      // check if the node belongs to the same graph and has the same uid
      if (this->graph_ == n.graph_ && this->uid_ == n.uid_)
        return true;

      return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const {

      // check if the node belongs to the same graph and if uid is less
      if (this->graph_ == n.graph_ && this->uid_ < n.uid_)
        return true;

      // compare graph pointers next
      if (this->graph_ < n.graph_)
        return true;

      return false;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer to original graph
    graph_type *graph_;

    // The index of the node
    size_type uid_;

    // Private constructor
    Node(const graph_type *graph, size_type uid)
        : graph_(const_cast<graph_type *>(graph)), uid_(uid) {}

    // Define get function
    Point &get() const { return graph_->points[uid_].p; }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { return no_nodes; }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position,
                const node_value_type &nvt = node_value_type()) {

    // initialize new node
    node_wrapper new_node;
    new_node.p = position;
    new_node.val = nvt;
    new_node.adj = std::vector<edge_wrapper>();
    new_node.idx = no_nodes;

    // add element to vector
    points.push_back(new_node);
    i2u.push_back(points.size() - 1);

    // increment number of nodes
    no_nodes += 1;

    // return node
    return Node(this, points.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {

    // check if node is in the graph
    if (n.uid_ < points.size() && this == n.graph_ && n.is_valid()) {
      return true;
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

    if (i < no_nodes)
      return Node(this, i2u[i]);

    // if i is bigger than the size then return invalid node
    return Node();
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
    Edge() {}

    /** Computes the length of initial edge
     *
     */
    double length() const {
      Point n1 = node1().position();
      Point n2 = node2().position();
      return norm(n1 - n2);
    }

    /** Return a node of this Edge */
    Node node1() const { return graph_->node(nid1_); }

    /** Return the other node of this Edge */
    Node node2() const { return graph_->node(nid2_); }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge 0.01;
     */
    bool operator==(const Edge &e) const {

      // check if two edges belong to the same graph and have the same euid
      if (this->graph_ == e.graph_ &&
          std::min(nid1_, nid2_) == std::min(e.nid1_, e.nid2_) &&
          std::max(nid1_, nid2_) == std::max(e.nid1_, e.nid2_)) {
        return true;
      }

      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {

      // a global comparison between two edges
      if (std::tie(graph_, std::min(nid1_, nid2_), std::max(nid1_, nid2_)) <
          std::tie(e.graph_, std::min(e.nid1_, e.nid2_),
                   std::max(e.nid1_, e.nid2_)))
        return true;

      return false;
    }

    /** Returns the value of the current edge */
    edge_value_type &value() {
      edge_value_type to_return;
      for (size_type i = 0; i < graph_->points[graph_->i2u[nid1_]].adj.size();
           i++) {
        if (graph_->points[graph_->i2u[nid1_]].adj[i].id ==
            graph_->i2u[nid2_]) {
          return graph_->points[graph_->i2u[nid1_]].adj[i].eval;
        }
      }

      // Default return
      return graph_->points[0].adj[0].eval;
    }

    /** Returns the value of the current edge */
    const edge_value_type &value() const {
      edge_value_type to_return;
      for (size_type i = 0; i < graph_->points[graph_->i2u[nid1_]].adj.size();
           i++) {
        if (graph_->points[graph_->i2u[nid1_]].adj[i].id ==
            graph_->i2u[nid2_]) {
          return graph_->points[graph_->i2u[nid1_]].adj[i].eval;
        }
      }

      // Default return
      return graph_->points[0].adj[0].eval;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // pointer to the graph
    graph_type *graph_;

    // node id 1
    size_type nid1_;

    // node id 2;
    size_type nid2_;

    // Private constructor
    Edge(const graph_type *graph, size_type nid1, size_type nid2)
        : graph_(const_cast<graph_type *>(graph)),
          nid1_(graph_->points[nid1].idx), nid2_(graph_->points[nid2].idx) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return no_edges; }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(), i); // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {

    for (size_type i = 0; i < this->points[a.uid_].adj.size(); i++) {
      if (this->points[a.uid_].adj[i].id == b.uid_) {
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
  Edge add_edge(const Node &a, const Node &b) {
    // check if edge exists
    if (has_edge(a, b)) {
      return Edge(this, a.uid_, b.uid_);
    }

    // add new edge
    edge_wrapper new_edge1;
    new_edge1.id = b.uid_;
    new_edge1.eval = edge_value_type();
    edge_wrapper new_edge2;
    new_edge2.id = a.uid_;
    new_edge2.eval = edge_value_type();
    points[a.uid_].adj.push_back(new_edge1);
    points[b.uid_].adj.push_back(new_edge2);

    // increment number of edges
    no_edges++;

    // return edge constructed
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points.clear();
    i2u.clear();
    no_edges = 0;
    no_nodes = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Dereference the iterator to return the node it points to
     * @return an object of type node_type in the graph object
     *
     * @pre iterator != graph.node_end() && iterator != NodeIterator()
     */
    Node operator*() const { return this->gp_->node(this->nuid_); }

    /** Increment the iterator so that it points to the next node
     * @return reference to incremented node iterator
     *
     * @pre iterator != graph.node_end() && iterator != NodeIterator()
     */
    NodeIterator &operator++() {
      this->nuid_++;
      return *this;
    }

    /** Compare equality for two node iterators
     * @param[in] @a ni is the node_iterator to be compared
     * @return a boolean that is true if the two node_iterators are true
     */
    bool operator==(const NodeIterator &ni) const {
      return (this->gp_ == ni.gp_ && this->nuid_ == ni.nuid_);
    }

    /** Compare not equality for two node iterators
     * @param[in] @a ni is the node_iterator to be compared
     * @return !(*this == ni)
     */
    bool operator!=(const NodeIterator &ni) const { return !(*this == ni); }

  private:
    friend class Graph;

    /** A variable that stores node uid */
    size_type nuid_;

    /** Variable to store pointer the graph */
    graph_type *gp_;

    /** Construct a node_iterator */
    NodeIterator(size_type nuid, const graph_type *g)
        : nuid_(nuid), gp_(const_cast<graph_type *>(g)) {}
  };

  /** Return a node iterator pointing to the beginning of the node structure
   * @return node_iterator to the beginning
   */
  node_iterator node_begin() const { return node_iterator(0, this); }

  /** Return a node iterator pointing to one past the end of the node structure
   * @return node_iterator to the end
   */
  node_iterator node_end() const { return node_iterator(num_nodes(), this); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** Dereference the iterator to return the edge it points to
     * @return an object of type edge_type in the graph object
     *
     * @pre iterator != node.edge_end() && iterator != IncidentIterator()
     */
    Edge operator*() const {
      size_type curr = gp_->points[gp_->i2u[node_id_]].adj[neigh_ind_].id;
      return Edge(gp_, gp_->i2u[node_id_], curr);
    }

    /** Increment the iterator so that it points to the next neighboring node
     * @return reference to incremented incident iterator
     *
     * @pre iterator != node.edge_end() && iterator != IncidentIterator()
     */
    IncidentIterator &operator++() {
      this->neigh_ind_++;
      return *this;
    }

    /** Compare equality for two incident iterators
     * @param[in] @a n is the incident_iterator to be compared
     * @return a boolean that is true if the two incident_iterators are true
     */
    bool operator==(const IncidentIterator &n) const {
      return (this->gp_ == n.gp_ && this->node_id_ == n.node_id_ &&
              this->neigh_ind_ == n.neigh_ind_);
    }

    /** Compare not equality for two incident iterators
     * @param[in] @a n is the incident_iterator to be compared
     * @return !(*this == n)
     */
    bool operator!=(const IncidentIterator &n) const { return !(*this == n); }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    /** Id of the node being acted upon. */
    size_type node_id_;

    /** Pointer to the graph object. */
    graph_type *gp_;

    /** Index of the neighbors that we have shifted through */
    size_type neigh_ind_;

    IncidentIterator(size_type ni, const graph_type *g, size_type nei)
        : node_id_(ni), gp_(const_cast<graph_type *>(g)), neigh_ind_(nei) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Dereference the iterator to return the edge it points to
     * @return an object of type edge_type in the graph object
     *
     * @pre iterator != graph.edge_end() && iterator != EdgeIterator()
     */
    Edge operator*() const { return *iit_; }

    /** Increment the iterator so that it points to the next edge
     * @return reference to incremented edge iterator
     *
     * @pre iterator != graph.edge_end() && iterator != EdgeIterator()
     */
    EdgeIterator &operator++() {
      ++iit_;
      fix();
      return *this;
    }

    /** Compare equality for two edge iterators
     * @param[in] @a ei is the edge_iterator to be compared
     * @return a boolean that is true if the two edge_iterators are true
     */
    bool operator==(const EdgeIterator &ei) const {
      return (this->gp_ == ei.gp_ && this->ni_ == ei.ni_ &&
              this->iit_ == ei.iit_);
    }

    /** Compare not equality for two iedge iterators
     * @param[in] @a ei is the edge_iterator to be compared
     * @return !(*this == n)
     */
    bool operator!=(const EdgeIterator &ei) const { return !(*this == ei); }

  private:
    friend class Graph;

    // node iterator
    node_iterator ni_;

    // incident edge iterator
    incident_iterator iit_;

    // variable to store a pointer to the graph
    graph_type *gp_;

    /** Private constructor */
    EdgeIterator(node_iterator ni, incident_iterator iit, const graph_type *gp)
        : ni_(ni), iit_(iit), gp_(const_cast<graph_type *>(gp)) {
      fix();
    }

    void fix() {

      while (ni_ != gp_->node_end()) {
        while (iit_ != (*ni_).edge_end()) {
          if ((*iit_).node1() < (*iit_).node2()) {
            return;
          } else {
            ++iit_;
          }
        }
        ++ni_;
        if (ni_ == gp_->node_end()) {
          break;
        }
        iit_ = (*ni_).edge_begin();
      }

      iit_ = (*(gp_->node_begin())).edge_begin();

      return;
    }
  };

  /** Return a edge iterator pointing to the beginning of the edge structure
   * @return edge_iterator to the beginning
   */
  edge_iterator edge_begin() const {
    return edge_iterator(this->node_begin(),
                         (*(this->node_begin())).edge_begin(), this);
  }

  /** Return a edge iterator pointing to one past the end of the edge structure
   * @return edge_iterator to one past the edge
   */
  edge_iterator edge_end() const {
    return edge_iterator(this->node_end(), (*(this->node_begin())).edge_begin(),
                         this);
  }

  /*
     REMOVE FUNCTIONS
   */

  /** Remove a single node from the graph and edges associated with it
   * @param[in,out] @a n reference to the node to be removed from the graph
   * @return        true (non-zero) value of size_type if the node is successfully 
   *                removed from the graph. false (0) otherwise.
   *
   * @post          All nodes with indices i < @a n.index() are not effected 
   *                (their indices remain unchanged). All nodes with index 
   *                i > @a n.index() are moved to index i - 1.  (if remove is success) 
   * @post          new size() = old size() - 1 (if remove is success)
   * @post          Degree of neighbors of n is reduced by 1. (if remove is success)
   *
   * @note          Invalidates the node @a n. 
   *                At most O(num_nodes()) operations for sparse graphs.
   */
  size_type remove_node(const Node &n) {

    if (!has_node(n)) {
      return false;
    }

    size_type init_ind = n.index();
    // update indices or nodes
    auto it = std::next(node_begin(), n.index());
    points[(*it).uid_].idx = -1;
    ++it;
    for (; it != node_end(); ++it) {
      points[(*it).uid_].idx--;
    }

    // iterate over elements in adj
    for (size_type i = 0; i < n.degree(); i++) {
      size_type neigh_id = points[n.uid_].adj[i].id;
      size_type sz = points[neigh_id].adj.size();

      for (size_type j = 0; j < sz; j++) {
        if (points[neigh_id].adj[j].id == n.uid_) {
          points[neigh_id].adj[j] = points[neigh_id].adj[sz - 1];
          points[neigh_id].adj.pop_back();
          break;
        }
      }
    }

    // erase element from i2u
    i2u.erase(i2u.begin() + init_ind);

    // update number of nodes in graph
    no_nodes--;
    no_edges -= n.degree();
    points[n.uid_].adj.clear();

    return true;
  }

  /** Remove a single node from the graph and edges associated with it
   * @param[in,out] @a n_it pointer to the node to be removed from the graph
   * @return        pointer to the element at the removed node's position if 
   *                the node was removed successfully. else return @a n_it. 
   *
   * @pre           @pre old node_begin() <= @a it < old node_end().
   *
   * @post          All nodes with indices i < @a (*n_it).index() are not 
   *                effected (their indices remain unchanged). All nodes with 
   *                index i > @a (*n_it).index() are moved to index i - 1. 
   *                (if remove is success)
   * @post          new size() = old size() - 1 (if remove is success)
   * @post          Degree of neighbors of (*n_it) is reduced by 1. 
   *                (if remove is success)
   *
   * @note          Invalidates all iterators and @a (*n_it). 
   *                At most O(num_nodes()) operations for sparse graphs.
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  /** Remove a single edge from the graph
   * @param[in,out] @a n1 reference to the node whose edge needs to be removed
   * @param[in,out] @a n2 reference to the node whose edge needs to be removed
   * @return        true (non-zero) value of size_type if the edge is successfully 
   *                removed from the graph. false (0) otherwise. 
   *
   * @post          new num_edges() = old num_edges() - 1 (if remove is success)
   * @post          Degree of neighbors of @a n1 and @a n2 is reduced by 1. 
   *                (if remove is success)
   *
   * @note          At most O(num_nodes()) operations for graphs.
   */
  size_type remove_edge(const Node &n1, const Node &n2) {

    if (!has_edge(n1, n2)) {
      return false;
    }

    // remove from n1
    auto inner_begin = points[n1.uid_].adj.begin();
    auto inner_end = points[n1.uid_].adj.end();
    for (; inner_begin != inner_end; ++inner_begin) {
      if ((*inner_begin).id == n2.uid_) {
        (*inner_begin) = (*(--inner_end));
        points[n1.uid_].adj.pop_back();
        inner_end = points[n1.uid_].adj.end();
        break;
      }
    }

    // remove from n2
    inner_begin = points[n2.uid_].adj.begin();
    inner_end = points[n2.uid_].adj.end();
    for (; inner_begin != inner_end; ++inner_begin) {
      if ((*inner_begin).id == n1.uid_) {
        (*inner_begin) = *(--inner_end);
        points[n2.uid_].adj.pop_back();
        inner_end = points[n2.uid_].adj.end();
        break;
      }
    }

    // update number of edges
    no_edges--;
    return true;
  }

  /** Remove a single edge from the graph
   * @param[in,out] @a e reference to the edge that needs to be removed
   * @return        true (non-zero) value of size_type if the edge is successfully 
   *                removed from the graph. false (0) otherwise. 
   *
   * @post          new num_edges() = old num_edges() - 1 (if remove is success)
   * @post          Degree of neighbors of @a e.node1() and 
   *                @a e.node2() is reduced by 1. 
   *                (if remove is success)
   *
   * @note          At most O(num_nodes()) operations for graphs.
   */
  size_type remove_edge(const Edge &e) {
    node_type n1 = e.node1();
    node_type n2 = e.node2();
    return remove_edge(n1, n2);
  }

  /** Remove a single edge from the graph
   * @param[in,out] @a e_it pointer to the edge to be removed from the graph
   * @return        pointer to the edge at the removed edge's position if 
   *                the edge was removed successfully. else return @a e_it. 
   *
   * @pre           @pre old edge_begin() <= @a it < old edge_end().
   *
   * @post          new num_edges() = old num_edges() - 1 (if remove is success)
   * @post          Degree of neighbors (*e_it).node1() and 
   *                (*e_it).node2() is reduced by 1. 
   *                (if remove is success)
   *
   * @note          Invalidates all iterators except node iterator 
   *                At most O(num_nodes()) operations for graphs.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    edge_type e = *e_it;
    remove_edge(e);
    return e_it;
  }

private:

  // Internal types
  struct edge_wrapper {
    size_type id;
    edge_value_type eval;
  };

  struct node_wrapper {
    Point p;
    node_value_type val;
    std::vector<edge_wrapper> adj;
    size_type idx;
  };

  std::vector<node_wrapper> points;
  std::vector<size_type> i2u;

  size_type no_nodes;
  size_type no_edges;
};

#endif // CME212_GRAPH_HPP