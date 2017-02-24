#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iterator>
#include <vector>
#include <iostream>

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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this
   * graph. */
  using node_value_type = V;
  using edge_value_type = E;
  using graph_type = Graph;

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
  Graph(): edge_count(0), nodes_(), adjacency_(), adj_edge_val_(), i2u_() {
    // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].p;
    }

    Point& position() {
      return graph_->nodes_[uid_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].uid;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Test whether this node and @a n are equal.
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_) and (uid_ == n.uid_);
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
      if (graph_ == n.graph_) {
        return uid_ < n.uid_;
      }
      return (graph_ < n.graph_);
    }

    /** Get the value of this node where user is allowed to modify @a value in this node
     * @return node_value_type value of this node as a reference.
     *
     * Complexity O(1).
     */
    node_value_type& value() {
      return graph_->nodes_[uid_].val;
    }

    /** Get the value of this node but user is not allowed to modify it
     * @return node_value_type value of this node as a const reference.
     *
     * Complexity O(1).
     */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].val;
    }

    /** Get the degree of the node
     * @return the number of incident edges connected to the current node
     * @post 0 <= degree() < num_nodes().
     *
     * Complexity O(1)
     */
    size_type degree() const {
      return graph_->adjacency_[uid_].size();
    }

    /** Return an Incident Iterator pointing to the first edge of this graph
     * If the container is empty, the returned value shall not be dereferenced.
     *
     * Complexity O(1)
     */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, uid_, 0);
    }

    /** Return an IncidentIterator that points to the past-the-end element in adjacency container of this graph.
     * If the container is empty, the returned value compares equal to the one returned by begin with the same argument.
     *
     * Complexity O(1)
     */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, uid_, degree());
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    Node(const Graph* graph, size_type uid) : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    node_info new_node;
    new_node.uid = i2u_.size();
    new_node.p = position;
    new_node.val = val;
    nodes_.push_back(new_node);
    i2u_.push_back(nodes_.size() - 1);
    adjacency_.resize(nodes_.size());
    adj_edge_val_.resize(nodes_.size());
    return Node(this, nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    return (this == n.graph_) and (n.uid_ < nodes_.size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //assert(i < node_count);
    (void) i;             // Quiet compiler warning
    return Node(this, i2u_[i]);
  }

  /** Remove the specified node as well as all associated edges with this node from the graph
   * param[in] n Node of this Graph to be removed.
   * @pre has_node(@a n) == true.
   * @pre g.node(i).index() == i for all i with 0 ≤ i < g.num nodes().
   * @post g.node(i).index() == i for all i with 0 ≤ i < g.num nodes().
   * @post g.node(@a n.index()) == @a n.
   * @post new num_nodes() = old num_nodes() - 1.
   * @post has_node(@a n) == false.
   * @post All edges @a e that have e.node1() == @a n or e.node2() == @a n become invalid and get removed.
   * @post new num_edges() = old num_edges() - n.degree().
   * @post Invalidates all outstanding Iterators (NodeIter, EdgeIter, IncidentIter)
   * @post Invalidates IncidentIterators that iterate on incident edges of node @a n.
   * @return The index of the removed node @a n
   *
   * Complexity: O(degree()^2) on average, O(num_nodes()^2) in the worst case.
   */
  size_type remove_node(const Node& n) {
    size_type ind;
    for(ind = 0; ind < i2u_.size(); ++ind) {
      if(i2u_[ind] == n.uid_) break;
    }
    i2u_.erase(i2u_.begin() + ind);

    while(adjacency_[n.uid_].size() > 0) {
      remove_edge(n, Node(n.graph_, adjacency_[n.uid_][0]));
    }

    for(size_type i = ind; i < i2u_.size(); ++i) {
      (nodes_.at(i2u_[i])).uid = i;
    }
    return ind;
  }

  /** Remove the specified node as well as all associated edges with this node from the graph using NodeIter
   * param[in] it NodeIterator of this Graph pointing to the node to be removed.
   * @pre has_node(@a n) == true.
   * @pre g.node(i).index() == i for all i with 0 ≤ i < g.num nodes().
   * @post g.node(i).index() == i for all i with 0 ≤ i < g.num nodes().
   * @post g.node(@a n.index()) == @a n.
   * @post new num_nodes() = old num_nodes() - 1.
   * @post has_node(@a n) == false.
   * @post All edges @a e that have e.node1() == @a n or e.node2() == @a n become invalid and get removed.
   * @post new num_edges() = old num_edges() - n.degree().
   * @post Invalidates all outstanding Iterators (NodeIter, EdgeIter, IncidentIter)
   * @post Invalidates IncidentIterators that iterates on incident edges of node @a n.
   * @return The NodeIter pointing to the index where the node was before it is removed in the i2u_.
   *
   * Complexity: O(degree()^2) on average, O(num_nodes()^2) in the worst case.
   */

  node_iterator remove_node(node_iterator it) {
    int ind = remove_node(*it);
    return NodeIterator(this, ind);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool same_graph = (graph_ == e.graph_);
      bool direction1 = ((node1_uid_ == e.node1_uid_) and (node2_uid_ == e.node2_uid_));
      bool direction2 = ((node1_uid_ == e.node2_uid_) and (node2_uid_ == e.node1_uid_));
      return same_graph and (direction1 or direction2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_) {
        if (node1_uid_ == e.node1_uid_) {
          return (node2_uid_ < e.node2_uid_);
        }
        return (node1_uid_ < e.node1_uid_);
      }
      return (graph_ < e.graph_);
    }

    double length() const {
      return norm(node1().position() - node2().position());
    }

    edge_value_type& value() {
      Edge edge = Edge(graph_, fmin(node1_uid_, node2_uid_), fmax(node1_uid_, node2_uid_));
      size_type e_ind = 0;
      for(auto it = edge.node1().edge_begin(); it != edge.node1().edge_end(); ++it) {
        if((*it).node2_uid_ == edge.node2_uid_) {
          e_ind = it.edge_index_;
          break;
        }
      }
      return graph_->adj_edge_val_[edge.node1_uid_][e_ind];
    }

    const edge_value_type& value() const {
      Edge edge = Edge(graph_, fmin(node1_uid_, node2_uid_), fmax(node1_uid_, node2_uid_));
      size_type e_ind = 0;
      for(auto it = edge.node1().edge_begin(); it != edge.node1().edge_end(); ++it) {
        if((*it).node2_uid_ == edge.node2_uid_) {
          e_ind = it.edge_index_;
        }
      }
      return graph_->adj_edge_val_[edge.node1_uid_][e_ind];
    }


  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type node1_uid_;
    size_type node2_uid_;

    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
        : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    //assert(i < edge_count);
    (void) i;             // Quiet compiler warning
    EdgeIterator ei = edge_begin();
    while (i-- > 0)
      ++ei;
    return *ei;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    for (size_type j = 0; j < adjacency_[a.uid_].size(); j++) {
      if (adjacency_[a.uid_][j] == b.uid_) {
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
    // HW0: YOUR CODE HERE
    (void) a, (void) b;   // Quiet compiler warning
    if (!has_edge(a, b)) {
      adjacency_[a.uid_].push_back(b.uid_);
      adjacency_[b.uid_].push_back(a.uid_);
      adj_edge_val_[a.uid_].push_back(edge_value_type());
      adj_edge_val_[b.uid_].push_back(edge_value_type());
      ++edge_count;
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove the specified edge from this Graph.
   * @param[in] n1 Nodes on one side of the Edge to be removed.
   * @param[in] n2 Nodes on the other side of an Edge to be removed.
   * @pre Nodes @a n1 and @a n2 are valid.
   * @post If has_Edge(n1, n2) == True, new num_edge() = old num_edge() - 1.
   * @post If has_Edge(n1, n2) == False, new num_edge() = old num_edge().
   * @post Invalidates all outstanding Iterators (EdgeIter, IncidentIter)
   * @return The index where Node @ n2 was on the adjacency list before it was removed if the has_Edge(n1, n2) == False
   * return @a n1's old num_degree()
   *
   * Complexity: O(degree()).
   */

  size_type remove_edge(const Node& n1, const Node& n2) {
    size_type ind = n1.degree();
    for(auto it = n1.edge_begin(); it != n1.edge_end(); ++it) {
      if((*it).node2_uid_ == n2.uid_) {
        ind = it.edge_index_;
        adjacency_[n1.uid_].erase(adjacency_[n1.uid_].begin() + ind);
        adj_edge_val_[n1.uid_].erase(adj_edge_val_[n1.uid_].begin() + ind);
        --edge_count;
        break;
      }
    }
    for(auto it = n2.edge_begin(); it != n2.edge_end(); ++it) {
      if((*it).node2_uid_ == n1.uid_) {
        adjacency_[n2.uid_].erase(adjacency_[n2.uid_].begin() + it.edge_index_);
        adj_edge_val_[n2.uid_].erase(adj_edge_val_[n2.uid_].begin() + it.edge_index_);
        break;
      }
    }
    return ind;
  }

  /** Remove the specified edge from this Graph.
  * @param[in] e Edge of this graph to be removed.
  * @pre Edge @a e is valid
  * @post If has_Edge(@a e.node1(), @a e.node2()) == True, new num_edge() = old num_edge() - 1.
  * @post If has_Edge(@a e.node1(), @a e.node2()) == False, new num_edge() = old num_edge().
  * @post Invalidates all outstanding Iterators (EdgeIter, IncidentIter)
  * @return The index where Node @a e.node2() was on the adjacency list before it was removed if the has_Edge(n1, n2) == False
  * return @a e.node1()'s old num_degree()
  *
  * Complexity: O(degree()).
  */

  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove the Edge that the specified EdgeIter points to from this Graph.
 * @param[in] e_it EdgeIterator pointing to the Edge of this graph we want to remove
 * @pre EdgeIter @a e_it is valid
 * @post If has_Edge(@a (*e_it).node1(), @a (*e_it).node2()) == True, new num_edge() = old num_edge() - 1.
 * @post If has_Edge(@a (*e_it).node1(), @a (*e_it).node2()) == False, new num_edge() = old num_edge().
 * @post Invalidates all outstanding Iterators (EdgeIter, IncidentIter)
 * @return The EdgeIter pointing to where Node @a (*e_it).node2()  was on the adjacency list before it was removed if has_Edge(@a (*e_it).node1(), @a (*e_it).node2()) == False
 * return @a (*e_it).node1() 's old num_degree()
 *
 * Complexity: O(degree()).
 */

  edge_iterator remove_edge(edge_iterator e_it) {
    return EdgeIterator(this, (*e_it).node1().uid_, remove_edge((*e_it).node1(), (*e_it).node2()));
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    adjacency_.clear();
    adj_edge_val_.clear();
    i2u_.clear();
    edge_count = 0;

  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

      /** Dereference the NodeIterator to get the Node that the iterator points to
       * @pre 0 <= index_ < graph_->size(), NodeIterator != node_end().
       * @return the valid Node this iterator poiting to
       *
       * Complexity: O(1).
       */
      Node operator*() const{
        return graph_->node(index_);
      }

      /** Move the index where the current NodeIterator is pointing at forward by 1
       * @post NodeIterator points to the next node
       * @return the reference of the NodeIterator pointing to the node that it would be pointing to if advanced 1 position.
       *
       * Complexity: O(1).
       */
      NodeIterator& operator++(){
        if (index_ < graph_->size())
          ++index_;
        return *this ;
      }

      /** Test whether this NodeIterator and NodeIterator @a i are equal.
       * Equal NodeIterators represent the same graph and the same indices they are pointing at.
       * @param[in] i NodeIterator in graph.
       * @return True if this NodeIterator and @a i are the same, otherwise False
       *
       * Complexity: O(1).
       */
      bool operator==(const NodeIterator& i) const{
        return ((graph_ == i.graph_) && (index_ == i.index_));
      }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type index_;
    NodeIterator(const Graph* g, size_type idx) : graph_(const_cast<Graph*>(g)), index_(idx) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
    /** Return a NodeIterator pointing to the first node of this graph
     * If the node container is empty, the returned value shall not be dereferenced.
     * It returns node_end() if num_nodes() == 0.
     */
    node_iterator node_begin() const{
      return node_iterator(this, 0);
    }
    /** Return a NodeIterator that points to the past-the-end element in the node container of this graph.
     * If the node container is empty, the returned value compares equal to the one returned by begin with the same argument.
     */
    node_iterator node_end() const{
      return node_iterator(this, size());
    }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>  {
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference the IncidentIterator to get the Edge that the iterator points to
     * @pre 0 <= node_uid_ < graph_->adjacency_.size(), 0 <= edge_index_ < graph_->adjacency_[node_uid_].size()
     * @return the incident Edge this iterator pointing to
     *
     * Complexity: O(1).
     */
    Edge operator*() const{
      size_type node2 = graph_->adjacency_[node_uid_][edge_index_];
      return Edge(graph_, node_uid_ , node2);
    }

    /** Move the index where the current IncidentIterator is pointing at forward by 1
      * @post IncidentIterator points to the next valid one.
      * @return the IncidentIterator pointing to the edge connected to the same node that it would be
      * pointing to if advanced 1 position.
      *
      * Complexity: O(1).
      */
    IncidentIterator& operator++(){
      if (edge_index_ < graph_->adjacency_[node_uid_].size())
        ++edge_index_;
      return *this;
    }

    /** Test whether this IncidentIterator and IncidentIterator @a i are equal.
     * Equal IncidentIterators represent the same graph and the same node ids and edges indices they are pointing at.
     * @param[in] i IncidentIterator in graph.
     * @return True if this IncidentIterator and @a i are the same, otherwise False
     *
     * Complexity: O(1).
     */
    bool operator==(const IncidentIterator& i) const{
      return ((graph_ == i.graph_) && (node_uid_ == i.node_uid_) && (edge_index_ == i.edge_index_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_uid_;  //which node we are currently iterating over
    size_type edge_index_; //the edge index we are at from that node
    IncidentIterator(const Graph* g, size_type node_uid, size_type index) :
        graph_(const_cast<Graph*>(g)), node_uid_(node_uid), edge_index_(index) {
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference the IncidentIterator to get the Edge that the iterator points to
     * @pre 0 <= src_uid_ < graph_->adjacency_.size(), 0 <= dst_ind_ < graph_->adjacency_[src_uid_].size()
     *       graph_->adjacency_[src_uid_][dst_ind_] < src_uid_
     * @return the Edge that this iterator currently pointing at
     *
     * Complexity: O(1).
     */
    Edge operator*() const{
      return Edge(graph_, src_uid_, graph_->adjacency_[src_uid_][dst_ind_]);
    }

    /** Move the index where the current EdgeIterator is pointing at forward by 1
     * @post EdgeIterator points to the next valid edge ordered bu src_uid then dst_ind.
     * @return the EdgeIterator pointing to the next edge that it would be pointing to if advanced 1 position.
     *
     * Complexity: O(1).
     */
    EdgeIterator& operator++(){
      if (src_uid_ < graph_->adjacency_.size()){
        ++dst_ind_;
      }
      get_valid_edge();
      return *this;
    }

    /** Test whether this EdgeIterator and EdgeIterator @a i are equal.
     * Equal EdgeIterators represent the same graph and the same src_uid and dst_ind_.
     * @param[in] i EdgeIterator in graph.
     * @return True if this EdgeIterator and @a i are the same, otherwise False
     *
     * Complexity: O(1).
     */
    bool operator==(const EdgeIterator& i) const{
      return ((graph_ == i.graph_) && (src_uid_ == i.src_uid_) && (dst_ind_ == i.dst_ind_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type src_uid_;
    size_type dst_ind_;
    EdgeIterator(const Graph *g, size_type uid, size_type index) : graph_(const_cast<Graph*>(g)), src_uid_(uid), dst_ind_(index) {
      get_valid_edge();
    }

    /** Move the current EdgeIterator forward until it hits the next valid edge and make sure that we do not double count any
     * of them by pointing to only the edges that has the id of the source node higher than the dst node.
     */
    void get_valid_edge(){
      long num_adj = graph_->adjacency_.size();
      if (src_uid_ >= num_adj)
        return;
      if (graph_->adjacency_[src_uid_].size() <= dst_ind_) {
        ++src_uid_;
        dst_ind_ = 0;
      }
      while((src_uid_ < num_adj) && (graph_->adjacency_[src_uid_].size() == 0)){
        ++src_uid_;
      }
      while((src_uid_ < num_adj) && (graph_->adjacency_[src_uid_][dst_ind_] < src_uid_)) {
        ++dst_ind_;
        if (dst_ind_ >= graph_->adjacency_[src_uid_].size()) {
          ++src_uid_;
          dst_ind_ = 0;
        }
        while((src_uid_< num_adj) && (graph_->adjacency_[src_uid_].size() == 0)){
          ++src_uid_;
        }
      }
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return an EdgeIterator pointing to the first edge of this graph
   * If the container is empty, the returned value shall not be dereferenced.
   */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0, 0);
  }

  /** Return an EdgeIterator that points to the past-the-end element in the adjacency container of this graph.
   * If the adjacency is empty, the returned value compares equal to the one returned by begin with the same argument.
   */
  edge_iterator edge_end() const {
    return edge_iterator(this, adjacency_.size(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct node_info
  {
    Point p;
    size_type uid;
    node_value_type val;
  };


  size_type edge_count;
  std::vector<node_info> nodes_;
  std::vector<std::vector<size_type>> adjacency_;
  std::vector<std::vector<edge_value_type>> adj_edge_val_;
  std::vector<size_type> i2u_;

};

#endif // CME212_GRAPH_HPP
