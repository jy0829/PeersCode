#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <utility>
#include <cassert>

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
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  /** Type of the value stored in each node. */
  using node_value_type = V;
  using edge_value_type = E;

  /** Type of this graph. */
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

  /** Type of uid_ of nodes. */
  using uid_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
      : nodes_(), edges_(), i2u_() 
  {
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
      // Nothing needed, the node is invalid
    }

    /** Return this node's position as a const. */
    const Point& position() const {
      return graph_->nodes_[uid_].p_;
    }

    /** Return this node's position. */
    Point& position() {
      return graph_->nodes_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].idx_;
    }

   /** Return this node's value, which can be modified.
    * @return the node's value
    *
    * Complexity: O(1)
    */
    node_value_type& value(){
      return graph_->nodes_[uid_].v_;
    }

   /** Return this node's value, as a const value.
    * @param[out] value  The node's value (const)
    *
    * Complexity: O(1)
    */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].v_;
    }

   /** Return this node's degree in the graph.
    * @return node's degree
    *
    * Complexity: O(1)
    */
    size_type degree() const {
      return graph_->edges_[uid_].size();
    }

   /** Creates an iterator going through all the neighbours of the node.
    * @return an incident_iterator object pointing at the first edge of the node.
    *
    * Complexity: O(1)
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

   /** Returns an iterator pointing at the end of the edges from this node.
    * @return an incident_iterator object 
    *
    * Complexity: O(1)
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Remove this node from the graph.
     */
    size_type remove_node() {
      return graph_->remove_node(*this);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (uid_ == n.uid_);
    }

    /** Test whether this node is less than @a n in a global order.
     * The order is a lexicographic order on (graph_, uid_)
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
      else if (n.graph_ < graph_) {
          return false;
      }
      else {
          return (uid_ < n.uid_);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    /* */
    bool valid() const {
      return uid_ >= 0 && uid_ < graph_->nodes_.size()
          && graph_->nodes_[uid_].idx_ < graph_->i2u_.size()
          && graph_->i2u_[graph_->nodes_[uid_].idx_] == uid_;
    }

    // Pointer back to the Graph.
    Graph *graph_;
    // This node's unique identification number.
    uid_type uid_;
    /** Private Constructor */
    Node(const Graph* graph, uid_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
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
   * @param[in] value    The new node's value, of type node_value_type
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // Add node to nodes_, and return the corresponding Node 
    size_type idx = num_nodes();
    uid_type uid = nodes_.size();

    nodes_.push_back(node_info{idx, position, value});
    edges_.push_back(std::vector<edge_info>());
    i2u_.push_back(uid);

    return Node(this, uid); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // First check if a node is valid, and then if its graph is the same as this
    return this == n.graph_ && n.valid();
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

  /** Remove the node @a n from the graph if it exists.
   * @param[in] n The node to delete
   * @return the number of deleted nodes (0 or 1)
   *
   * @pre old num_nodes() >= 1
   * @pre 0 <= n.index() < num_nodes()
   * @post new num_nodes() == old num_nodes() - 1 
   *
   * @validity: invalidates any iterator on the nodes or edges
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  size_type remove_node(const Node& n) {
    if (has_node(n)) {
      // First remove every edge attached to node n
      while (n.edge_begin() != n.edge_end()) {
        remove_edge(*(n.edge_begin()));
      }

      // We need to invalidate n, by removing its uid from i2u_
      auto it = std::find(i2u_.begin(), i2u_.end(), n.uid_);
      assert(it != i2u_.end());  // should always find n.uid_
      i2u_.erase(it);

      // Now, for every uid after it, we decrement the index
      for (; it != i2u_.end(); ++it) {
        --nodes_[*it].idx_;
      }
      return 1;
    }
    else {
      return 0;
    }
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
      // Nothing to do here, invalid edge
    }

   /** Return this edge's length
    * @return the L2 distance between the two nodes of this edge.
    *
    * Complexity: O(1) amortized operations.
    */
    double length() const{
      return norm(node1().position() - node2().position());
    }

   /** Return this edge value, which can be modified.
    * @return the edge value, of type edge_value_type
    *
    * Complexity: O(1) amortized operations.
    */
    edge_value_type& value() {
      auto comp = [&] (edge_info e_info) {
        return uid2_ == e_info.uid_;
      };
      auto it = std::find_if(graph_->edges_[uid1_].begin(), graph_->edges_[uid1_].end(), comp);
      return it->v_;
    }

   /** Return this edge value, which can't be modified.
    * @return the edge value (const), of type edge_value_type
    *
    * Complexity: O(1) amortized operations.
    */
    const edge_value_type& value() const {
      auto comp = [&] (edge_info e_info) {
        return uid2_ == e_info.uid_;
      };
      auto it = std::find_if(graph_->edges_[uid1_].begin(), graph_->edges_[uid1_].end(), comp);
      return it->v_;
    }

    /** Return the first node of this Edge */
    Node node1() const {
      return Node(graph_, uid1_);
    }

    /** Return the second node of this Edge */
    Node node2() const {
      return Node(graph_, uid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges are in the same graph and
     * represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph_ == e.graph_
          && ((uid1_ == e.uid1_ && uid2_ == e.uid2_) 
           || (uid2_ == e.uid1_ && uid1_ == e.uid2_));
    }

    /** Test whether this edge is less than @a e in a global order.
     * The ordering is a lexicographic order on (graph, node1().index(), node2().index())
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ < e.graph_) {
        return true;
      }
      else if (e.graph_ < graph_) {
        return false;
      }
      else {
        if (uid1_ < e.uid1_) {
          return true;
        }
        else if (e.uid1_ < uid1_) {
          return false;
        }
        else {
          return (uid2_ < e.uid2_);
        }
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph 
    Graph *graph_;
    // The node 1 (origin node) uid
    uid_type uid1_;
    // The node 2 (destination node) uid
    uid_type uid2_;

    /** Private Constructor */
    Edge(const Graph* graph, uid_type uid1, uid_type uid2)
        : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(num_nodes())
   */
  size_type num_edges() const {
    size_type count = 0;
    for (uid_type uid : i2u_) {
      count += edges_[uid].size();
    }
    assert(count % 2 == 0);
    return count / 2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  Edge edge(size_type i) const {
    EdgeIterator it_start = edge_begin();
    for (size_type j = 0; j < i; j++) {
      ++it_start;
    }
    return *it_start;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(a.degree())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Only check for (a, b), because (b, a) exists iff (a, b) exists
    assert(a.valid());
    assert(b.valid());
    auto comp = [&b] (edge_info e_info) {
      return e_info.uid_ == b.uid_;
    };
    return (std::find_if(edges_[a.uid_].begin(), edges_[a.uid_].end(), comp) != edges_[a.uid_].end());
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
   * Complexity: O(a.degree()) amortized operations
   */
  Edge add_edge(const Node& a, const Node& b) {
    // We will add both (a, b) and (b, a) in the graph.
    if (has_edge(a, b)) {
      // We found the edge
      return Edge(this, a.uid_, b.uid_);
    }
    else {
      edges_[a.uid_].push_back(edge_info{b.uid_, edge_value_type()});
      edges_[b.uid_].push_back(edge_info{a.uid_, edge_value_type()});  // also add the edge (b, a), which SHOULD NOT exist
      return Edge(this, a.uid_, b.uid_);
    }
  }

  /** Remove the edge (@a n1, @a n2) from the graph.
   * @param[in] n1 The first node of the edge to delete
   * @param[in] n2 The second node of the edge to delete
   * @return the number of deleted edges (0 or 1)
   *
   * @pre old num_edges() >= 1
   * @post new num_edges() == old num_edges() - 1 
   *
   * @validity: invalidates any iterator on the nodes or edges
   *
   * Complexity: O(n1.degree() + n2.degree())
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (has_edge(n1, n2)) {
      auto comp1 = [&n2] (edge_info e_info) {
        return e_info.uid_ == n2.uid_;
      };
      auto n2_it = std::find_if(edges_[n1.uid_].begin(), edges_[n1.uid_].end(), comp1);
      edges_[n1.uid_].erase(n2_it);

      auto comp2 = [&n1] (edge_info e_info) {
        return e_info.uid_ == n1.uid_;
      };
      auto n1_it = std::find_if(edges_[n2.uid_].begin(), edges_[n2.uid_].end(), comp2);
      edges_[n2.uid_].erase(n1_it);

      return 1;
    }
    else {
      return 0;
    }
  }

  /** Remove the edge @a e from the graph.
   * @param[in] e The edge to delete
   * @return the number of deleted edges (0 or 1)
   *
   * @pre old num_edges() >= 1
   * @post new num_edges() == old num_edges() - 1 
   *
   * @validity: invalidates any outstanding node or edge
   *            invalidates any iterator on the nodes or edges
   *
   * Complexity: O(n1.degree() + n2.degree())
   */
  size_type remove_edge(const Edge& e) {
    Node n1 = e.node1();
    Node n2 = e.node2();

    return remove_edge(n1, n2);
  }

  /** Remove the edge where @a e_it is pointing at from the graph.
   * @param[in] e An edge_iterator pointing to the edge to delete.
   * @return An edge_iterator pointing to the previous next edge.
   *
   * @pre old num_edges() >= 1
   * @post new num_edges() == old num_edges() - 1 
   *
   * @validity: invalidates any outstanding node or edge
   *            invalidates any iterator on the nodes or edges
   *
   * Complexity: O(n1.degree() + n2.degree())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    size_type e_id = e.index();
    remove_edge(e);
    return edge_begin() + e_id;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    i2u_.clear();
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

   /** Deferencement operator for NodeIterator
    * @return Node to which the NodeIterator points to
    *
    * Complexity: 0(1)
    */
    Node operator*() const
    {
      return graph_->node(idx_);
    }

   /** Increment operator for NodeIterator
    * @return NodeIterator pointing to the next node
    *
    * Complexity: 0(1)
    */
    NodeIterator& operator++()
    {
      ++idx_;
      return *this;
    }

   /** Equality operator for two NodeIterators
    * @return True if the iterators refer to the same value
    *
    * Complexity: 0(1)
    */
    bool operator==(const NodeIterator& it) const
    {
      return (graph_ == it.graph_) && (idx_ == it.idx_);
    }

   private:
    friend class Graph;
    
    // Pointer back to the Graph 
    Graph *graph_;
    // index of the current Node we are pointing to
    size_type idx_;

    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {
    }
    
  };


 /** Creates an iterator from the beginning of the nodes
  * @return a node_iterator object pointing to the first Node
  *
  * Complexity: 0(1)
  */
  node_iterator node_begin() const
  {
      return NodeIterator(this, 0);
  }

 /** Creates an iterator at the end of the nodes
  * @return a node_iterator object pointing to the first Node
  *
  * Complexity: 0(1)
  */
  node_iterator node_end() const
  {
      return NodeIterator(this, num_nodes());
  }

  /** Erases the node where the iterator is pointing at
   *  and return an iterator to the next node.
   * @param[in] n_it Node iterator pointing to the node to delete 
   * @return a node_iterator to the new location of the next node
   *
   * @pre old num_nodes() >= 1
   * @post new num_nodes() == old num_nodes() - 1 
   *
   * @validity: invalidates any outstanding edge or node object
   *            invalidates any iterator on the nodes or edges
   *
   * Complexity: O(sum degrees of neighbours of n + num_nodes())
   *             <= O(num_nodes() + num_edges())
   */
  node_iterator remove_node(node_iterator n_it) {
    node_type n = *n_it;
    remove_node(n);
    return NodeIterator(this, n.index());
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


   /** Deferencement operator for IncidentIterator
    * @return Edge to which the IncidentIterator points to
    *
    * Complexity: 0(1)
    */
    Edge operator*() const
    {
      return Edge(graph_, uid1_, graph_->edges_[uid1_][idx2_].uid_);
    }

   /** Increment operator, goes to the next node in the adjacency list of node 1
    * @return an IncidentIterator object pointing to the next edge
    *
    * Complexity: 0(1)
    */
    IncidentIterator& operator++()
    {
      ++idx2_;
      return *this;
    }

   /** Equality operator for two NodeIterators
    * @return True if the iterators refer to the same value
    *
    * Complexity: 0(1)
    */
     bool operator==(const IncidentIterator& it) const
     {
       return (graph_ == it.graph_) && (uid1_ == it.uid1_) && (idx2_ == it.idx2_);
     }

   private:
    friend class Graph;
      // Pointer back to the Graph 
      Graph *graph_;
      // uid of the node1 of the edge
      uid_type uid1_;
      // index of the node2 of the edge in graph_->edges_[uid1_]
      size_type idx2_;

      /** Private Constructor */
      IncidentIterator(const Graph* graph, uid_type uid1, size_type idx2)
          : graph_(const_cast<Graph*>(graph)), uid1_(uid1), idx2_(idx2) {
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

   /** Dereferences the EdgeIterator and returns the current Edge
    * @pre the iterator points to a valid edge, i.e. an edge with (node1.index() < node2.index())
    *
    * @return the Edge at the EdgeIterator location
    */
    Edge operator*() const
    {
      return Edge(graph_, uid1(), uid2());
    }

   /** Increment the iterator until we find an edge (node1, node2)
    *  with node1.index() < node2.index()
    *
    * @pre the iterator is not at the end of the end of @a edges_
    *      i.e. @a node1_ != @a graph_->edges_.size() 
    * @pre if the condition above is satisfied, @a node2_id_ != @a graph_->edges_[@a node1_].size()
    * @post *this points to a valid edge (node1 < node2) or to the end of the edges_
    */
    EdgeIterator& operator++() {
      do {
        ++idx2_;
        if (idx2_ >= graph_->edges_[uid1()].size()) {
          ++idx1_;
          idx2_ = 0;
        }
      } while ((idx1_ != graph_->num_nodes())
            && ((graph_->edges_[uid1()].size() == 0)
             || (uid1() >= uid2())));
      return *this;
    }

   /** Check if two EdgeIterators are pointing at the same location
    * @a x: EdgeIterator to compare with
    *
    */
    bool operator==(const EdgeIterator& x) const {
      return (graph_ == x.graph_) && (idx1_ == x.idx1_) && (idx2_ == x.idx2_);
    }

   private:
    friend class Graph;

    uid_type uid1() const {
      return graph_->i2u_[idx1_];
    }
    
    uid_type uid2() const {
      return graph_->edges_[uid1()][idx2_].uid_;
    }
    
    // Pointer back to the Graph 
    Graph *graph_;
    // index of the first node we are pointing to
    size_type idx1_;
    // index of the second node in the neighbours of node_1
    // i.e. an id in graph_->edges_[i2u_[idx1_]]
    size_type idx2_;;

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type idx1, size_type idx2)
        : graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2) {
    }

  };

 /** Creates an EdgeIterator pointing to the first (valid) edge.
  *
  */
  edge_iterator edge_begin() const
  {
    if (num_edges() == 0) return edge_end();
    else {
      for (auto it = edges_.begin(); it != edges_.end(); ++it) {
        if (it->size() > 0) {
          return EdgeIterator(this, std::distance(edges_.begin(), it), 0);
        }
      }
      assert(false);
      return EdgeIterator(this, 0, 0);  // should not get here
    }
  }

 /** Creates an EdgeIterator pointing to the end of edges.
  *
  */
  edge_iterator edge_end() const
  {
    return EdgeIterator(this, num_nodes(), 0);
  }

 private:
  // contains info for nodes
  struct node_info {
    size_type idx_;      // index in nodes_ or edges_
    Point p_;            // position
    node_value_type v_;  // value
  };

  // contains info for edges: uid_ of destination node and value
  struct edge_info {
    uid_type uid_;
    edge_value_type v_;
  };

  // nodes_ contains the Points and values corresponding to every node
  std::vector <node_info> nodes_;
  // edges_ contains the neighbours of each node, and values for the edges
  std::vector <std::vector <edge_info> > edges_;

  // i2u_: size num_nodes(), gives the uid of a node given its index.
  std::vector <uid_type> i2u_;

};


#endif // CME212_GRAPH_HPP
