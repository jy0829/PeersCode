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
template <typename V>
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
  
  /** Type of the value stored in each node. */
  using node_value_type = V;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
      : nodes_(), nodes_values_(), edges_(), node_size_(0)
  {
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
      // Nothing needed, the node is invalid
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1 #1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for those functions.

   /** Return this node's value, which can be modified.
    * @return the node's value
    *
    * Complexity: O(1) amortized operations.
    */
    node_value_type& value(){
      return graph_->nodes_values_[uid_];
    }

   /** Return this node's value, as a const value.
    * @param[out] value  The node's value (const)
    *
    * Complexity: O(1) amortized operations.
    */
    const node_value_type& value() const {
      return graph_->nodes_values_[uid_];
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

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_) && (uid_ == n.uid_);
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
      // HW0: YOUR CODE HERE
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
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer back to the Graph.
    Graph *graph_;
    // This node's unique identification number.
    size_type uid_;
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value    The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Add node to nodes_, and return the corresponding Node while incrementing node_size_
    nodes_.push_back(position);
    nodes_values_.push_back(value);
    edges_.push_back(std::vector<size_type>());
    return Node(this, node_size_++); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // Nodes are never deleted here, so a node is in the graph if its id
    // is less than the node size
    return (this == n.graph_) && (n.uid_ < node_size_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
      // Nothing to do here, invalid edge
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((graph_ == e.graph_) && (node1_ == e.node1_) && (node2_ == e.node2__)) ||
             ((graph_ == e.graph_) && (node1_ == e.node2_) && (node1_ == e.node2__));
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
      else if (e.graph_ < graph_) {
        return false;
      }
      else {
        if (node1_ < e.node1_) {
          return true;
        }
        else if (e.node1_ < node1_) {
          return false;
        }
        else {
          return (node2_ < e.node2_);
        }
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    
    // Pointer back to the Graph 
    Graph *graph_;
    // The node 1 (origin node) uid_
    size_type node1_;
    // The node 2 (destination node) uid_
    size_type node2_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type node1, size_type node2)
        : graph_(const_cast<Graph*>(graph)), node1_(node1), node2_(node2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    EdgeIterator it_start = edge_begin();
    EdgeIterator it_end = edge_end();
    for (size_type j = 0; j < i; j++) {
      ++it_start;
    }
    return *it_start;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    return (std::find(edges_[a.uid_].begin(), edges_[a.uid_].end(), b.uid_) != edges_.end());
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
    // We will add both (a, b) and (b, a) in the graph.
    auto it = std::find(edges_[a.uid_].begin(), edges_[a.uid_].end(), b.uid_);
    if (it != edges_[a.uid_].end()) {
      // We found the edge
      return Edge(this, a.uid_, b.uid_);
    }
    else {
      edges_[a.uid_].push_back(b.uid_);
      edges_[b.uid_].push_back(a.uid_);  // also add the edge (b, a), which SHOULD NOT exist
      ++edge_size_;
      return Edge(this, a.uid_, b.uid_);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    nodes_values_.clear();
    edges_.clear();
    node_size_ = 0;
    edge_size_ = 0;
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

   /** Deferencement operator for NodeIterator
    * @return Node to which the NodeIterator points to
    *
    * Complexity: 0(1)
    */
    Node operator*() const
    {
      return graph_->node(uid_);
    }

   /** Increment operator for NodeIterator
    * @return NodeIterator pointing to the next node
    *
    * Complexity: 0(1)
    */
    NodeIterator& operator++()
    {
      ++uid_;
      return *this;
    }

   /** Equality operator for two NodeIterators
    * @return True if the iterators refer to the same value
    *
    * Complexity: 0(1)
    */
    bool operator==(const NodeIterator& it) const
    {
      return (graph_ == it.graph_) && (uid_ == it.uid_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    
    // Pointer back to the Graph 
    Graph *graph_;
    // id of the current Node we are pointing to
    size_type uid_;

    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

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
      return NodeIterator(this, node_size_);
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

   /** Deferencement operator for IncidentIterator
    * @return Edge to which the IncidentIterator points to
    *
    * Complexity: 0(1)
    */
    Edge operator*() const
    {
      return Edge(graph_, node1_, graph_->edges_[node1_][node2_id_]);
    }

   /** Increment operator, goes to the next node in the adjacency list of node 1
    * @return an IncidentIterator object pointing to the next edge
    *
    * Complexity: 0(1)
    */
    IncidentIterator& operator++()
    {
      ++node2_id_;
      return *this;
    }

   /** Equality operator for two NodeIterators
    * @return True if the iterators refer to the same value
    *
    * Complexity: 0(1)
    */
     bool operator==(const IncidentIterator& it) const
     {
       return (graph_ == it.graph_) && (node1_ == it.node1_) && (node2_id_ == it.node2_id_);
     }

   private:
    friend class Graph;
      // HW1 #3: YOUR CODE HERE
      // Pointer back to the Graph 
      Graph *graph_;
      // id of the node1 of the edge
      size_type node1_;
      // id of the node2 of the edge in graph_->edges_[node1_]
      size_type node2_id_;

      /** Private Constructor */
      IncidentIterator(const Graph* graph, size_type node1, size_type node2_id)
          : graph_(const_cast<Graph*>(graph)), node1_(node1), node2_id_(node2_id) {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

   /** Dereferences the EdgeIterator and returns the current Edge
    * @pre the iterator points to a valid edge, i.e. an edge with (node1.index() < node2.index())
    *
    * @return the Edge at the EdgeIterator location
    */
    Edge operator*() const
    {
      return Edge(graph_, node1_, graph_->edges_[node1_][node2_]);
    }

   /** Increment the iterator until we find an edge (node1, node2)
    *  with node1.index() < node2.index()
    *
    * @pre the iterator is not at the end of the end of @a edges_
    *      i.e. @a node1_ != @a graph_->edges_.size() 
    * @pre if the condition above is satisfied, @a node2_ != @a graph_->edges_[@a node1_].size()
    * @post *this points to a valid edge (node1 < node2) or to the end of the edges_
    */
    EdgeIterator& operator++() {
      do {
        ++node2_;
        if (node2_ == graph_->edges_[node1_].size()) {
          ++node1_;
          node2_ = 0;
        }
      } while ((node1_ != graph_->edges_.size()) && (node1_ >= graph_->edges_[node1_][node2_]));
      return *this;
    }

   /** Check if two EdgeIterators are pointing at the same location
    * @a x: EdgeIterator to compare with
    *
    */
    bool operator==(const EdgeIterator& x) const {
      return (graph_ == x.graph_) && (node1_ == x.node1_) && (node2_ == x.node2_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
    // Pointer back to the Graph 
    Graph *graph_;
    // id of the first node we are pointing to
    size_type node1_;
    // id of the second node in the neighbours of node_1
    // i.e. an id in graph_->edges_[node1_]
    size_type node2_;

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type node1, size_type node2)
        : graph_(const_cast<Graph*>(graph)), node1_(node1), node2_(node2) {
      // we want to point to the first valid edge
      if ((node1_ != graph_->edges_.size()) && (node1_ >= graph_->edges_[node1_][node2_])) {
        ++(*this);
      }

    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

 /** Creates an EdgeIterator pointing to the first (valid) edge.
  *
  */
  edge_iterator edge_begin() const
  {
    return EdgeIterator(this, 0, 0);
  }

 /** Creates an EdgeIterator pointing to the end of edges.
  *
  */
  edge_iterator edge_end() const
  {
    return EdgeIterator(this, edges_.size(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // nodes_ contains the Points corresponding to every node
  std::vector <Point> nodes_;
  // nodes_values_ contains the values corresponding to every node
  std::vector <node_value_type> nodes_values_;
  // edges_ contains the neighbours of each node
  std::vector <std::vector <size_type> > edges_;

  size_type node_size_;
  size_type edge_size_;

};


#endif // CME212_GRAPH_HPP
