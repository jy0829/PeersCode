#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
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

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Allow user-specified node value. */
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
    // HW0: YOUR CODE HERE
    graph_node_index_ = 0;
    graph_edge_index_ = 0;
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
      graph_ = NULL;
      node_index_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> nodes[node_index_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return this node's value*/
    node_value_type& value(){
      assert(graph_);
      return graph_ -> nodes[node_index_].value_;
    }

    /** Return this const node's const value*/
    const node_value_type& value() const{
      assert(graph_);
      return graph_ -> nodes[node_index_].value_;     
    }

    /** Return this node's degree*/
    size_type degree() const {
      assert(graph_);
      return graph_ -> adjacents[node_index_].size();
    }

    /** Return this node's begin incident iterator*/
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, node_index_, 0);
    }

    /** Return this node's end incident iterator*/
    incident_iterator edge_end() const {
      return incident_iterator(graph_, node_index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // return true when the nodes have the same graph and the same index;
      if (graph_ && n.graph_) {
        return graph_ == n.graph_ && node_index_ == n.node_index_;
      }
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
    bool operator<(const Node& n) const {
      return std::tie(graph_, node_index_) < std::tie(n.graph_, n.node_index_);

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type node_index_;
    Node(const graph_type* graph, size_type node_index) {
      graph_ = const_cast<graph_type*>(graph);
      node_index_ = node_index;
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes.size();
    // return 0;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {

    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[] node_value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    
    nodes.push_back(node_info(position, node_value));
    adjacents.push_back(std::vector<size_type>());
    return node(nodes.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (!n.graph_) return true;
    return (this == n.graph_ && n.node_index_ < nodes.size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i);        // Invalid node
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
      graph_ = NULL;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // return true when graph_, node1_, and node2_ are the same;
      if (graph_ && e.graph_) {
        return graph_ == e.graph_ && node1_ == e.node1_ && node2_ == e.node2_;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // same as node
      return std::tie(graph_, node1_, node2_) < std::tie(e.graph_, e.node1_, e.node2_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type node1_;
    size_type node2_;
    Edge(const graph_type* gragh, size_type node1, size_type node2) {
      graph_ = const_cast<graph_type*>(gragh);
      node1_ = node1;
      node2_ = node2;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return graph_edge_index_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  // Edge edge(size_type i) const __attribute__((deprecated)) {
  //   // HW0: YOUR CODE HERE
  //   return *std::next(edge_begin(), i);
  // }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (graph_edge_index_ == 0) return false;
    for (size_type i = 0; i < adjacents[a.node_index_].size(); ++i) {
      if (adjacents[a.node_index_][i] == b.node_index_) return true;
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
    if (has_edge(a, b)) return edge_type(this, a.node_index_, b.node_index_);
    adjacents[a.node_index_].push_back(b.node_index_);
    adjacents[b.node_index_].push_back(a.node_index_);
    graph_edge_index_++;
    return edge_type(this, a.node_index_, b.node_index_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear();
    adjacents.clear();
    graph_edge_index_ = 0;
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
    /** Return the node that the iterator points to. */
    Node operator*() const {
      return value_type(graph_, node_index_);
    }

    /** Point to next node and return the node. */
    NodeIterator& operator++() {
      node_index_++;
      return *this;
    }
    /** Return whether two nodes are same or not: graph and index should be the same. */
    bool operator==(const NodeIterator& nit) const {
      if (graph_ && nit.graph_) {
        return graph_ == nit.graph_ && node_index_ == nit.node_index_;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_index_;
    NodeIterator (const Graph* graph, size_type node_index):graph_(const_cast<Graph*>(graph)), node_index_(node_index){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return the begin iterator of nodes*/
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /** Return the end iterator of nodes*/
  node_iterator node_end() const {
    return NodeIterator(this, nodes.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
    /** Return the edge that the incidentIterator points to. */
    Edge operator*() const {
      return edge_type(graph_, node_index_, graph_ -> adjacents[node_index_][adjacents_index_]);
    }

    /** Return the next incidentIterator. */
    IncidentIterator& operator++() {
      ++adjacents_index_;
      return *this;
    }
    /** Return if two interators are point to the same thing or not*/
    bool operator==(const IncidentIterator& iit) const{
      if (graph_ && iit.graph_) {
        return graph_ == iit.graph_ && node_index_ == iit.node_index_ && adjacents_index_ == iit.adjacents_index_;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_index_;
    size_type adjacents_index_;
    IncidentIterator(const Graph* graph, size_type node_index, size_type adjacents_index){
       graph_ = const_cast<Graph*>(graph);
       node_index_ = node_index;
       adjacents_index_ = adjacents_index;
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      graph_ = NULL;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the edge that the iterator points to*/
    Edge operator*() const {
      assert(node_index_ < graph_ -> num_nodes() && adjacents_index_ < graph_ -> adjacents[node_index_].size());
      return edge_type(graph_, node_index_, graph_ -> adjacents[node_index_][adjacents_index_]);
    }

    /** Increase iterator and return the new iterator*/
    EdgeIterator& operator++() {
      ++adjacents_index_;
      while (node_index_ < graph_ -> num_nodes()) {
        while (adjacents_index_ < graph_ -> adjacents[node_index_].size()) {
          if (node_index_ < graph_ -> adjacents[node_index_][adjacents_index_])
            return *this;
          ++adjacents_index_;
        }
        ++node_index_;
        adjacents_index_ = 0;
      }
      return *this;
    }

    /** Return if two edge iterator point to the same thing or not*/
    bool operator==(const EdgeIterator& eit) const {
      if (graph_ && eit.graph_) {
        return graph_ == eit.graph_ && node_index_ == eit.node_index_ && adjacents_index_ == eit.adjacents_index_;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type node_index_;
    size_type adjacents_index_; 

    EdgeIterator(const graph_type* graph, size_type node_index, size_type adjacents_index) {
      graph_ = const_cast<graph_type*> (graph);
      node_index_ = node_index;
      adjacents_index_ = adjacents_index;
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return the begin iterator of edges*/
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }
  /** Return the end iterator of edges*/
  edge_iterator edge_end() const {
    return EdgeIterator(this, nodes.size(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct node_info {
    Point position_;
    node_value_type value_;
    node_info(const Point& position, const node_value_type& value) {
      position_ = position;
      value_ = value;
    }
  };
  std::vector<node_info> nodes;
  std::vector<std::vector<size_type>> adjacents;
  size_type graph_node_index_;
  size_type graph_edge_index_;

};

#endif // CME212_GRAPH_HPP
