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
class Graph {
 private:
  
  // Internal types

  struct edge_wrapper
  {
    unsigned nid1;
    unsigned nid2;
  };

  std::vector<Point> points;
  std::vector<edge_wrapper> edges;

  unsigned no_nodes;
  unsigned no_edges;

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

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    : points(), edges(), no_nodes(0), no_edges(0) {
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
  class Node {
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
      
      return get();

    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      // ensure index is in the correct range
      if (uid_ < graph_->no_nodes)
        return uid_;

      // this code should not be executed
      assert(false);
      return size_type(-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
    
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
    bool operator<(const Node& n) const {
      
      // check if the node belongs to the same graph and if uid is less
      if (this->graph_ == n.graph_ && this->uid_ < n.uid_)
        return true;
      
      if (this->graph_ < n.graph_) 
        return true;
      
      return false;

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    // Pointer to original graph 
    graph_type* graph_;

    // The index of the node
    size_type uid_;

    // Private constructor
    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

    // Define get function
    Point& get() const {

      return graph_->points[uid_];

    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {

    return no_nodes;

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
  Node add_node(const Point& position) {

    // add element to vector
    points.push_back(position);

    no_nodes += 1;

    // return node
    return Node(this, no_nodes - 1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    
    // check if node is in the graph
    if (n.uid_ < no_nodes && this == n.graph_) {
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
      return Node(this, i);

    // if i is bigger than the size then return invalid node

    return Node();        // Invalid node
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      
    }

    /** Return a node of this Edge */
    Node node1() const {
      
      // call the node function. the order takes care
      // of the issue raised on canvas
      if (order_) {
        return graph_->node(graph_->edges[euid_].nid1);
      }
      else {
        return graph_->node(graph_->edges[euid_].nid2);
      }
      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      
      // very similar to the above implementation
      if (order_) {
        return graph_->node(graph_->edges[euid_].nid2);
      }
      else {
        return graph_->node(graph_->edges[euid_].nid1);
      }
      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      // check if two edges belong to the same graph and have the same euid
      if(this->graph_ == e.graph_ && this->euid_ == e.euid_) {
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

      // a global comparison between two edges
      if(this->graph_ == e.graph_ && this->euid_ < e.euid_) 
        return true;
      
      if(this->graph_ < e.graph_)
        return true;
      
      return false;

    }

   private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // pointer to the graph
    graph_type *graph_;

    // edge id
    size_type euid_;

    // order in which nodes are stored. by default order is true, 
    // meaning the struct always contains the smaller node id first.
    bool order_;

    // Private constructor
    Edge(const graph_type* graph, size_type euid, bool order)
        : graph_(const_cast<graph_type*>(graph)), euid_(euid), order_(order) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    
    return no_edges;

  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {

    // check that i is less than vector size
    if (i < no_edges) {

      return Edge(this, i, true);

    }

    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    
    // iterate through the edges vector
    for (size_type i = 0; i < no_edges; i++) {

      edge_wrapper curr_edge = edges[i];

      // check if the edge pair exists
      if ((a.uid_ == curr_edge.nid1 && b.uid_ == curr_edge.nid2) || 
          (a.uid_ == curr_edge.nid2 && b.uid_ == curr_edge.nid1))
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {

    // check if edge exists
    bool present = false;
    size_type index = 0;

    // iterate through the vector
    for (size_type i = 0; i < no_edges; i++) {

      edge_wrapper curr_edge = edges[i];

      // compare current edge nodes to a and b
      if ((a.uid_ == curr_edge.nid1 && b.uid_ == curr_edge.nid2) || 
          (a.uid_ == curr_edge.nid2 && b.uid_ == curr_edge.nid1)) {
        present = true;
        index = i;
        break;
      }

    }    

    // if edge is present return a proxy with order_ initialized
    // to fulfil the return conditions
    if (present) {
      return Edge(this, index, a.uid_ <= b.uid_);
    }
    else {
      // create new edge
      // by default we store nodes with lower uids first
      edge_wrapper new_edge;
      if (a.uid_ <= b.uid_)
        new_edge = {.nid1 = a.uid_, .nid2 = b.uid_};
      else
        new_edge = {.nid1 = b.uid_, .nid2 = a.uid_};

      // add edge to vector and return it
      edges.push_back(new_edge);
      no_edges += 1;

      return Edge(this, no_edges - 1, true);
    }

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    
    points.clear();
    edges.clear();
    no_edges = 0;
    no_nodes = 0;

  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
