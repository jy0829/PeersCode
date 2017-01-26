#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <utility>
#include <vector>

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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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

  /** Vector of Points that represent each node in graph */
  std::vector<Point> nodes;
  /** Vector containing indices of nodes that form each edge */
  std::vector<std::pair<size_type, size_type>> edges;
  /** Map containing strings of adjacent nodes */
  std::map<size_type, std::string> connections;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {

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
    Node(const graph_type* g, size_type nn) 
      : graph_(const_cast<graph_type*>(g)), node_index_(nn) {

    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes[node_index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_index_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ( (this->node_index_ == n.node_index_) && (n.graph_ == this->graph_) ) {
        return true;
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
      if (this->node_index_ < n.node_index_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  
    /** Pointer back to the Graph container */
    graph_type *graph_;
    /** This node's index in Graph container */
    size_type node_index_;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
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
    nodes.push_back(position); 
    return Node(this, nodes.size() - 1);    
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /** Using the node invariant property node(i) => node.node_index_ = i */
    size_type idx = n.node_index_;
    /** Check if graph has at least this many indices and check node equality */
    if ( (num_nodes() > idx) && (node(idx) == n) ) {
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge(const graph_type* g, size_type na, size_type nb) 
      : graph_(const_cast<graph_type*>(g)), 
        uid_a(na),
        uid_b(nb) {

    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->graph_, uid_a);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(this->graph_, uid_b);   
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_) {
        return false;
      }
      Node a = node1();
      Node b = node2();

      Node e1 = e.node1();
      Node e2 = e.node2();
      if ( (a == e1 && b == e2) || (a == e2 && b == e1) ) {
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
      if (this->uid_a < e.uid_a && this->uid_b < e.uid_b) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type *graph_;
    size_type uid_a;
    size_type uid_b;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, edges[i].first, edges[i].second);    
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // iterate through edges and check if there exists an edge
    // that connects @a a and @a b.

    /*
    for (size_type it = 0; it < edges.size(); it++) {
      if ( (edge(it).node1() == a && edge(it).node2() == b) 
        || (edge(it).node1() == b && edge(it).node2() == a) ) {
        return true;
      }
    }
    */

    /** Use Connections map for complexity: O(num_adjacent_nodes) */
    if (connections.find(a.node_index_) == connections.end() 
      || connections.find(b.node_index_) == connections.end() ) {
      return false;
    }
    std::string edges1 = connections.at(a.node_index_);
    std::string edges2 = connections.at(b.node_index_);
    /** add delimiters to strings to ensure we match
    *   the vertices we are searching for            */
    std::string a_ = "," + std::to_string(a.node_index_) + ",";
    std::string b_ = "," + std::to_string(b.node_index_) + ",";
    if ( (edges1.find(b_) != std::string::npos)
        && (edges2.find(a_) != std::string::npos) ) {
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
    if ( has_edge(a,b) ) {
      return Edge(this, a.node_index_, b.node_index_);
    }
    std::pair<size_type, size_type> ab = std::make_pair(a.node_index_, b.node_index_);
    edges.push_back(ab);
    /** Use delimiters in string for parsing */
    connections[a.node_index_] += "," + std::to_string(b.node_index_) + ",";
    connections[b.node_index_] += "," + std::to_string(a.node_index_) + ",";
    return edge(edges.size() - 1);  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges.clear();
    nodes.clear();
    connections.clear();
   }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
