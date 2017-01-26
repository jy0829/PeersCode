#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

class Graph {
private:
  struct proxy_nodes; // proxy to store Points for nodes
  struct proxy_edges; // proxy to storenode uids for edges
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
  Graph() : nodes(), edges() {
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point &position() const {

      // HW0: YOUR CODE HERE
      return graphs->nodes[uid].points;
      // return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid;
      return size_type(-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      // HW0: YOUR CODE HERE
      if (this->graphs == n.graphs && this->uid == n.uid) {
        return true;
      }
      (void)n; // Quiet compiler warning
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
      // HW0: YOUR CODE HERE
      // if both nodes belong to same graph we compare uids of the nodes.
      // if the graphs are not the same, we compare the graph pointers
      //		(mentioned in class)
      if (this->graphs == n.graphs) {
        if (this->uid < n.uid) {
          return true;
        } else {
          return false;
        }
      } else if (this->graphs > n.graphs) {
        return false;
      } else if (this->graphs < n.graphs) {
        return true;
      }
      (void)n; // Quiet compiler warning
      return false;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph *graphs;
    size_type uid;
    Node(const Graph *graphs, size_type uid)
        : graphs(const_cast<Graph *>(graphs)), uid(uid) {}
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes.size();
    return 0;
  }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position) {
    // HW0: YOUR CODE HERE
    proxy_nodes P;
    P.points = position;
    nodes.emplace_back(P);
    return Node(this, this->size() - 1);
    (void)position; // Quiet compiler warning
    return Node();  // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    if (n.graphs == this && n.uid < this->size()) {
      return true;
    }
    // HW0: YOUR CODE HERE

    (void)n; // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < this->size()) {
      return Node(this, i);
    }
    (void)i;       // Quiet compiler warning
    return Node(); // Invalid node
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graphs, graphs->edges[index].uid1);
      return Node(); // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graphs, graphs->edges[index].uid2);
      return Node(); // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      if (this->graphs == e.graphs && this->index == e.index) {
        return true;
      }
      (void)e; // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      // compared graphs as mentioned in class
      if (this->graphs == e.graphs) {
        if (this->index < e.index) {
          return true;
        } else {
          return false;
        }
      } else if (this->graphs < e.graphs) {
        return true;
      } else {
        return false;
      }
      (void)e; // Quiet compiler warning
      return false;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph *graphs;
    size_type index;
    Edge(const Graph *graphs, size_type index)
        : graphs(const_cast<Graph *>(graphs)), index(index) {}
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
    return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_edges()) {
      return Edge(this, i);
    }
    (void)i;       // Quiet compiler warning
    return Edge(); // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    // HW0: YOUR CODE HERE
    for (int i = 0; size_type(i) < num_edges(); i++) {
      if ((edges[i].uid1 == a.uid && edges[i].uid2 == b.uid) ||
          (edges[i].uid1 == b.uid && edges[i].uid2 == a.uid)) {
        return true;
      }
    }
    (void)a;
    (void)b; // Quiet compiler warning
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
    // HW0: YOUR CODE HERE
    for (int i = 0; size_type(i) < num_edges(); i++) {
      if ((edges[i].uid1 == a.uid && edges[i].uid2 == b.uid)) {
        return Edge(this, i);
      }
    }
    proxy_edges E;
    E.uid1 = a.uid;
    E.uid2 = b.uid;
    edges.emplace_back(E);
    return Edge(this, this->size() - 1);
    (void)a, (void)b; // Quiet compiler warning
    return Edge();    // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges.clear();
    nodes.clear();
    // HW0: YOUR CODE HERE
  }

private:
  struct proxy_nodes {
    Point points;
  };

  std::vector<proxy_nodes> nodes; // vector to store the nodes through the proxy

  struct proxy_edges {
    size_type uid1;
    size_type uid2;
  };
  std::vector<proxy_edges> edges; // vector to store the edges through the proxy

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
