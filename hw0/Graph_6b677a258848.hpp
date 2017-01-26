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
  Graph() : nodes(), edges() {}

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
      // do nothing- invalid node
    }

    /** Return this node's position. */
    const Point &position() const {

      return (graph->nodes[uid]);
      // this might be invalid in the case the graph pointer is not pointing to
      // a valid graph
      // but the specifications don't seem clear on how to deal with it.
      // exception handling will increase run time etc.
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      if (uid < graph->nodes.size())
        return uid;
      // there is no clear way to return an invalid node's size, like -1
      // since the return type is unsigned.
      return 0;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {

      if ((this->graph == n.graph) and (this->uid == n.uid)) {
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
    bool operator<(const Node &n) const {
      // if graphs are the same, check for node ids otherwise compare graph
      // pointers
      if ((this->graph == n.graph) and (this->uid < n.uid)) {
        return true;
      } else if (this->graph < n.graph) {
        return true;
      }
      return false;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // graph pointer tracks the Graph to which the node belongs
    graph_type *graph;
    // unique id for a Node
    size_type uid;
    // Valid private constructor for a Node
    Node(const graph_type *g, size_type id)
        : graph(const_cast<graph_type *>(g)), uid(id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { return nodes.size(); }

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

    size_type id = nodes.size();
    // insert node at the end of the vector
    nodes.push_back(position);
    // return node_type
    node_type node = node_type(this, id);
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    // node belongs to a graph, if its uid is in the correct range
    // and it has a pointer to this graph
    if ((n.graph == this) and (n.uid < nodes.size()))
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { return Node(this, i); }

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
      // do nothing - invalid edge
    }

    /** Return a node of this Edge */
    Node node1() const {
      // this might be invalid in the case the graph pointer is not pointing to
      // a valid graph
      // but the specifications don't seem clear on how to deal with it.
      // exception handling will increase run time etc.
      return graph->node(graph->edges[uid].u);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // this might be invalid in the case the graph pointer is not pointing to
      // a valid graph
      // but the specifications don't seem clear on how to deal with it.
      // exception handling will increase run time etc.
      return graph->node(graph->edges[uid].v);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      // edges with same uid and same graph are equal
      if ((this->uid == e.uid) && (this->graph == e.graph))
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      // if the graphsa are different compare uids else compare the graph
      // pointers
      if ((this->uid < e.uid) && (this->graph == e.graph)) {
        return true;
      } else if (this->graph < e.graph) {
        return true;
      }
      return false;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // graph pointer to the Graph to which the Edge belongs
    graph_type *graph;
    // unique id for the Edge
    size_type uid;
    // Private constructor for valid edge objects
    Edge(const graph_type *g, size_type id)
        : graph(const_cast<graph_type *>(g)), uid(id) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return edges.size(); }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { return Edge(this, i); }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    for (auto e : edges) {
      if ((e.u == a.uid and e.v == b.uid) || (e.v == a.uid and e.u == b.uid))
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
  Edge add_edge(const Node &a, const Node &b) {
    // check if the edge exists, if not add it
    for (size_type i = 0; i < edges.size(); i++) {
      if ((edges[i].u == a.uid and edges[i].v == b.uid) ||
          (edges[i].v == a.uid and edges[i].u == b.uid)) {
        return edge(i);
      }
    }

    size_type id = edges.size();
    Edge edge = edge_type(this, id);
    internal_edge e;
    e.u = a.uid;
    e.v = b.uid;
    edges.emplace_back(e);
    return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    edges.clear();
    nodes.clear();
  }

private:
  // struct used for the proxy for edges
  struct internal_edge {
    // u,v denote the unique ids for the two nodes for the edge
    size_type u, v;
  };
  // vector of nodes keeps track of the nodes belonging to a Graph
  std::vector<Point> nodes;
  // vector of edges keeps track the edges belonging to a Graph
  std::vector<internal_edge> edges;
};

#endif // CME212_GRAPH_HPP
