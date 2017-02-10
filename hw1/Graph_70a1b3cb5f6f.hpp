#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"
#include <algorithm>
#include <cassert>
#include <tuple>
#include <vector>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V> class Graph {
public:
  using size_type = unsigned;

private:
  struct proxy_nodes; // proxy to store Points for nodes
  size_type nedges = 0;

public:
  using node_value_type = V;
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes(), adjacency() {}

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
    const Point &position() const { return graphs->nodes[uid].points; }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      if (uid < graphs->nodes.size())
        return uid;
      return 0;
    }

    /**   @return @a val by refernce of type V for the node.
    *@pre @a uid < number of nodes in the graph
    *@pre must be called on valid node.
    */
    node_value_type &value() { return graphs->nodes[uid].val; }
    /**  @return @a val for by refernce of type const V for the node
    *    @pre @a uid < number of nodes in the graph
    *    @pre must be called on valid node.
    */
    const node_value_type &value() const { return graphs->nodes[uid].val; }
    /** degree() returns the number of neigbours to the given node
    *@return of type size_type
    *    @pre must be called on valid node.
    */
    size_type degree() const { return graphs->adjacency[uid].size(); }
    /** edge_begin() returns an incident iterator for given node which points to
    *the first neighbour of the node
    * @return of type incident_iterator
    * @pre must be called on a valid node.
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graphs, uid, 0);
    }
    /** edge_end() returns an incident iterator for given node which points to
    *one
    *past the last neighbour of the node
    * @return of type incident_iterator
    * @pre must be called on a valid node.
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graphs, uid, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      if (this->graphs == n.graphs && this->uid == n.uid) {
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
      return false;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph *graphs;
    size_type uid;
    Node(const Graph *graphs, size_type uid)
        : graphs(const_cast<Graph *>(graphs)), uid(uid) {}
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
  Node add_node(const Point &position,
                const node_value_type &node_val = node_value_type()) {
    proxy_nodes P;
    P.points = position;
    P.val = node_val;
    nodes.emplace_back(P);
    adjacency.emplace_back(std::vector<size_type>());
    return Node(this, this->size() - 1);
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
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < this->size()) {
      return Node(this, i);
    }
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
  class Edge : private totally_ordered<Edge> {
  public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const { return graphs->node(uid1); }

    /** Return the other node of this Edge */
    Node node2() const { return graphs->node(uid2); }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      if (std::tie(graphs, std::min(uid1, uid2), std::max(uid1, uid2)) ==
          std::tie(e.graphs, std::min(e.uid1, e.uid2),
                   std::max(e.uid1, e.uid2))) {
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
      // compared graphs as mentioned in class
      if (std::tie(graphs, std::min(uid1, uid2), std::max(uid1, uid2)) <
          std::tie(e.graphs, std::min(e.uid1, e.uid2),
                   std::max(e.uid1, e.uid2))) {
        return true;
      }
      if (this->graphs < e.graphs) {
        return true;
      }
      return false;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph *graphs;
    size_type uid1;
    size_type uid2;
    Edge(const Graph *graphs, size_type uid1, size_type uid2)
        : graphs(const_cast<Graph *>(graphs)), uid1(uid1), uid2(uid2) {}
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return nedges; }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i < num_edges()) {
      return *std::next(edge_begin(), i);
    }
    return Edge(); // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    for (int i = 0; size_type(i) < adjacency[a.uid].size(); i++) {
      if (b.uid == adjacency[a.uid][i]) {
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
    if (has_edge(a, b)) {
      return Edge(this, a.uid, b.uid);
    }
    adjacency[a.uid].push_back(b.uid);
    nedges = nedges + 1;
    adjacency[b.uid].push_back(a.uid);
    return Edge(this, a.uid, b.uid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    adjacency.clear();
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

    /** operator*() returns the node the node pointer is pointing to.
    * @return of type Node.
    * @pre !(NodeIterator == graph.node_end())
    */
    Node operator*() const { return graphs->node(curr); }
    /** operator++() increments the iterator by one position.
    * @return a pointer to a node of type @a NodeIterator
    * @pre !(NodeIterator == graph.node_end())*/
    NodeIterator &operator++() {
      this->curr = curr + 1;
      return *this;
    }
    /** operator==() compares the current @NodeIterator and given @NodeIterator.
    * @param[in] @a n is the input @a NodeIterator for comparison
    *@return type bool
     */
    bool operator==(const NodeIterator &n) const {
      if (this->graphs == n.graphs && n.curr == this->curr) {
        return true;
      }
      return false;
    }
    /** operator!=() is true if !(operator==())
    * @param[in] @a n is the input @a NodeIterator for comparison
    *@return type bool
     */

    bool operator!=(const NodeIterator &n) const { return !(*this == n); }

  private:
    friend class Graph;
    Graph *graphs;
    size_type curr;

    NodeIterator(const Graph *graphs, size_type p)
        : graphs(const_cast<Graph *>(graphs)), curr(p) {}
  };

  /** node_begin() returns a pointer to the first node in graph
*@return type node_iterator
 */
  node_iterator node_begin() const { return NodeIterator(this, 0); }

  /** node_end() returns a pointer to one past the last node in graph
  *@return type node_iterator
   */
  node_iterator node_end() const { return NodeIterator(this, size()); }

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

    /** operator*() returns the edge the incident iterator pointer is pointing
    * to.
    * @return of type Edge.
    * @pre !(IncidentIterator == node.edge_end())
    */
    Edge operator*() const {
      return Edge(graphs, node_a, graphs->adjacency[node_a][incident_curr]);
    }
    /** operator++() increments the iterator by one position.
    * @return a pointer to a edge of type @a IncidentIterator
    * @pre !(IncidentIterator == node.node_end())*/
    IncidentIterator &operator++() {
      incident_curr = incident_curr + 1;
      return *this;
      // return IncidentIterator(graphs,incident_curr+1,node_a);
    }
    /** operator==() compares the current @IncidentIterator and given
    *@IncidentIterator.
    * @param[in] @a Inc is the input @a IncidentIterator for comparison
    *@return type bool
     */
    bool operator==(const IncidentIterator &Inc) const {
      if (this->graphs == Inc.graphs && this->node_a == Inc.node_a &&
          this->incident_curr == Inc.incident_curr) {
        return true;
      }
      return false;
    }
    /** operator!=() is true if !(operator==())
    * @param[in] @a Inc is the input @a IncidentIterator for comparison
    *@return type bool
     */
    bool operator!=(const IncidentIterator &Inc) const {
      return !(*this == Inc);
    }

  private:
    friend class Graph;
    Graph *graphs;
    size_type incident_curr;
    size_type node_a;

    IncidentIterator(const Graph *graphs, size_type node_a, size_type curr_p)
        : graphs(const_cast<Graph *>(graphs)), incident_curr(curr_p),
          node_a(node_a) {}
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
    /** operator*() returns the edge the edge iterator pointer is pointing to.
    * @return of type Edge.
    * @pre !(EdgeIterator == graph.edge_end())
    */
    Edge operator*() const {
      return Edge(graphs, (*node_curr).index(), (*inc_curr).node2().index());
    }
    /** operator++() increments the iterator by one position.
    * @return a pointer to a edge of type @a EdgeIterator
    * @pre !(EdgeIterator == graph.edge_end())*/
    EdgeIterator &operator++() {
      ++inc_curr;
      fix();
      return *this;
    }
    /** operator==() compares the current @EdgeIterator and given @EdgeIterator.
    * @param[in] @a e is the input @a EdgeIterator for comparison
    *@return type bool
     */
    bool operator==(const EdgeIterator &e) const {
      if (graphs == e.graphs && node_curr == e.node_curr &&
          inc_curr == e.inc_curr) {
        return true;
      }
      return false;
    }
    /** operator!=() is true if !(operator==())
    * @param[in] @a e is the input @a EdgeIterator for comparison
    *@return type bool
     */
    bool operator!=(const EdgeIterator &e) const { return !(*this == e); }

  private:
    friend class Graph;
    Graph *graphs;
    NodeIterator node_curr = (*graphs).node_begin();
    IncidentIterator inc_curr = (*node_curr).edge_begin();

    EdgeIterator(const Graph *graphs, NodeIterator p1, IncidentIterator p2)
        : graphs(const_cast<Graph *>(graphs)), node_curr(p1), inc_curr(p2) {
      fix();
    }

    void fix() {
      while (node_curr != graphs->node_end()) {
        while (inc_curr != (*node_curr).edge_end()) {
          if ((*node_curr) < (*inc_curr).node2()) {
            ++inc_curr;
          } else {
            break;
          }
        }
        if (inc_curr != (*node_curr).edge_end()) {
          break;
        } else {
          ++node_curr;
          inc_curr = (*node_curr).edge_begin();
        }
      }
      if (node_curr == graphs->node_end()) {
        node_curr = (*graphs).node_begin();
        inc_curr = (*node_curr).edge_begin();
      }
    }
  };

  /** edge_begin() returns a pointer to the first edge in graph
  *@return type edge_iterator
   */

  edge_iterator edge_begin() const {
    return EdgeIterator(this, node_begin(), (*node_begin()).edge_begin());
    // return EdgeIterator(this,0,0);
  }
  /** edge_end() returns a pointer to one past the last edge in graph
  *@return type edge_iterator
   */
  edge_iterator edge_end() const {

    return EdgeIterator(this, node_end(), (*node_end()).edge_end());
  }

private:
  struct proxy_nodes {
    Point points;
    node_value_type val;
  };

  std::vector<proxy_nodes> nodes; // vector to store the nodes through the proxy
  std::vector<std::vector<size_type>> adjacency;
};

#endif // CME212_GRAPH_HPP
