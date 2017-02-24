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
template <typename V, typename E> class Graph {
public:
  using size_type = unsigned;

private:
  struct proxy_nodes; // proxy to store Points for nodes
  size_type nedges = 0;

public:
  using node_value_type = V;
  using edge_value_type = E;
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
  Graph() : nodes() {}

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

    Point &position() { return graphs->nodes[uid].points; }
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      if (uid < graphs->nodes.size())
        return graphs->nodes[uid].idx;
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
    size_type degree() const { return graphs->nodes[uid].adjacency.size(); }
    /** edge_begin() returns an incident iterator for given node which points to
    *the first neighbour of the node
    * @return of type incident_iterator
    * @pre must be called on a valid node.
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graphs, graphs->nodes[uid].idx, 0);
    }
    /** edge_end() returns an incident iterator for given node which points to
    *one
    *past the last neighbour of the node
    * @return of type incident_iterator
    * @pre must be called on a valid node.
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graphs, graphs->nodes[uid].idx, degree());
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
  size_type size() const { return i2u.size(); }

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
    P.adjacency = std::vector<proxy_edges>();
    P.idx = i2u.size();
    nodes.emplace_back(P);
    i2u.push_back(nodes.size() - 1);

    return Node(this, nodes.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    if (n.graphs == this && n.uid < nodes.size()) {
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
  Node node(size_type i) const { return Node(this, i2u[i]); }

  /** Removes specified node and incident edges from the graph.
  * @param[in,out] @a n is the node to be removed.
  * @return 1 if node is deleted and 0 if it is not.
  *
  * @post All elements with position i < @a n.index(), maintain same position.
  * if returns 1, All elements with position i > @a n.index() are moved to
  * position i-1.
  * @post if returns 1, new num_nodes() = old num_nodes() - 1
  * @post if returns 1,  neigbours to node n, has new degree()= old degree() -1
  * @note Invalidates the node.
  *  At most O(num_nodes()) operations for sparse graphs.
  */
  size_type remove_node(const Node &n) {

    if (!has_node(n)) {
      return 0;
    }
    for (size_type i = 0; i < i2u.size(); i++) {
      if (nodes[i2u[i]].idx > n.index()) {
        --nodes[i2u[i]].idx;
      }
    }
    for (size_type j = 0; j < nodes[n.uid].adjacency.size(); ++j) {
      auto neigh = nodes[n.uid].adjacency[j].node2_id;
      for (size_type i = 0; i < nodes[neigh].adjacency.size(); ++i) {
        if (nodes[neigh].adjacency[i].node2_id == n.uid) {
          nodes[neigh].adjacency[i] =
              nodes[neigh].adjacency[nodes[neigh].adjacency.size() - 1];
          nodes[neigh].adjacency.pop_back();
          --nedges;
          break;
        }
      }
    }

    for (auto k = i2u.begin(); k != i2u.end(); ++k) {
      if (n.uid == (*k)) {
        i2u.erase(k);
        break;
      }
    }
    nodes[n.uid].idx = -1;
    return 1;
  }
  /** Removes specified node and incident edges from the graph.
  * @param[in] n_it Iterator pointing to the node to be removed.
  * @return @a n_it iterator if node is deleted and @a n_it iterator if it is
  * not.
  * @pre old node_begin() <= @a n_it < old node_end().
  * @post All elements with position i < @a (*n_it).index(), maintain same
  * position.
  * if returns 1, All elements with position i > @a (*n_it).index() are moved to
  * position i-1.
  * @post if returns 1, new num_nodes() = old num_nodes() - 1
  * @post if returns 1, neigbours to node (*n_it), has new degree()= old
  * degree() -1
  * @note Invalidates all iterators.
  *  At most O(num_nodes()) operations for sparse graphs.
  */
  node_iterator remove_node(node_iterator n_it) {
    auto n = *n_it;
    remove_node(n);
    return n_it;
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
    edge_value_type &value() {
      for (size_type i = 0; i < graphs->nodes[node1().uid].adjacency.size();
           i++) {
        if (node2().uid == graphs->nodes[node1().uid].adjacency[i].node2_id) {
          return graphs->nodes[node1().uid].adjacency[i].e_val;
        }
      }
      return graphs->nodes[node1().uid].adjacency[0].e_val;
    }

    const edge_value_type &value() const {
      for (size_type i = 0; i < graphs->nodes[uid1].adjacency.size(); i++) {
        if (uid2 == graphs->nodes[uid1].adjacency[i].node2_id) {
          return graphs->nodes[uid1].adjacency[i].e_val;
        }
      }
      return graphs->nodes[uid1].adjacency[0].e_val;
    }
    double length() const {
      return norm(node1().position() - node2().position());
    }
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
        : graphs(const_cast<Graph *>(graphs)), uid1(graphs->nodes[uid1].idx),
          uid2(graphs->nodes[uid2].idx) {}
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
    for (int i = 0; size_type(i) < nodes[a.uid].adjacency.size(); i++) {
      if (b.uid == nodes[a.uid].adjacency[i].node2_id) {
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
    proxy_edges e;
    e.node2_id = b.uid;
    e.e_val = edge_value_type();
    nodes[a.uid].adjacency.push_back(e);
    nedges = nedges + 1;
    e.node2_id = a.uid;
    nodes[b.uid].adjacency.push_back(e);
    return Edge(this, a.uid, b.uid);
  }

  /** Removes specified edge from the graph.
  * @param[in,out] @a n1 is the node incident to the edge.
  * @param[in,out] @a n2 is the node incident to the edge.
  * @return 1 if edge is deleted and 0 if it is not.
  * @post if returns 1, new num_edges() = old num_edges() - 1
  * @post if returns 1,  node @a n1 and node @a n2, has new degree()= old
  * degree() -1
  *  At most O(num_nodes()) operations for sparse graphs.
  */
  size_type remove_edge(const Node &n1, const Node &n2) {
    if (!has_edge(n1, n2)) {
      return 0;
    }
    for (size_type j = 0; j < nodes[n1.uid].adjacency.size(); ++j) {
      if (nodes[n1.uid].adjacency[j].node2_id == n2.uid) {
        nodes[n1.uid].adjacency[j] =
            nodes[n1.uid].adjacency[nodes[n1.uid].adjacency.size() - 1];
        nodes[n1.uid].adjacency.pop_back();
        --nedges;
        break;
      }
    }
    for (size_type j = 0; j < nodes[n2.uid].adjacency.size(); ++j) {
      if (nodes[n2.uid].adjacency[j].node2_id == n1.uid) {
        nodes[n2.uid].adjacency[j] =
            nodes[n2.uid].adjacency[nodes[n2.uid].adjacency.size() - 1];
        nodes[n2.uid].adjacency.pop_back();
        break;
      }
    }
    return 1;
  }
  /** Removes specified edge from the graph.
  * @param[in,out] @a e is the edge incident to be removed.
  * @return 1 if edge is deleted and 0 if it is not.
  * @post if returns 1, new num_edges() = old num_edges() - 1
  * @post if returns 1,  nodes incident to e node @a n1 and node @a n2, has new
  * degree()= old degree() -1
  *  At most O(num_nodes()) operations for sparse graphs.
  */
  size_type remove_edge(const Edge &e) {
    auto n1 = e.node1();
    auto n2 = e.node2();
    return remove_edge(n1, n2);
  }
  /** Removes specified edge from the graph.
  * @param[in] e_it Iterator pointing to the edge to be removed.
  * @return @a e_it iterator if edge is deleted and @a e_it iterator if it is
  * not.
  * @pre old edge_begin() <= @a e_it < old edge_end().
  * @post if returns 1, new num_edges() = old num_edges() - 1
  * @post if returns 1, nodes incident  to edge (*e_it), has new degree()= old
  * degree() -1
  * @note Invalidates all edge and incident iterators.
  *  At most O(num_nodes()) operations for sparse graphs.
  */

  edge_iterator remove_edge(edge_iterator e_it) {
    auto e = *e_it;
    remove_edge(e);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    i2u.clear();
    nedges = 0;
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

      return Edge(
          graphs, graphs->i2u[node_a],
          graphs->nodes[graphs->i2u[node_a]].adjacency[incident_curr].node2_id);
    }
    /** operator++() increments the iterator by one position.
    * @return a pointer to a edge of type @a IncidentIterator
    * @pre !(IncidentIterator == node.node_end())*/
    IncidentIterator &operator++() {

      incident_curr = incident_curr + 1;

      return *this;
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
    size_type node_a;
    size_type incident_curr;

    IncidentIterator(const Graph *graphs, size_type node_a, size_type curr_p)
        : graphs(const_cast<Graph *>(graphs)), node_a(node_a),
          incident_curr(curr_p) {}
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
      // return Edge(graphs, (*node_curr).index(), (*inc_curr).node2().index());
      return (*inc_curr);
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

      if ((graphs == e.graphs) && (node_curr == e.node_curr) &&
          (inc_curr == e.inc_curr)) {
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
      int flag = 0;
      while (node_curr != graphs->node_end()) {
        while (inc_curr != (*node_curr).edge_end()) {
          if ((*node_curr) < (*inc_curr).node2()) {

            flag = 1;
            break;
          }
          ++inc_curr;
        }
        if (flag == 1) {
          break;
        }
        ++node_curr;

        if (node_curr == graphs->node_end()) {
          inc_curr = (*(graphs->node_begin())).edge_begin();
          break;
        }

        inc_curr = (*(node_curr)).edge_begin();
      }
    }
  };

  /** edge_begin() returns a pointer to the first edge in graph
  *@return type edge_iterator
   */

  edge_iterator edge_begin() const {
    auto n = *node_begin();
    return EdgeIterator(this, node_begin(), (n).edge_begin());
  }
  /** edge_end() returns a pointer to one past the last edge in graph
  *@return type edge_iterator
   */
  edge_iterator edge_end() const {
    auto n = *node_begin();
    return EdgeIterator(this, node_end(), (n).edge_begin());
  }

private:
  struct proxy_edges {
    size_type node2_id;
    edge_value_type e_val;
  };
  struct proxy_nodes {
    Point points;
    node_value_type val;
    std::vector<proxy_edges> adjacency;
    size_type idx;
  };

  std::vector<proxy_nodes> nodes; // vector to store the nodes through the proxy
  std::vector<size_type> i2u;
};

#endif // CME212_GRAPH_HPP
