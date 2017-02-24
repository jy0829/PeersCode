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
private:
 

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

  using node_value_type = V;
  using edge_value_type = E;

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
    Node() {
      // do nothing- invalid node
    }

    /** Return this node's position. */
    const Point &position() const {
      
      return (graph->nodes[uid].point);
      // this might be invalid in the case the graph pointer is not pointing to
      // a valid graph
      // but the specifications don't seem clear on how to deal with it.
      // exception handling will increase run time etc.
    }
    
    /** Return this node's position. */
    Point &position() {
      return (graph->nodes[uid].point);
      // this might be invalid in the case the graph pointer is not pointing to
      // a valid graph
      // but the specifications don't seem clear on how to deal with it.
      // exception handling will increase run time etc.
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      if (uid < graph->nodes.size())
        return graph->nodes[uid].index;
      // there is no clear way to return an invalid node's size, like -1
      // since the return type is unsigned.
      return 0;
    }

    /** returns a reference to the value of node, of the type V
*@return a reference to the value of node, of the type V
*@pre The Node should be a valid Node object
*/
    node_value_type &value() { return graph->nodes[uid].val; }
    /** Returns a reference to the value of a const node, of the type V
*@return a reference to the value of const node, of the type V
*@pre The Node should be a valid Node object
*/
    const node_value_type &value() const {
      return ((const node_value_type &)&(graph->nodes[uid].val));
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
    /**Returns the number of neighbors of the given node
*@return type size_type
*@pre The Node should be a valid Node object
*/
    size_type degree() const { return graph->nodes[uid].edges.size(); }
    /**Returns an incident iterator, pointing to the first neighbor of the node
*@return an incident iterator, pointing to the first neighbor of the node
*@pre The Node should be a valid Node object
*/
    incident_iterator edge_begin() const {
      return incident_iterator(graph, graph->nodes[uid].index, 0);
    }
    /**Returns an incident iterator, pointing to one past the last neighbor of
*the node
*@return an incident iterator, pointing to the last neighbor of the node
*@pre The Node should be a valid Node object
*/
    incident_iterator edge_end() const {
      return incident_iterator(graph, graph->nodes[uid].index, degree());
    }

    bool valid() const {
      // std::cout<<uid<<" "<<graph->nodes.size()<<"\n";
      return (uid >= 0 && uid < graph->nodes.size() &&
              graph->nodes[uid].index < (long)graph->i2u.size() &&
              graph->i2u[graph->nodes[uid].index] == uid);
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
                const node_value_type &n_value = node_value_type()) {

    size_type id = nodes.size();
    // insert node at the end of the vector
    internal_node n;
    n.point = position;
    n.val = n_value;
    n.edges = std::vector<internal_edge>();
    n.index = i2u.size();
    nodes.push_back(n);
    // populating i2u
    i2u.push_back(id);
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
    if ((n.graph == this) and (n.uid < nodes.size()) and n.valid())
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { return Node(this, i2u[i]); }

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
      // do nothing - invalid edge
    }

    /** Return a node of this Edge */
    Node node1() const {
      // this might be invalid in the case the graph pointer is not pointing to
      // a valid graph
      // but the specifications don't seem clear on how to deal with it.
      // exception handling will increase run time etc.
      // std::cout<<"here in node1()";
      return graph->node(node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // this might be invalid in the case the graph pointer is not pointing to
      // a valid graph
      // but the specifications don't seem clear on how to deal with it.
      // exception handling will increase run time etc.
      return graph->node(node2_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      // edges with same pair of node uids and same graph are equal

      return std::tie(this->graph, std::min(this->node1_id, this->node2_id),
                      std::max(this->node1_id, this->node2_id)) ==
             std::tie(e.graph, std::min(e.node1_id, e.node2_id),
                      std::max(e.node1_id, e.node2_id));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      // if the graphsa are different compare node uids else compare the graph
      // pointers
      return std::tie(this->graph, std::min(this->node1_id, this->node2_id),
                      std::max(this->node1_id, this->node2_id)) <
             std::tie(e.graph, std::min(e.node1_id, e.node2_id),
                      std::max(e.node1_id, e.node2_id));
    }
    /** Calculates the length of the edge */
    double length() const {
      Point p1 = graph->node(node1_id).position();
      Point p2 = graph->node(node2_id).position();
      Point p = p1 - p2;
      return norm(p);
    }
    
    /** Returns a reference to the edge value */
    edge_value_type &value() {
      for (unsigned int i = 0;
           i < graph->nodes[graph->i2u[node1_id]].edges.size(); i++) {
        if (graph->nodes[graph->i2u[node1_id]].edges[i].node2_id ==
            graph->i2u[node2_id]) {
          return graph->nodes[graph->i2u[node1_id]].edges[i].val;
        }
      }
      return graph->nodes[0].edges[0].val;
    }

    /** Returns a const reference to the edge value */
    const edge_value_type &value() const {
      for (unsigned int i = 0;
           i < graph->nodes[graph->i2u[node1_id]].edges.size(); i++) {
        if (graph->nodes[graph->i2u[node1_id]].edges[i].index ==
            graph->i2u[node2_id]) {
          return graph->nodes[graph->i2u[node1_id]].edges[i].val;
        }
      }
      return graph->nodes[0].edges[0].val;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // graph pointer to the Graph to which the Edge belongs
    graph_type *graph;
    // node indexes for the nodes which the edge is defined on
    size_type node1_id;
    size_type node2_id;
    // Private constructor for valid edge objects
    Edge(const graph_type *g, size_type node1id, size_type node2id)
        : graph(const_cast<graph_type *>(g)) {
      node1_id = graph->nodes[node1id].index;
      node2_id = graph->nodes[node2id].index;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return number_of_edges; }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i < num_edges())
      return *std::next(edge_begin(), i);
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {

    for (auto e : nodes[a.uid].edges) {
      if (e.node2_id == b.uid)
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

    if (has_edge(a, b) == true) {
      return edge_type(this, a.uid, b.uid);
    }
    internal_edge e1;
    e1.node2_id = b.uid;
    e1.val = edge_value_type();
    internal_edge e2;
    e2.node2_id = a.uid;
    e2.val = edge_value_type();
    nodes[a.uid].edges.push_back(e1);
    nodes[b.uid].edges.push_back(e2);
    number_of_edges++;

    return edge_type(this, a.uid, b.uid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    number_of_edges = 0;
    nodes.clear();
    i2u.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** If the iterator is valid, this returns the node to which the iterator is
*pointing to
*@return node to which the iterator is pointing to
*@pre iterator should not be equal to the graph.node_end() and should be valid
*/
    Node operator*() const { return g->node(p); }

    /**For a valid iterator, Increments the interator to point to the next node
*in the graph
   *@return the reference to an incremented node iterator
*@pre iterator should not be equal to the graph.node_end() and should be valid
*/
    NodeIterator &operator++() {
      ++p;
      return *this;
    }
    /**Compares the given node iterator to the this iterator. Returns true if
*both are equal, else returns false
*@param[in] @n , of NodeIterator type.
*@return a boolean, true if the iterators are equal, false otherwise
*/
    bool operator==(const NodeIterator &n) const { return p == n.p; }

  private:
    friend class Graph;
    // p is an index keeping track of which node the iterator is pointing to.
    size_type p;
    // pointer to the graph from which the node iterator has been called.
    graph_type *g;
    NodeIterator(const graph_type *g1, size_type p1) {

      g = const_cast<graph_type *>(g1);
      p = p1;
    }
  };

  // Supply definitions AND SPECIFICATIONS for:
  /**Returns a NodeIterator pointing to the first node of the graph
*/
  node_iterator node_begin() const { return node_iterator(this, 0); }
  /**Returns a NodeIterator pointing to one past the last node of the graph
*/
  node_iterator node_end() const { return node_iterator(this, num_nodes()); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    // Supply definitions AND SPECIFICATIONS for:
    /**Returns the edge the Incident Iterator is pointing to
*@return type Edge,the edge the Incident Iterator is pointing to
*@pre Edge is not equal to node.edge_end() and is valid
*/
    Edge operator*() const {
      // std::cout<<"node_id from incident iterator"<<node_id<<"\n";
      return Edge(g, g->i2u[node_id],
                  g->nodes[g->i2u[node_id]].edges[curr].node2_id);
    }
    /**For a valid iterator, Increments the interator to point to the next edge
    *incident on the node
       *@return the reference to an incremented incident iterator
    *@pre iterator should not be equal to the node.edge_end() and should be
    *valid
    */
    IncidentIterator &operator++() {
      ++curr;
      return *this;
    }
    /**Compares the given incident iterator to the this iterator. Returns true
    *if both are equal, else returns false
    *@param[in] @n , of IncidentIterator type.
    *@return a boolean, true if the iterators are equal, false otherwise
    */
    bool operator==(const IncidentIterator &i) const {
      if ((i.node_id == node_id) && (i.curr == curr))
        return true;
      else
        return false;
    }

  private:
    friend class Graph;
    graph_type *g;
    size_type node_id;
    size_type curr;
    IncidentIterator(const graph_type *g1, size_type nid, size_type curr1) {
      g = const_cast<graph_type *>(g1);
      node_id = nid;
      curr = curr1;
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
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /**Returns the edge the Edge Iterator is pointing to
*@return type Edge,the edge the Edge Iterator is pointing to
*@pre Edge is not equal to node.edge_end() and is valid
*/
    Edge operator*() const { return (*incident_iter); }
    /**For a valid iterator, Increments the interator to point to the next edge
    *of the graph
       *@return the reference to an incremented edge iterator
    *@pre iterator should not be equal to the node.edge_end() and should be
    *valid
    */
    EdgeIterator &operator++() {
      ++incident_iter;
      fix();
      return *this;
    }
    /**Compares the given edge iterator to the this iterator. Returns true if
    *both are equal, else returns false
    *@param[in] @n , of EdgeIterator type.
    *@return a boolean, true if the iterators are equal, false otherwise
    */
    bool operator==(const EdgeIterator &e) const {

      return ((e.g == g) && (e.node_iter == node_iter) &&
              (e.incident_iter == incident_iter));
    }

  private:
    friend class Graph;
    graph_type *g;
    NodeIterator node_iter;
    IncidentIterator incident_iter;
    EdgeIterator(const graph_type *g1, NodeIterator n1, IncidentIterator i1) {
      g = const_cast<graph_type *>(g1);
      node_iter = n1;
      incident_iter = i1;

      fix();
    }

    void fix() {
      int flag = 0;
      while (node_iter != g->node_end()) {
        while (incident_iter != (*node_iter).edge_end()) {
          if ((*node_iter) < (*incident_iter).node2()) {

            flag = 1;
            break;
          }
          ++incident_iter;
        }
        if (flag == 1) {
          break;
        }
        ++node_iter;

        if (node_iter == g->node_end()) {
          incident_iter = (*(g->node_begin())).edge_begin();
          break;
        }

        incident_iter = (*(node_iter)).edge_begin();
      }
    }
  };

  /**Returns a EdgeIterator pointing to the first edge of the graph
*/
  edge_iterator edge_begin() const {

    auto n = *node_begin();
    return edge_iterator(this, node_begin(), n.edge_begin());
  }
  /**Returns a EdgeIterator pointing to one past the last node of the graph
  */
  edge_iterator edge_end() const {
    auto n = *node_begin();

    return edge_iterator(this, node_end(), n.edge_begin());
  }
  /**Remove a single node and all its incident edges from the graph
  * @param[in,out] @a n, the node to be removed
  * @return 1(successful) if the node was deleted, 0(fail) otherwise
  * @post If remove was successful, Indices of all nodes with index < @a n are 
  *        not effected
  *       Indices of all nodes with index > @a n are moved to position index-1
  * @post new size() = old size() - 1, if remove was successful
  * @post neighbors of @a n have their degree reduced by 1, if remove was successful
  * @note invalidates @a n
  * At most O(num_nodes()) operations for sparse graphs
  */
  size_type remove_node(const Node &n) {
 
    if (!has_node(n)) {
      return 0;
    } 
    for (unsigned int i = 0; i < i2u.size(); i++) {
      if (nodes[i2u[i]].index > n.index()) {
        --nodes[i2u[i]].index;
      }
    }
    for (unsigned int i = 0; i < nodes[n.uid].edges.size(); ++i) {
      size_type neighbor = nodes[n.uid].edges[i].node2_id;
      for (unsigned int j = 0; j < nodes[neighbor].edges.size(); ++j) {
        if (nodes[neighbor].edges[j].node2_id == n.uid) {
          nodes[neighbor].edges[j] = nodes[neighbor].edges.back();
          nodes[neighbor].edges.pop_back();
          --number_of_edges;
          break;
        }
      }
    }
    
    i2u.erase(i2u.begin() + n.index());
    nodes[n.uid].index = -1;
    nodes[n.uid].edges.clear();
    return 1;
  }

  /**Remove a single node and all its incident edges from the graph
  * @param[in] @a n_it, an iterator node to be removed
  * @return iterator to the node following the removed node if the remove was   successful
  *         else return @a n_it
  * @pre old node_begin() <= @n_it < old node_end()
  * @post If remove is successful, Indices of all nodes with index < that of 
  *       node pointed by iterator are not effected
  *       Indices of all nodes with index > that of node pointed by iterator 
  *       are moved to position index-1
  * @post new size() = old size() - 1, if remove is successful
  * @post neighbors of @a *n_it have their degree reduced by 1
  * @note invalidates all iterators
  * At most O(num_nodes()) operations for sparse graphs
  */
  node_iterator remove_node(node_iterator n_it) {
    auto n = *n_it;
    remove_node(n);
    return n_it;
  }
  
  /**Remove a single edge 
  * @param[in,out] @a n1 @ n2, the nodes incident to the edge to be removed
  * @return 1(successful) if the edge was deleted, 0(fail) otherwise
  * @post new num_edges() = old num_edges() - 1, if remove was successful
  * At most O(num_nodes()) operations
  */
  size_type remove_edge(const Node &n1,
                        const Node &n2) { 

    if (!has_edge(n1, n2)) {
      return 0;
    }
    
    for (unsigned int i = 0; i < nodes[n1.uid].edges.size();
         ++i) { 
      if (nodes[n1.uid].edges[i].node2_id == n2.uid) {
        nodes[n1.uid].edges[i] =
            nodes[n1.uid].edges[nodes[n1.uid].edges.size() - 1];
        nodes[n1.uid].edges.pop_back();
        
        break;
      }
    }

    for (unsigned int i = 0; i < nodes[n2.uid].edges.size(); ++i) {
      if (nodes[n2.uid].edges[i].node2_id == n1.uid) {
        nodes[n2.uid].edges[i] =
            nodes[n2.uid].edges[nodes[n2.uid].edges.size() - 1];
        nodes[n2.uid].edges.pop_back();
        
        break;
      }
    }

    
      --number_of_edges;
      return 1;
   
  }
  /**Remove a single edge 
  * @param[in,out] @a e, the edge to be removed
  * @return 1(successful) if the edge was deleted, 0(fail) otherwise
  * @post new num_edges() = old num_edges() - 1, if remove was successful
  * At most O(num_nodes()) operations
  */
  size_type remove_edge(const Edge &e) {
    return remove_edge(e.node1(), e.node2());
  }
  
  /**Remove a single edge 
  * @param[in] @a e_it, iterator to the edge to be removed
  * @return iterator to the node following the removed node if the remove was successful
  *         else return @a e_it
  * @pre old edge_begin() <=@a e_it < edge_end()
  * @post new num_edges() = old num_edges() - 1, if remove was successful
  * @note invalidates all edge and incident iterators
  * At most O(num_nodes()) operations
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    auto e = *e_it;
    remove_edge(e);
    return e_it;
  }

private:
  //struct used for the proxy for edges
  struct internal_edge {
    size_type node2_id;
    edge_value_type val;
  };
  // struct used for the proxy for nodes
  struct internal_node {
    Point point;
    node_value_type val;
    long index;
    std::vector<internal_edge> edges;
  };

  // vector of nodes keeps track of the nodes belonging to a Graph
  std::vector<internal_node> nodes;
  // for edge sizes
  size_type number_of_edges = 0;
  std::vector<size_type> i2u;
};

#endif // CME212_GRAPH_HPP
