#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <tuple>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph : private totally_ordered<Graph<V,E>> {
 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

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

  /** Predeclaration of Node value type from template. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Predeclaration of Node value type from template. */
  using edge_value_type = E;

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

  struct edge_info {
    size_type neigh_;
    edge_value_type value_;
    edge_info(size_type neigh, edge_value_type value
      = edge_value_type()) : neigh_(neigh), value_(value) {}
  };

  struct node_info {
    Point p_;
    node_value_type value_;
    std::vector<edge_info> adj_;
    size_type idx_;
    node_info(Point p, size_type idx, node_value_type value = 
      node_value_type()) : p_(p), value_(value), idx_(idx) {}
  };

  /** Construct an empty graph. */
  Graph() {
    next_id_ = 0;
    size_edges_ = 0;
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[uid_].p_;
    }

    Point& position() {
      return graph_->nodes_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].idx_;
    }

    /** return a reference to the value of this (non-constant) node.  
        The value can be reassigned through this reference.
    */
    node_value_type& value() {
      return graph_->nodes_[uid_].value_;
    }

    /** return a reference to the value of this (constant) node.
        The value cannot be reassigned through this reference. 
    */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].value_;
    }

    /** return the number of edges adjacent to this node */
    size_type degree() const {
      return graph_->nodes_[uid_].adj_.size();
    }

    /** return incident iterator iin such that if this node has adjacent edges,
     * *iin is one of these edges; if not, iin is the null pointer. Also, any
     * other edge adjacent to this node can be given by incrementing iin
     * some finite number of times and dereferencing the result */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    /** return incident iterator iin_end such that iin_end is the null pointer,
     * and if we have an incident iterator iin of this node such that *iin = e
     * then incrementing iin some number of times will yield iin_end. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.graph_ == graph_) {
        if (n.uid_ == uid_) {
          return true;
        }
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
      return (std::tie(graph_, uid_) < std::tie(n.graph_, n.uid_));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type uid_;
    /** Construct a valid node.
     */
    Node(const graph_type* graph, size_type uid) {
      graph_ = const_cast<graph_type*>(graph);
      uid_ = uid; 
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = 
   node_value_type()) {
    Node n(this, next_id_);
    nodes_.push_back(node_info(position, i2u_.size(), value));
    i2u_.push_back(next_id_);
    ++next_id_;
    return n;     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this) {
      if (n.graph_->nodes_[n.uid_].idx_ < i2u_.size()) {
        return true;
      }
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
    return Node(this, i2u_[i]);
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, end_points_.first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, end_points_.second);
    }

    /** Return the value of this Edge. */
    edge_value_type& value() {
      unsigned index = 0;
      for (unsigned i = 0; i < graph_->nodes_[end_points_.first].adj_.size(); 
        ++i) {
        if (graph_->nodes_[end_points_.first].adj_[i].neigh_ == 
          end_points_.second) {
          index = i;
          break;
        }
      }
      return graph_->nodes_[end_points_.first].adj_[index].value_;
    }

    /** Return the value of this const Edge. */
    const edge_value_type& value() const {
      unsigned index = 0;
      for (unsigned i = 0; i < graph_->nodes_[end_points_.first].adj_.size(); 
        ++i) {
        if (graph_->nodes_[end_points_.first].adj_[i].neigh_ == 
          end_points_.second) {
          index = i;
          break;
        }
      }
      return graph_->nodes_[end_points_.first].adj_[index].value_;
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ != e.graph_) {
        return false;
      }
      if ((end_points_.first == e.end_points_.first) and 
        (end_points_.second == e.end_points_.second)) {
        return true;
      }
      if ((end_points_.first == e.end_points_.second) and 
        (end_points_.second == e.end_points_.first)) {
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
      return (std::tie(graph_, std::min(end_points_.first, end_points_.second),        std::max(end_points_.first, end_points_.second)) < 
        std::tie(e.graph_, std::min(e.end_points_.first, e.end_points_.second),        std::max(e.end_points_.first, e.end_points_.second)));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    std::pair<size_type, size_type> end_points_;
    /** Create valid edge using nodes. */
    Edge(const graph_type* graph, Node a, Node b) {
      graph_ = const_cast<graph_type*>(graph);
      end_points_.first = a.uid_;
      end_points_.second = b.uid_;
    }
    /** Create valid edge using uids. */
    Edge(const graph_type* graph, size_type aid, size_type bid) {
      graph_ = const_cast<graph_type*>(graph);
      end_points_.first = aid;
      end_points_.second = bid;
    }
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return size_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned i = 0; i < nodes_[a.uid_].adj_.size(); ++i) {
      if (nodes_[a.uid_].adj_[i].neigh_ == b.uid_) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = 
    edge_value_type()) {
    if (has_edge(a, b)) {
      return Edge(this, a, b);
    }
    else {
      nodes_[a.uid_].adj_.push_back(edge_info(b.uid_, value));
      nodes_[b.uid_].adj_.push_back(edge_info(a.uid_, value));
      ++size_edges_; 
      return Edge(this, a, b);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    next_id_ = 0;
    size_edges_ = 0;
    nodes_.clear();
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

    /** Dereference NodeIterator. 
     * @pre id_ < graph_.size()
     * @return node
     * @post node is a node of graph_. */
    Node operator*() const {
      return graph_->node(idx_);
    }

    /** Increment node Iterator. */
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /** if the private elements id_ and graph_ are
     the same, return true. Otherwise return false*/
    bool operator==(const NodeIterator& ni) const {
      return ((ni.graph_ == graph_) && (ni.idx_ == idx_));
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type idx_;
    /** Construct a valid Node Iterator. */
    NodeIterator(const graph_type* graph, size_type idx) {
      idx_ = idx;
      graph_ = graph;
    }
  };

    /** returns node iterator ni such that if graph_.size() > 0 then *ni is a
     * node; otherwise, ni is the nullpointer. Also, for every node in the graph     * incrementing ni some finite number of times and dereferencing yields 
     * this node. */
    node_iterator node_begin() const {
      return NodeIterator(this, 0);
    }

    /** returns node iterator ni_end such that ni_end is the null pointer, and 
      * if ni is a node iterator such that *ni = n, for some node n, then 
      * incrementing ni some finite number of times yields ni_end. */
    node_iterator node_end() const {
      return NodeIterator(this, size());
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

    /** @pre Node n1(graph_, node_id_).degree() < neigh_id_. 
     * @return Edge e. 
     * @post e is an edge incident on n1.  */
    Edge operator*() const {
      Node n1(graph_, node_id_);
      Node n2(graph_, graph_->nodes_[node_id_].adj_[neigh_id_].neigh_);
      return Edge(graph_, n1, n2);
    }

    /** Increment incident iterator */
    IncidentIterator& operator++() {
      ++neigh_id_;
      return *this;
    }
 
    /** return true if graph and node id are the same, false otherwise*/
    bool operator==(const IncidentIterator& iit) const {
      return ((iit.graph_ == graph_) && (iit.node_id_ == node_id_) && 
        (iit.neigh_id_ == neigh_id_)); 
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type node_id_;
    size_type neigh_id_;
    /** Construct valid Incident Iterator */
    IncidentIterator(const graph_type* graph, size_type node_id,
      size_type neigh_id) {
      graph_ = graph;
      node_id_ = node_id;
      neigh_id_ = neigh_id;
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

    /** @pre @a e_id_ < @a graph_->edges_.size(). 
     * @return Edge. 
     * @post return edge in the graph*/
    Edge operator*() const {
      return Edge(graph_, graph_->i2u_[node_idx_], 
        graph_->nodes_[graph_->i2u_[node_idx_]].adj_[neigh_idx_].neigh_);
    }

    /** Increment edge iterator */
    EdgeIterator& operator++() {
      ++neigh_idx_;
      while (node_idx_ < graph_->size()) {
        while (neigh_idx_ < 
          graph_->nodes_[graph_->i2u_[node_idx_]].adj_.size()) {
          if (graph_->nodes_[graph_->i2u_[node_idx_]].adj_[neigh_idx_].neigh_ 
            > graph_->i2u_[node_idx_]) {
            return *this;
          }
          ++neigh_idx_;
        }
        neigh_idx_ = 0;
        ++node_idx_;
      }
      return *this;
    }

    /** return true if graph and id are equal, false otherwise*/
    bool operator==(const EdgeIterator& ei) const {
      return ((ei.graph_ == graph_) && (ei.node_idx_ == node_idx_) && 
        (ei.neigh_idx_ == neigh_idx_));
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type node_idx_;
    size_type neigh_idx_;
    /** Construct valid Edge iterator*/
    EdgeIterator(const graph_type* graph, size_type node_idx, size_type 
      neigh_idx) : graph_(graph), node_idx_(node_idx), neigh_idx_(neigh_idx) {}
    
  };

  /** @return edge_iterator ei.
   * @post *ei == e for some edge e such that for each edge e' in the edges of
   * the graph there exists a nonnegative integer n such that incrementing ei 
   * n times will yield *ei == e'
   * */
  edge_iterator edge_begin() const {
    auto ei = EdgeIterator(this, 0, 0);
    while (ei.node_idx_ < ei.graph_->size()) {
      while (ei.neigh_idx_ < 
        ei.graph_->nodes_[ei.graph_->i2u_[ei.node_idx_]].adj_.size()) {
        if (ei.graph_->nodes_[ei.graph_->i2u_[ei.node_idx_]].adj_[ei.neigh_idx_].neigh_ > ei.graph_->i2u_[ei.node_idx_]) {
          return ei;
        }
        ++ei.neigh_idx_;
      }
      ei.neigh_idx_ = 0;
      ++ei.node_idx_;
    }
    return ei;
  }

  /** @return edge_iterator ei.
   * @post *ei is not a reference to an edge, but if ei' is the edge iterator
   * that has been incremented the most while still remaining an edge when 
   * dereferenced, then ++ei' == ei*/
  edge_iterator edge_end() const {
    return EdgeIterator(this, size(), 0);
  }

  /** Function to remove a given node and all of its incident edges
   *  from the graph.
   * @param n    Node to be removed.
   * @return     The uid that is now invalidated.
   * @pre @a n is a valid node of the graph.
   * @post all edges for which @a n was one of the end points are invalidated.
   * @post @a n is not a valid node.
   * @post some node indices may be changed.  We will still have the property
   *       that n1.index() < num_nodes() if n1 is a valid node.
   * @post any node iterators that were pointing to the last element in the 
   *       graph will be invalidated, as will the end iterator.
   * @post All edge iterators are invalidated.
   * Assuming the maximum degree k of the graph is small (meaning we can 
   * take k^2 as a constant), this function has a runtime of O(n).
   */
  size_type remove_node(const Node& n) {
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      Node n2 = e.node2();
      for (unsigned i = 0; i < nodes_[n2.uid_].adj_.size(); ++i) {
        if (nodes_[n2.uid_].adj_[i].neigh_ == n.uid_) {
          --size_edges_;
          nodes_[n2.uid_].adj_.erase(nodes_[n2.uid_].adj_.begin() + i);
        }
      }
    }
    i2u_.erase(i2u_.begin() + nodes_[n.uid_].idx_);
    for (auto ni = node_begin(); ni != node_end(); ++ni) {
      Node n1 = *ni;
      if (n1.uid_ > n.uid_) {
        --nodes_[n1.uid_].idx_;
      }
    }
    return n.uid_;
  }

  /** Function to remove a given node and all of its incident edges
   *  from the graph.
   * @param n_it Iterator such that *n_it is the node to be removed.
   * @return     An iterator pointing to the node that was after *@a n_it in
   *             the global ordering.  If *@a n_it was the last node, this is
   *             node_end().
   * @pre @a *n_it is a valid node of the graph.
   * @post all edges for which @a *n_it was an end point are invalidated.
   * @post *@a n_it is not a valid node.
   * @post some node indices may be changed.  We will still have the property
   *       that n1.index() < num_nodes() if n1 is a valid node.
   * @post any node iterators that were pointing to the last element in the 
   *       graph will be invalidated, as will the end iterator.
   * @post All edge iterators are invalidated.
   * Assuming the maximum degree k of the graph is small (meaning we can 
   * take k^2 as a constant), this function has a runtime of O(n).
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  /** Function to remove an edge with the given endpoints from the graph. 
   * @param n1   Endpoint of edge to be removed.
   * @param n2   Endpoint of edge to be removed.
   * @return     The number of edges that are now in the graph.
   * @pre @a n1 and @a n2 are valid nodes of the graph.
   * @post If there is an edge between @a n1 and @a n2, it will be invalidated.
   * @post All edge iterators are invalidated. 
   * This function takes O(k) time, where k is the maximum degree of the graph.
   */ 
  size_type remove_edge(const Node& n1, const Node& n2) {
    bool test = false;
    for (unsigned i = 0; i < nodes_[n1.uid_].adj_.size(); ++i) {
      if (n2.uid_ == nodes_[n1.uid_].adj_[i].neigh_) {
        test = true;
        nodes_[n1.uid_].adj_.erase(nodes_[n1.uid_].adj_.begin() + i);
        break;
      }
    } 
    for (unsigned i = 0; i < nodes_[n2.uid_].adj_.size(); ++i) {
      if (n1.uid_ == nodes_[n2.uid_].adj_[i].neigh_) {
        nodes_[n2.uid_].adj_.erase(nodes_[n2.uid_].adj_.begin() + i);
        break;
      }
    } 
    if (test) {
      --size_edges_;
    }
    return size_edges_;
  }

  /** Function to remove a given edge from the graph. 
   * @param e    Edge to be removed.
   * @return     The number of edges that are now in the graph.
   * @pre @a e is a valid edge of the graph. 
   * @post The edge @a e is invalidated.
   * @post All edge iterators are invalidated. 
   * This function takes O(k) time, where k is the maximum degree of the graph.
   */ 
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }
  /** Function to remove an edge with the given iterator from the graph. 
   * @param e_it Iterator such that *@a e_it is the edge to be removed. 
   * @return     An iterator pointing to the edge that was following *@a e_it
   *             In the graph.  If *@a e_it was the last edge, return
   *             edge_end().
   * @pre *@a e_it is a valid edge of the graph.
   * @post All edges equal to *@a e_it are invalidated.
   * @post All edge iterators are invalidated. 
   * This function takes O(k) time, where k is the maximum degree of the graph.
   */ 

  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

 private:
  size_type next_id_;
  size_type next_e_id_;
  size_type size_edges_ = 0;
  std::vector<node_info> nodes_;
  std::vector<size_type> i2u_;
  /** Find the index of an edge in a graph.
  * @pre @a a and @a b are distinct valid nodes of this graph for which 
  * has_edge(@a a, @a b) == true 
  * @return @a i, where edges_[i] is a pair containing @a a.id_ and @a b.id_ 
  */
   size_type find_edge_index(const Node& a, const Node& b) {
     int index = 0;
     Edge e1(this, a, b);
     for (auto ei = edge_begin(); ei != edge_end(); ++ei) {
       Edge e2 = *ei;
       if (e1 == e2) {
         return index;
       }
       ++index;
     }
     return size_edges_ + 1;
   }
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
