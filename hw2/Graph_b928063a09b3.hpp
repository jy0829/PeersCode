#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <functional>
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
class Graph {
  
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
  /** Synonym for Node value. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Synonym for Node value. */
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

 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  
  /** @nodeinfo encapsulate the paramaters included in @nodes_. */
  struct nodeinfo {
    size_type node_idx_;
    Point position_;
    node_value_type value_;
    nodeinfo(size_type node_idx, Point position, node_value_type value)
        : node_idx_(node_idx), position_(position), value_(value) {}
  };
 
  /** @edgeinfo encapsulate the paramaters included in @edges_. */
  struct edgeinfo {
    size_type node1_uid_;
    size_type node2_uid_;
    edge_value_type value_;
    edgeinfo(size_type node1_uid, size_type node2_uid, edge_value_type value)
        : node1_uid_(node1_uid), node2_uid_(node2_uid), value_(value) {}
  };
  
 public:
 
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
      : nodes_(), nodes_i2u_(), adj_() {
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[node_uid_].position_;
    }
    
    /** Return this node's modifiable position. */
    Point& position() {
      return graph_->nodes_[node_uid_].position_;
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[node_uid_].node_idx_;
    }

    /** Return this node's modifiable value. */
    node_value_type& value() {
      return graph_->nodes_[node_uid_].value_;
    }
    
    /** Return this node's value. */
    const node_value_type& value() const {
      return graph_->nodes_[node_uid_].value_;
    }
    
    /** Return the number of edges incident to this node. */
    size_type degree() const {
      return graph_->adj_[node_uid_].size();
    }
    
    /** Return incident iterator to the beginning of edges incident to this node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, node_uid_, 0);
    }
    
    /** Return incident iterator to the end of edges incident to this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, node_uid_, degree());
    }
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ and node_uid_ == n.node_uid_;
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
      std::less<Graph*> compare;
      return (graph_ == n.graph_ and node_uid_ < n.node_uid_)
             or compare(graph_, n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    /** Pointer back to the Graph container. */
    Graph* graph_;
    /** This node's unique identification number. */
    size_type node_uid_;
    /** Private Constructor. */
    Node(const Graph* graph, size_type node_uid)
        : graph_(const_cast<Graph*>(graph)), 
          node_uid_(node_uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_i2u_.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    nodes_.push_back(nodeinfo(num_nodes(), position, value));
    nodes_i2u_.push_back(nodes_.size() - 1);
    adj_.push_back(std::vector<edgeinfo>());
    return Node(this, nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_ and n.index() < num_nodes();
  }

  /** Remove a node from the graph.
   * @return  num_nodes() if !has_node(@a n), 0 otherwise.
   *
   * @pre   @a n is a valid node of the graph, whether in it or not.
   * @post  All edges linked to @a n are removed.
   * @post  new @a n.node_idx_ is set to num_nodes() so that !has_node(n).
   * @post  All nodes with node_idx_ > old @a n.node_idx_ have their indices
   *        reduced by 1.
   * @post  The size of new @nodes_i2u_ is reduced by 1.
   */
  size_type remove_node(const Node& n) {
    size_type node_idx = n.index();
    // Stop function if node isn't in graph
    if (!has_node(n))
      return num_nodes();
    // Remove the edges attached to the node
    for (incident_iterator iit = n.edge_begin() ; iit != n.edge_end() ; )
      remove_edge(*iit);
    adj_[n.node_uid_].clear();
    // Remove the node from nodes in graph
    nodes_i2u_.erase(nodes_i2u_.begin() + node_idx);
    // Update indices for nodes in graph, invalidating removed node
    for (auto it = nodes_.begin() ; it != nodes_.end() ; ++it) {
      if (it->node_idx_ == node_idx)
        it->node_idx_ = num_nodes();
      else if (!(it->node_idx_ < node_idx))
        it->node_idx_ -= 1;
    }
    return 0;
  }
  
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, nodes_i2u_[i]);
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
    Node node1() const {
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_uid_);
    }

    /** Return this edge's modifiable value. */
    edge_value_type& value() {
      for (auto it = graph_->adj_[node1_uid_].begin() ; it != graph_->adj_[node1_uid_].end() ; ++it)
        if (it->node2_uid_ == node2_uid_)
          return it->value_;
    }
    
    /** Return this edge's value. */
    const edge_value_type& value() const {
      for (auto it = adj_[node1_uid_].begin() ; it != adj_[node1_uid_].end() ; ++it)
        if (it->node2_uid_ == node2_uid_)
          return it->value_;
    }
    
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(e.node1() < e.node2());
      return (graph_ == e.graph_
              and node1_uid_ == e.node1_uid_
              and node2_uid_ == e.node2_uid_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(e.node1() < e.node2());
      std::less<Graph*> compare;
      return (graph_ == e.graph_
              and (node1_uid_ == e.node1_uid_ ? node2_uid_ < e.node2_uid_ : node1_uid_ < e.node1_uid_))
          or compare(graph_, e.graph_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    /** Pointer back to the Graph container. */
    Graph* graph_;
    /** This edge's nodes' identification numbers. */
    size_type node1_uid_;
    size_type node2_uid_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
        : graph_(const_cast<Graph*>(graph)), 
          node1_uid_(node1_uid),
          node2_uid_(node2_uid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type n = 0;
    for (auto it = adj_.begin() ; it != adj_.end() ; ++it)
      n += it->size();
    return n/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    edge_iterator eit = edge_begin();
    for (size_type count = 0 ; count < i ; ++count)
      ++eit;
    return (*eit);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto it = adj_[a.node_uid_].begin() ; it != adj_[a.node_uid_].end() ; ++it)
      if (it->node2_uid_ == b.node_uid_)
        return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    if (!has_edge(a, b)) {
      adj_[a.node_uid_].push_back(edgeinfo(a.node_uid_, b.node_uid_, value));
      adj_[b.node_uid_].push_back(edgeinfo(b.node_uid_, a.node_uid_, value));
    }
    return Edge(this, a.node_uid_, b.node_uid_);
  }
  
  /** Remove an edge from the graph.
   * @return  num_edges() if !has_edge(@a a, @a b), 0 otherwise.
   *
   * @pre   The edge e between @a a and @a b is a valid edge of the graph,
   *        whether in it or not.
   * @post  The edges edge(@a a, @a b) and edge(@a b, @a a) are removed from adj_.
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // Stop function if edge isn't in graph
    if (!has_edge(a, b))
      return num_edges();
    // Remove edge from adj_
    for (auto it = adj_[a.node_uid_].begin() ; it != adj_[a.node_uid_].end() ; ++it) {
      if (it->node2_uid_ == b.node_uid_) {
        adj_[a.node_uid_].erase(it);
        break;
      }  
    }
    for (auto it = adj_[b.node_uid_].begin() ; it != adj_[b.node_uid_].end() ; ++it) {
      if (it->node2_uid_ == a.node_uid_) {
        adj_[b.node_uid_].erase(it);
        break;
      }  
    }
    return 0;
  }
  
  /** Remove an edge from the graph.
   * @return  num_edges() if !has_edge(@a a, @a b), 0 otherwise.
   *
   * @pre   The edge e between @a e.node1() and @a e.node2() is a valid
   *        edge of the graph, whether in it or not.
   * @post  The edges edge(@a e.node1(), @a e.node2()) and
   *        edge(@a e.node2(), @a e.node1*() are removed from adj_.
   */
  size_type remove_edge(const Edge& e) {
    // Stop function if edge isn't in graph
    if (!has_edge(e.node1(), e.node2()))
      return num_edges();
    // Remove edge from adj_
    for (auto it = adj_[e.node1_uid_].begin() ; it != adj_[e.node1_uid_].end() ; ++it) {
      if (it->node2_uid_ == e.node2_uid_) {
        adj_[e.node1_uid_].erase(it);
        break;
      }  
    }
    for (auto it = adj_[e.node2_uid_].begin() ; it != adj_[e.node2_uid_].end() ; ++it) {
      if (it->node2_uid_ == e.node1_uid_) {
        adj_[e.node2_uid_].erase(it);
        break;
      }  
    }
    return 0;
  }
  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    nodes_i2u_.clear();
    adj_.clear();
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
    NodeIterator() {}

    /** Dereferences node iterator.
     * @pre 0 <= @a node_it_ < num_nodes()
     * @post result_node.index() == @a node_it_
     */
    Node operator*() const {
      return graph_->node(node_it_);
    }
    
    /** Increment node iterator. */
    NodeIterator& operator++() {
      node_it_++;
      return *this;
    }
    
    /** Compare two node iterators. */
    bool operator==(const NodeIterator& nit) const {
      return graph_ == nit.graph_ and node_it_ == nit.node_it_;
    }

   private:
    friend class Graph;
    /** Pointer back to the Graph container. */
    Graph* graph_;
    /** ID of the node the iterator is looking at. */
    size_type node_it_;
    /** Private constructor. */
    NodeIterator(const Graph* graph, size_type node_it)
      : graph_(const_cast<Graph*>(graph)),
        node_it_(node_it) {}
  };

  /** Return node iterator to the beginning of nodes_. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  
  /** Return node iterator to the end of nodes_. */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  /** Remove the node pointed to from the graph.
   * @return  node_end() if !has_node(@a (*nit)), pointer to the position of
   *          old @a nit otherwise.
   *
   * @pre  @a (*nit) is a valid node of the graph, whether in it or not.
   */
  node_iterator remove_node(node_iterator nit) {
    size_type removed = remove_node((*nit));
    if (removed == num_nodes())
      return NodeIterator(this, num_nodes());
    return nit;
  }
  
  //
  // Incident Iterator
  //
  
  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** Dereferences incident iterator. */
    Edge operator*() const {
      return Edge(graph_, node_uid_, graph_->adj_[node_uid_][incident_it_].node2_uid_);
    }
    
    /** Increment incident iterator. */
    IncidentIterator& operator++() {
      ++incident_it_;
      return *this;
    }
    
    /** Compare two incident iterators. */
    bool operator==(const IncidentIterator& iit) const {
      return graph_ == iit.graph_ 
             and node_uid_ == iit.node_uid_ 
             and incident_it_ == iit.incident_it_;
    }

   private:
    friend class Graph;
    /** Pointer back to the Graph container. */
    Graph* graph_;
    /** ID of the node the iterator is looking at. */
    size_type node_uid_;
    /** Position of the edge the iterator is looking at amoung edges incident to the node. */
    size_type incident_it_;
    /** Private constructor. */
    IncidentIterator(const Graph* graph, size_type node_uid, size_type incident_it)
      : graph_(const_cast<Graph*>(graph)),
        node_uid_(node_uid),
        incident_it_(incident_it) {}
  };

  //
  // Edge Iterator
  //
  
  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Dereferences edge iterator.
     * @pre 0 <= @a edge_it_ < num_edges()
     * @post result_edge.index() == @a edge_it_
     */
    Edge operator*() const{
      return Edge(graph_, node_uid_, graph_->adj_[node_uid_][incident_it_].node2_uid_);
    }
    
    /** Increment edge iterator. */
    EdgeIterator& operator++() {
      do {
        if (graph_->adj_[node_uid_].size() != 0) {
          for (++incident_it_ ; incident_it_ < graph_->adj_[node_uid_].size() ; ++incident_it_)
            if (!(graph_->adj_[node_uid_][incident_it_].node2_uid_ < node_uid_))
              return *this;
        }
        incident_it_ = 0;
        node_uid_++;
      } while (node_uid_ != graph_->num_nodes());
      return *this;
      
    }
    
    /** Compare two edge iterators. */
    bool operator==(const EdgeIterator& eit) const {
      return graph_ == eit.graph_ 
             and node_uid_ == eit.node_uid_
             and incident_it_ == eit.incident_it_;
    }
    
   private:
    friend class Graph;
    /** Pointer back to the Graph container. */
    Graph* graph_;
    /** ID of the first node of the edge the iterator is looking at. */
    size_type node_uid_;
    /** Position of the edge the iterator is looking at amoung edges 
     *  incident to its first node. */
    size_type incident_it_;
    /** Private constructor. */
    EdgeIterator(const Graph* graph, size_type node_uid, size_type incident_it)
      : graph_(const_cast<Graph*>(graph)),
        node_uid_(node_uid),
        incident_it_(incident_it) {}
  };

  /** Return edge iterator to the beginning of edges_. */
  edge_iterator edge_begin() const {
    for (size_type i = 0 ; i < num_nodes() ; ++i)
      if (adj_[i].size() != 0)
        return EdgeIterator(this, i, 0);
    return EdgeIterator(this, num_nodes(), 0);
  }
  
  /** Return edge iterator to the end of edges_. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_nodes(), 0);
  }
  
  /** Remove the edge pointed to from the graph.
   * @return  edge_end() if !has_edge(@a (*nit).node1(), @a (*nit).node2()),
   *          pointer to the position of old @a eit otherwise.
   *
   * @pre  @a (*eit) is a valid edge of the graph, whether in it or not.
   */
  edge_iterator remove_edge(edge_iterator eit) {
    size_type removed = remove_edge((*eit));
    if (removed == num_edges())
      return EdgeIterator(this, num_nodes(), 0);
    return eit;
  }
  
 private:
 
  /** @nodes_ is a vector of pairs containing the Point refered by the node
   * and the value of typevalue V corresponding to the node, both called by proxy.
   */
  std::vector<nodeinfo> nodes_;
  
  /** @nodes_i2u_ is a vector mapping nodes indices to their uid. */
  std::vector<size_type> nodes_i2u_;
  
  /** @adj_ is a vector containing vector of pairs. It contains
   * at position i the pairs of the id of each node adjacent to node i and the
   * respective edge id that connects node i to that node.
   */
  std::vector<std::vector<edgeinfo>> adj_;
  
  
};

#endif // CME212_GRAPH_HPP
