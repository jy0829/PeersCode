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
template <typename V>
class Graph {
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
  /** Synonym for Node value. */
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
  Graph() 
      : nodes_(), edges_(), adj_() {
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
      return graph_->nodes_[node_id_].first;
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_id_;
    }

    /** Default node value */
    node_value_type& value() {
      return graph_->nodes_[node_id_].second;
    }
    
    /** Return this node's value. */
    const node_value_type& value() const {
      return graph_->nodes_[node_id_].second;
    }
    
    /** Return the number of edges incident to this node. */
    size_type degree() const {
      return graph_->adj_[node_id_].size();
    }
    
    /** Return incident iterator to the beginning of edges incident to this node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, node_id_, 0);
    }
    
    /** Return incident iterator to the end of edges incident to this node. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, node_id_, degree());
    }
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ and node_id_ == n.index();
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
      return (graph_ == n.graph_ and node_id_ < n.index())
             or compare(graph_, n.graph_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    /** Pointer back to the Graph container. */
    Graph* graph_;
    /** This node's identification number. */
    size_type node_id_;
    /** Private Constructor. */
    Node(const Graph* graph, size_type node_id)
        : graph_(const_cast<Graph*>(graph)), 
          node_id_(node_id) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
    nodes_.push_back(std::make_pair(position, value));
    adj_.push_back(std::vector<std::pair<size_type, size_type>>());
    return Node(this, nodes_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_);
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(e.node1() < e.node2());
      return graph_ == e.graph_
             and ((node1_id_ == e.node1() and node2_id_ == e.node2())
                 or (node2_id_ == e.node1() and node1_id_ == e.node2()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(e.node1() < e.node2());
      std::less<Graph*> compare;
      return (graph_ == e.graph_ and edge_id_ < e.index())
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
    /** This edge's identification number. */
    size_type edge_id_;
    /** This edge's nodes' identification numbers. */
    size_type node1_id_;
    size_type node2_id_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type edge_id, size_type node1_id, size_type node2_id)
        : graph_(const_cast<Graph*>(graph)), 
          edge_id_(edge_id),
          node1_id_(node1_id),
          node2_id_(node2_id) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i, edges_[i].first, edges_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned int i = 0; i < adj_[a.index()].size(); i++)
      if (b.index() == adj_[a.index()][i].first)
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
  Edge add_edge(const Node& a, const Node& b) {
    for (unsigned int i = 0; i < adj_[a.index()].size(); i++)
      if (b.index() == adj_[a.index()][i].first)
        return Edge(this, adj_[a.index()][i].second, a.index(), b.index());
    if (a < b)
      edges_.push_back(std::make_pair(a.index(), b.index()));
    else
      edges_.push_back(std::make_pair(b.index(), a.index()));
    adj_[a.index()].push_back(std::make_pair(b.index(), edges_.size()-1));
    adj_[b.index()].push_back(std::make_pair(a.index(), edges_.size()-1));
    return Edge(this, edges_.size()-1, a.index(), b.index());
  }
  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
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
    NodeIterator() {
    }

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
        node_it_(node_it) {
    }
  };

  /** Return node iterator to the beginning of nodes_. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  
  /** Return node iterator to the end of nodes_. */
  node_iterator node_end() const {
    return NodeIterator(this, size());
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
    IncidentIterator() {
    }

    /** Dereferences incident iterator. */
    Edge operator*() const {
      Edge unordered_edge = graph_->edge(graph_->adj_[node_id_][incident_it_].second);
      if (unordered_edge.node1() == graph_->node(node_id_))
        return graph_->add_edge(graph_->node(node_id_), unordered_edge.node2());
      else
        return graph_->add_edge(graph_->node(node_id_), unordered_edge.node1());
    }
    
    /** Increment incident iterator. */
    IncidentIterator& operator++() {
      incident_it_++;
      return *this;
    }
    
    /** Compare two incident iterators. */
    bool operator==(const IncidentIterator& iit) const {
      return graph_ == iit.graph_ 
             and node_id_ == iit.node_id_ 
             and incident_it_ == iit.incident_it_;
    }

   private:
    friend class Graph;
    /** Pointer back to the Graph container. */
    Graph* graph_;
    /** ID of the node the iterator is looking at. */
    size_type node_id_;
    /** Position of the edge the iterator is looking at amoung edges incident to the node. */
    size_type incident_it_;
    /** Private constructor. */
    IncidentIterator(const Graph* graph, size_type node_id, size_type incident_it)
      : graph_(const_cast<Graph*>(graph)),
        node_id_(node_id),
        incident_it_(incident_it) {
    }
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
    EdgeIterator() {
    }

    /** Dereferences edge iterator.
     * @pre 0 <= @a edge_it_ < num_edges()
     * @post result_edge.index() == @a edge_it_
     */
    Edge operator*() const{
      return graph_->edge(edge_it_);
    }
    
    /** Increment edge iterator. */
    EdgeIterator& operator++() {
      edge_it_++;
      return *this;
    }
    
    /** Compare two edge iterators. */
    bool operator==(const EdgeIterator& eit) const {
      return graph_ == eit.graph_ and edge_it_ == eit.edge_it_;
    }
    
   private:
    friend class Graph;
    /** Pointer back to the Graph container. */
    Graph* graph_;
    /** ID of the edge the iterator is looking at. */
    size_type edge_it_;
    /** Private constructor. */
    EdgeIterator(const Graph* graph, size_type edge_it)
      : graph_(const_cast<Graph*>(graph)),
        edge_it_(edge_it) {
    }
  };

  /** Return edge iterator to the beginning of edges_. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  
  /** Return edge iterator to the end of edges_. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }
  
 private:

  /** @nodes_ is a vector of pairs containing the Point refered by the node
   * and the value of typevalue V corresponding to the node, both called by proxy.
   */
  std::vector<std::pair<Point, V>> nodes_;
  /** @edges_ is a vector of pairs containing the 2 nodes id of an edge.
   * We choose to add them with Node1 < Node2 with the < operator defined for nodes.
   */
  std::vector<std::pair<size_type, size_type>> edges_;
  /** @adj_ is a vector containing vector of pairs. It contains
   * at position i the pairs of the id of each node adjacent to node i and the
   * respective edge id that connects node i to that node.
   */
  std::vector<std::vector<std::pair<size_type, size_type>>> adj_;
  
};

#endif // CME212_GRAPH_HPP
