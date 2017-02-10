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
template <typename V>
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

  /** Type of value associated with node */
  using node_value_type = V;

 private:
  // Vector containing the points associated with the nodes
  std::vector<Point> points;

  // Vector containing the values associated with the nodes
  std::vector<node_value_type> values;

  // Vector containing the node pairs for each edge
  std::vector<std::pair<size_type,size_type>> edge_nodes;

  // Vector containing vectors of nodes adjacent to each node
  std::vector<std::vector<size_type>> adjacent_nodes;

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() :
    points(),
    values(),
    edge_nodes(),
    adjacent_nodes() {
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
      return graph_->points[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return a reference to the value associated with the node */
    node_value_type& value() {
      return (graph_->values[uid_]);
    }

    /** Return the value associated with the node */
    const node_value_type& value() const {
      return (graph_->values[uid_]);
    }

    /** Return the degree of a node, i.e. the number of edges
     *  incident to it */
    size_type degree() const {
      return (graph_->adjacent_nodes[uid_].size());
    }

    /** Iterator over edges incident to a given node, returns
     *  iterator associated with the first edge incident to the node */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,uid_,graph_->adjacent_nodes[uid_].begin());
    }

    /** Iterator over edges incident to a given node, returns
     *  iterator associated with the last edge incident to the node */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,uid_,graph_->adjacent_nodes[uid_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_==n.graph_ && uid_==n.uid_);
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
      if (graph_!=n.graph_) {
        return (graph_ < n.graph_);
      }
      else {
        return (uid_ < n.uid_);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions
    friend class Graph;
    // Pointer back to the graph container
    graph_type* graph_;
    // The node's index in the graph
    size_type uid_;
    // Private constructor of Node object
    Node(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return points.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    points.push_back(position);
    values.push_back(val);
    adjacent_nodes.push_back(std::vector<size_type>());
    return Node(this, points.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (points[n.index()] == n.position());
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
      return Node(graph_,node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,node2_);
    }

    // Return the index of an Edge
    size_type index() const {
      std::pair<size_type,size_type> pair = std::make_pair(this->node1().index(),
                                            this->node2().index());
      size_type edge_index = find(graph_->edge_nodes.begin(),graph_->
                             edge_nodes.end(),pair)-graph_->edge_nodes.begin();
      return edge_index;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_!=e.graph_) {
        return false;
      }
      else {
        return ((e.node1().index()==node1_ && e.node2().index()==node2_)
               ||(e.node1().index()==node2_ && e.node2().index()==node1_));
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_!=e.graph_) {
        return (graph_ < e.graph_);
      }
      else {
        return (this->index() < e.index());
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer to the graph that the Edge belongs to
    graph_type* graph_;
    // ID of the first node of the Edge
    size_type node1_;
    // ID of the second node of the Edge
    size_type node2_;
    // Private constructor for Edge
    Edge(const graph_type* graph,size_type node1,size_type node2)
      : graph_(const_cast<graph_type*>(graph)),node1_(node1),node2_(node2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_nodes.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this,edge_nodes[i].first,edge_nodes[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return adjacent_nodes[a.index()].find(b.index()) !=
           adjacent_nodes[a.index()].end();
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
    size_type node1 = a.index();
    size_type node2 = b.index();
    std::vector<size_type> search_array = adjacent_nodes[node1];
    size_type index = find(search_array.begin(),search_array.end(),
                      node2) - search_array.begin();
    if (index >= search_array.size()) {
      edge_nodes.push_back(std::make_pair(a.index(),b.index()));
      adjacent_nodes[node1].push_back(node2);
      adjacent_nodes[node2].push_back(node1);
    }
    return Edge(this,node1,node2);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points.clear();
    values.clear();
    edge_nodes.clear();
    adjacent_nodes.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Dereference operation, returns the node the iterator
     *  is pointing to */
    Node operator*() const {
      return Node(graph_,uid_);
    }

    /** Increment operation, increments uid by 1, returns
     *  a reference to the pointer to the new node */
    node_iterator& operator++() {
      uid_++;
      return *this;
    }

    /** Comparison operation, returns whether the iterator is
     *  currently pointing to the node we are interested in */
    bool operator==(const node_iterator& nit) const {
      return (graph_==nit.graph_ && uid_==nit.uid_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    // Pointer to the graph the node iterator belongs to
    graph_type* graph_;
    // Node's ID which the node iterator points to
    size_type uid_;
    // Private constructor called by the NodeIterator class
    NodeIterator(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Node iterator begin */
  node_iterator node_begin() const {
    return node_iterator(this,0);
  }

  /** Node iterator end */
  node_iterator node_end() const {
    return node_iterator(this,this->num_nodes());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Dereference operator, returns the edge the pointer is
     *  pointing to */
    Edge operator*() const {
      return Edge(graph_,uid_,*it_);
    }

    /** Increment operator, increments the position of the
     *  iterator of adjacent nodes by one */
    IncidentIterator& operator++() {
      it_++;
      return *this;
    }

    /** Comparison operator, compares whether two incident
     *  iterators are pointing to the same adjacent node */
    bool operator==(const IncidentIterator& iit) const {
      return (graph_==iit.graph_ && uid_==iit.uid_ &&
              it_==iit.it_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Pointer to the graph associated with the incident iterator
    graph_type* graph_;
    // ID of the node we are interested in iterating over its incident edges
    size_type uid_;
    // Iterator over the vector containing the adjacent nodes
    std::vector<size_type>::iterator it_;
    // Private constructor for IncidentIterator
    IncidentIterator(const graph_type* graph, size_type uid,
                     std::vector<size_type>::iterator it)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid), it_(it) {
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator>  {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Define edge dereference operator */
    Edge operator*() const {
      return Edge(graph_,graph_->edge_nodes[uid_].first,
                  graph_->edge_nodes[uid_].second);
    }

    /** Increment operator, returns a reference to the EdgeIterator
     *  and increments the uid_ the EdgeIterator refers to */
    EdgeIterator& operator++() {
      uid_++;
      return *this;
    }

    /** Comparison operator */
    bool operator==(const EdgeIterator& eit) const {
      return (graph_==eit.graph_ && uid_==eit.uid_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Pointer to the graph the edge iterator belongs to
    graph_type* graph_;
    // Edge's ID which the iterator points to
    size_type uid_;
    // Private constructor called by the EdgeIterator class
    EdgeIterator(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Define the start of the edge iterator */
  edge_iterator edge_begin() const {
    return edge_iterator(this,0);
  }

  /** Define the end of the edge iterator */
  edge_iterator edge_end() const {
    return edge_iterator(this,this->num_edges());
  }

};

#endif // CME212_GRAPH_HPP
