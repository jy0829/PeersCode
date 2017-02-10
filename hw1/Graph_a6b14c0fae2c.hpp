#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {}

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
      g = nullptr;
      i = -1;
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(g != nullptr && i < g->size());
      return g->point_list[i];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(g != nullptr && i < g->size());
      return i;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /* Return the value associated with this node. */
    node_value_type& value() {
      assert(g != nullptr && i < g->size());
      return g->value_list[i];
    }

    /* Return the value associated with this node by const reference. */
    const node_value_type&& value() const {
      assert(g != nullptr && i < g->size());
      return (const node_value_type&&)g->value_list[i];
    }

    /* Return the degree of this node. */
    size_type degree() const {
      assert(g != nullptr && i < g->size());
      return g->neighbor_list_list[i].size();
    }

    /* Return the begin iterator of the incident edges of this node. */
    incident_iterator edge_begin() const {
      return incident_iterator(g, i, g->neighbor_list_list[i].begin());
    }

    /* Return the end iterator of the incident edges of this node. */
    incident_iterator edge_end() const {
      return incident_iterator(g, i, g->neighbor_list_list[i].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return g == n.g && i == n.i;
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
      return g == n.g && i < n.i;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(const Graph* g, size_type i) : g(const_cast<Graph*>(g)), i(i) {}
    
    Graph* g;
    size_type i;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return point_list.size();
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
    point_list.push_back(position);
    value_list.push_back(value);
    neighbor_list_list.emplace_back();
    return Node(this, point_list.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.g;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < point_list.size());
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(g, n1_index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g, n2_index);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return g == e.g && ((n1_index == e.n1_index && n2_index == e.n2_index) ||
        (n1_index == e.n2_index && n2_index == e.n1_index));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type min_index = std::min(n1_index, n2_index);
      size_type max_index = std::max(n1_index, n2_index);
      size_type e_min_index = std::min(e.n1_index, e.n2_index);
      size_type e_max_index = std::max(e.n1_index, e.n2_index);
      
      return g == e.g && ((min_index < e_min_index) ||
        (min_index == e_min_index && max_index < e_max_index));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* g, size_type n1_index, size_type n2_index) :
      g(const_cast<Graph*>(g)), n1_index(n1_index), n2_index(n2_index) {}

    Graph* g;
    size_type n1_index;
    size_type n2_index;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < edge_list.size());
    size_type n1_index = edge_list[i].first;
    size_type n2_index = edge_list[i].second;
    return Edge(this, n1_index, n2_index);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (a.g != this || b.g != this) {
      return false;
    }
    int n1_index = std::min(a.i, b.i);
    int n2_index = std::max(a.i, b.i);
    return edge_set.find(std::make_pair(n1_index, n2_index)) != edge_set.end();
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
    assert(a.g == this && a.i < size() && b.g == this && b.i < size());
    Edge newEdge(this, a.i, b.i);
    if (!has_edge(a, b)) {
      int n1_index = std::min(a.i, b.i);
      int n2_index = std::max(a.i, b.i);
      auto index_pair = std::make_pair(n1_index, n2_index);
      edge_set.insert(index_pair);
      edge_list.push_back(index_pair);
      neighbor_list_list[a.i].push_back(b.i);
      if (a.i != b.i) {
        neighbor_list_list[b.i].push_back(a.i);
      }
    }
    return newEdge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    point_list.clear();
    value_list.clear();
    edge_list.clear();
    edge_set.clear();
    neighbor_list_list.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    // Return the dereferenced node.
    Node operator*() const {
      return g->node(i);
    }

    // Forward one step.
    NodeIterator& operator++() {
      ++i;
      return *this;
    }

    bool operator==(const NodeIterator& iter) const {
      return g == iter.g && i == iter.i;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    NodeIterator(const Graph* g, size_type i) : g(g), i(i) {}

    const Graph* g;
    size_type i;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  node_iterator node_end() const {
    return node_iterator(this, size());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    Edge operator*() const {
      return Edge(g, i, *iter);
    }

    IncidentIterator& operator++() {
      ++iter;
      return *this;
    }

    bool operator==(const IncidentIterator& a) const {
      return g == a.g && i == a.i && iter == a.iter;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    IncidentIterator(const Graph* g, const size_type i,
      typename std::vector<size_type>::const_iterator iter) :
      g(const_cast<Graph*>(g)), i(i), iter(iter) {}

    Graph* g;
    const size_type i;
    typename std::vector<size_type>::const_iterator iter;
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const {
      return Edge(g, (*iter).first, (*iter).second);
    }

    EdgeIterator& operator++() {
      ++iter;
      return *this;
    }

    bool operator==(const EdgeIterator& a) const {
      return g == a.g && iter == a.iter;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    EdgeIterator(const Graph* g,
      typename std::vector<std::pair<size_type, size_type>>::const_iterator iter) :
      g(const_cast<Graph*>(g)), iter(iter) {}

    Graph* g;
    typename std::vector<std::pair<size_type, size_type>>::const_iterator iter;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const {
    return edge_iterator(this, edge_list.begin());
  }

  edge_iterator edge_end() const {
    return edge_iterator(this, edge_list.end());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> point_list;
  std::vector<V> value_list;
  std::vector<std::pair<size_type, size_type>> edge_list;
  std::set<std::pair<size_type, size_type>> edge_set;
  std::vector<std::vector<size_type>> neighbor_list_list;
};

#endif // CME212_GRAPH_HPP
