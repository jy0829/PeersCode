#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <map>
#include <memory>

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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct NodeInfo;
  struct EdgeInfo;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
  using edge_value_type = E;

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

  using uid_type = unsigned;

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
      uid = -1;
    }

    /** Return this node's position. */
    const Point& position() const {
      return g->node_info_list[index()]->position;
    }

    Point& position() {
      return g->node_info_list[index()]->position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(g != nullptr && g->node_uid_to_node_idx_map.count(uid) == 1);
      return g->node_uid_to_node_idx_map[uid];
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
      return g->node_info_list[index()]->value;
    }

    /* Return the value associated with this node by const reference. */
    const node_value_type&& value() const {
      return (const node_value_type&&)g->node_info_list[index()]->value;
    }

    /* Return the degree of this node. */
    size_type degree() const {
      return g->node_info_list[index()]->neighbor_uid_set.size();
    }

    /* Return the begin iterator of the incident edges of this node. */
    incident_iterator edge_begin() const {
      return incident_iterator(g, uid,
          g->node_info_list[index()]->neighbor_uid_set.begin());
    }

    /* Return the end iterator of the incident edges of this node. */
    incident_iterator edge_end() const {
      return incident_iterator(g, uid,
          g->node_info_list[index()]->neighbor_uid_set.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return g == n.g && uid == n.uid;
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
      return g < n.g || (g == n.g && uid < n.uid);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(const Graph* g, uid_type uid) : g(const_cast<Graph*>(g)), uid(uid) {}
    
    Graph* g;
    uid_type uid;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_info_list.size();
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
  Node add_node(const Point& position,
        const node_value_type& value = node_value_type()) {
    uid_type node_uid = next_node_uid_to_use;
    ++next_node_uid_to_use;
    node_uid_to_node_idx_map[node_uid] = num_nodes();
    node_info_list.emplace_back(new NodeInfo(node_uid, position, value));
    return Node(this, node_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.g && node_uid_to_node_idx_map.count(n.uid) == 1;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, node_info_list[i]->uid);
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
      return Node(g, uid1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g, uid2);
    }

    /* Return the length of the edge. */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /* 
     * Return the index of this edge.
     */
    size_type index() const {
      auto edge_uid = std::make_pair(
          std::min(uid1, uid2), std::max(uid1, uid2));
      assert(g != nullptr && g->edge_uid_to_edge_idx_map.count(edge_uid) == 1);
      return g->edge_uid_to_edge_idx_map[edge_uid];
    }

    /* 
     * Return the edge value.
     */
    edge_value_type& value() {
      return g->edge_info_list[index()]->value;
    }

    /* 
     * Return the edge value by const reference.
     */
    const edge_value_type& value() const {
      return g->edge_info_list[index()]->value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return g == e.g && ((uid1 == e.uid1 && uid2 == e.uid2) ||
          (uid1 == e.uid2 && uid2 == e.uid1));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      uid_type min_uid = std::min(uid1, uid2);
      uid_type max_uid = std::max(uid1, uid2);
      uid_type e_min_uid = std::min(e.uid1, e.uid2);
      uid_type e_max_uid = std::max(e.uid1, e.uid2);
      
      return g < e.g || (g == e.g && ((min_uid < e_min_uid) ||
          (min_uid == e_min_uid && max_uid < e_max_uid)));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* g, uid_type uid1, uid_type uid2) :
        g(const_cast<Graph*>(g)), uid1(uid1), uid2(uid2) {}

    Graph* g;
    uid_type uid1;
    uid_type uid2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_info_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    size_type uid1 = edge_info_list[i]->uid1;
    size_type uid2 = edge_info_list[i]->uid2;
    return Edge(this, uid1, uid2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    uid_type uid1 = std::min(a.uid, b.uid);
    uid_type uid2 = std::max(a.uid, b.uid);
    return edge_uid_to_edge_idx_map.count(std::make_pair(uid1, uid2)) == 1;
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
    if (!has_edge(a, b)) {
      uid_type uid1 = std::min(a.uid, b.uid);
      uid_type uid2 = std::max(a.uid, b.uid);
      edge_uid_to_edge_idx_map[std::make_pair(uid1, uid2)] = num_edges();
      edge_info_list.emplace_back(new EdgeInfo(uid1, uid2));
      size_type index1 = node_uid_to_node_idx_map[uid1];
      size_type index2 = node_uid_to_node_idx_map[uid2];
      node_info_list[index1]->neighbor_uid_set.insert(uid2);
      node_info_list[index2]->neighbor_uid_set.insert(uid1);
    }
    return Edge(this, a.uid, b.uid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    next_node_uid_to_use = 0;
    node_info_list.clear();
    node_uid_to_node_idx_map.clear();
    edge_info_list.clear();
    edge_uid_to_edge_idx_map.clear();
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
      return Edge(g, uid, *iter);
    }

    IncidentIterator& operator++() {
      ++iter;
      return *this;
    }

    bool operator==(const IncidentIterator& a) const {
      return g == a.g && uid == a.uid && iter == a.iter;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    IncidentIterator(const Graph* g, const uid_type uid,
      typename std::set<uid_type>::const_iterator iter) :
          g(const_cast<Graph*>(g)), uid(uid), iter(iter) {}

    Graph* g;
    const uid_type uid;
    typename std::set<uid_type>::const_iterator iter;
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
      return g->edge(i);
    }

    EdgeIterator& operator++() {
      ++i;
      return *this;
    }

    bool operator==(const EdgeIterator& iter) const {
      return g == iter.g && i == iter.i;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    EdgeIterator(const Graph* g, size_type i) : g(g), i(i) {}

    const Graph* g;
    size_type i;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  edge_iterator edge_end() const {
    return edge_iterator(this, num_edges());
  }

  /* Remove the node @a n.
   * @return true if the @a n is in the graph, false otherwise.
   * @post new has_node(@a n) == false
   * @post new num_edges() == old num_edges() + old n.degree();
   * Can invalidate edge indices and iterators of edges coming after
   *     any of the removed edges associated with the node in the old
   *     order.
   * Can invalidate node indices and iterators coming after the
   *     deleted node in the old order.
   * Complexity: O(max_degree * log(old size())).
   */
  bool remove_node(const Node& n) {
    if (!has_node(n)) {
      return false;
    }
    while (n.degree() > 0) {
      remove_edge(*(n.edge_begin()));
    }
    size_type i = node_uid_to_node_idx_map[n.uid];
    node_uid_to_node_idx_map[node_info_list.back()->uid] = i;
    node_info_list[i] = std::move(node_info_list.back());
    node_uid_to_node_idx_map.erase(n.uid);
    node_info_list.pop_back();
    return true;
  }

  /* Remove the node by iterator.
   */
  bool remove_node(const node_iterator& iter) {
    return remove_node(*iter);
  }

  /* 
   * Remove edge connecting @a n1 and @a n2.
   * @return true if the edge is in the old graph, false otherwise.
   * @post new has_edge(@a n1, @a n2) == false.
   * @post new num_edges() == old num_edges() - int(result)
   * Can invalidate edge indices and iterators of edges coming after
   *     the deleted edge in the old order.
   * Complexity: O(log(old size())).
   */
  bool remove_edge(const Node& n1, const Node& n2) {
    if (!(has_node(n1) && has_node(n2) && has_edge(n1, n2))) {
      return false;
    }
    node_info_list[n1.index()]->neighbor_uid_set.erase(n2.uid);
    node_info_list[n2.index()]->neighbor_uid_set.erase(n1.uid);
    uid_type uid1 = std::min(n1.uid, n2.uid);
    uid_type uid2 = std::max(n1.uid, n2.uid);
    auto edge_uid = std::make_pair(uid1, uid2);
    auto edge_uid_back = std::make_pair(
        edge_info_list.back()->uid1, edge_info_list.back()->uid2);
    size_type i = edge_uid_to_edge_idx_map[edge_uid];
    edge_uid_to_edge_idx_map[edge_uid_back] = i;
    edge_info_list[i] = std::move(edge_info_list.back());
    edge_uid_to_edge_idx_map.erase(edge_uid);
    edge_info_list.pop_back();
    return true;
  }

  /* 
   * Remove edge by giving edge itself.
   */
  bool remove_edge(const Edge& edge) {
    return remove_edge(edge.node1(), edge.node2());
  }

  /* 
   * Remove edge by iterator.
   */
  bool remove_edge(const edge_iterator& iter) {
    return remove_edge(*iter);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct NodeInfo {
    uid_type uid;
    Point position;
    node_value_type value;
    std::set<uid_type> neighbor_uid_set;
    NodeInfo(const uid_type uid, const Point& position,
        const node_value_type& value = node_value_type()) :
        uid(uid), position(position), value(value) {}
  };

  struct EdgeInfo {
    uid_type uid1;
    uid_type uid2;
    edge_value_type value;
    EdgeInfo(const uid_type uid1, const uid_type uid2,
        const edge_value_type& value = edge_value_type()) :
        uid1(uid1), uid2(uid2), value(value) {}
  };

  uid_type next_node_uid_to_use = 0;

  std::vector<std::unique_ptr<NodeInfo>> node_info_list;
  std::map<uid_type, size_type> node_uid_to_node_idx_map;
  
  std::vector<std::unique_ptr<EdgeInfo>> edge_info_list;
  std::map<std::pair<uid_type, uid_type>, size_type> edge_uid_to_edge_idx_map;

};

#endif // CME212_GRAPH_HPP
