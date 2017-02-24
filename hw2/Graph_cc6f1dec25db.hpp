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
  /** Predeclarations. */
  struct node_data;
  struct edge_data;
  struct pair_data;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph<V>;

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
  Graph()
    // HW0: YOUR CODE HERE
    : nodes_(), edges_(), node_count(), edge_count() {
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return this node's value. */
    V& value() {
      return graph_->nodes_[uid_].type;
    }
    /** Return this node's value. */
    const V& value() const {
      return graph_->nodes_[uid_].type;
    }
    /** Return this node's degree. */
    size_type degree() const {
      return(graph_->nodes_[uid_].num_neighbors);
    }
    /** Return IncidentIterator. */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, graph_->nodes_[uid_].neighbors.begin());
    }
    /** Return IncidentIterator. */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, graph_->nodes_[uid_].neighbors.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_ && uid_ == n.uid_);
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
      // HW0: YOUR CODE HERE
      return ((graph_ < n.graph_) || ((graph_ == n.graph_) && (uid_ < n.uid_)));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Allow Node to access parent Graph object
    Graph* graph_;
    // Unique ID
    size_type uid_;
    // Private Constructor
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) { 
        assert(graph_ != nullptr);
        assert(uid_ < graph_->num_nodes()); 
     }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_count;
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
  Node add_node(const Point& position, const V& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Add a set pair_data
    std::set<pair_data> temp_set;
    nodes_.push_back({position, node_count, value, 0, temp_set});
    ++node_count;
    // Return Node
    return Node(this, node_count-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (nodes_[n.uid_].position == n.position());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < node_count);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->edges_[eid_].n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->edges_[eid_].n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE
      return((graph_ == e.graph_) && (eid_ == e.eid_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE
      return ((graph_ < e.graph_) || ((graph_ == e.graph_) && (eid_ < e.eid_)));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Allow Edge to access parent Graph object
    Graph* graph_;
    // Edge ID
    size_type eid_;
    // Private constructor
    Edge(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), eid_(uid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < edge_count);
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_ == this && b.graph_ == this);
    node_data a_inf = nodes_[a.uid_];
    node_data b_inf = nodes_[b.uid_];
    pair_data b_part = {0, b.uid_};
    typename std::set<pair_data>::iterator it = a_inf.neighbors.find(b_part);
    return it != a_inf.neighbors.end();
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
    // HW0: YOUR CODE HERE
    assert(a.uid_ != b.uid_);
    node_data a_inf = nodes_[a.uid_];
    node_data b_inf = nodes_[b.uid_];

    // 0 for comparison
    pair_data b_part = {0, b.uid_};

    // a's neighbors and b's ID
    typename std::set<pair_data>::iterator it = a_inf.neighbors.find(b_part);

    // Create new edge if it cannot find
    if(it == a_inf.neighbors.end()) {
      // Add new edge
      edge_data temp = {a, b};
      edges_.push_back(temp);

      // Store a and b
      nodes_[a.uid_].neighbors.insert({edge_count, b.uid_});
      nodes_[b.uid_].neighbors.insert({edge_count, a.uid_});

      // Increment count
      ++nodes_[a.uid_].num_neighbors;
      ++nodes_[b.uid_].num_neighbors;
      ++edge_count;

      // Return Edge
      return Edge(this, edge_count-1);
    } else {
      return Edge(this, (*it).edge_pair_id);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    edge_count = 0;
    node_count = 0;
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
    /** Return this node's position. */
    Node operator*() const {
      assert(node_position < graph->size());
      return Node(graph, node_position);
    }
    /** Return updated value. */
    NodeIterator& operator++(){
      if(node_position < graph->size()) {
        ++node_position;
      }
      return *this;
    }
    /** Return comparison. */
    bool operator==(const NodeIterator& n_it) const {
      return(graph == n_it.graph && node_position == n_it.node_position);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    // Points to graph
    const Graph* graph;
    // Node's position
    size_type node_position;
    // Private Constructor
    NodeIterator(const Graph* graph, size_type uid)
      : graph(const_cast<Graph*>(graph)), node_position(uid) {  
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return node iterator. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /** Return size of graph. */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
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
    /** Return incident iterator. */
    Edge operator*() const {
      size_type id = (*incident_position).edge_pair_id;
      return Edge(graph, id);
    }
    /** Return updated value. */
    IncidentIterator& operator++() {
      ++incident_position;
      return *this;
    }
    /** Return comparison. */
    bool operator==(const IncidentIterator& i_it) const {
      return(graph == i_it.graph && incident_position == i_it.incident_position);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Points to graph
    Graph* graph;
    // Set of pair_data
    typename std::set<pair_data>::iterator incident_position;
    // Private Constructor
    IncidentIterator(const Graph* graph, 
      typename std::set<pair_data>::iterator it)
        : graph(const_cast<Graph*>(graph)), incident_position(it) {}
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return proxy node. */
    Edge operator*() const {
      assert(edge_position < graph -> num_edges());
      return Edge(graph, edge_position);
    }
    /** Return updated value. */
    EdgeIterator& operator++() {
      if(edge_position < graph -> num_edges())
        ++edge_position;
      return *this;
    }
    /** Return comparison. */
    bool operator==(const EdgeIterator& e_id) const {
      return (graph == e_id.graph && edge_position == e_id.edge_position);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Points to graph
    Graph* graph;
    // Edge position
    size_type edge_position;
    // Private Constructor
    EdgeIterator(const Graph* graph, size_type uid)
      : graph(const_cast<Graph*>(graph)), edge_position(uid) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return iterator to first edge. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  /** Return iterator to memory. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  struct pair_data {
    size_type edge_pair_id;
    size_type partner_id;
     
    bool operator<(const pair_data& a) const {
      return(partner_id < a.partner_id);
    }
     
    bool operator==(const pair_data& a) const {
      return(partner_id == a.partner_id);
    }
  };   

  // Node's position with neighbors, type, and degree
  struct node_data {
    Point position;
    size_type uid;
    node_value_type type;
    size_type num_neighbors;
    std::set<pair_data> neighbors;
  };
   
  // Edge's two nodes
  struct edge_data{
    node_type n1;
    node_type n2;
  };

  // Vectors that node_data and edge_data
  std::vector<node_data> nodes_;
  std::vector<edge_data> edges_;

 // Count for nodes and edges
  size_type node_count;
  size_type edge_count;
};

#endif // CME212_GRAPH_HPP
