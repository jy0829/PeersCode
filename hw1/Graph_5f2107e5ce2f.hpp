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


template<typename T>
struct list_element {
  T val_;
  list_element <T>* next_;
  list_element(T v) : val_(v), next_(nullptr) {}
};
template<typename T>
struct list {
  list_element<T>* head_;
  list_element<T>* tail_;
  unsigned size;
  list() : head_(nullptr), tail_(nullptr), size(0) {};
  void push_back(T& lev) {
    list_element<T>* current = new list_element<T>(lev);
    if (head_ == nullptr) {
      head_ = current;
      tail_ = current;
    } else {
      tail_->next_ = current;
      tail_ = current;
    }
    size++;
  }
  T pop_head() {
    assert(head_ != nullptr);
    T ret = head_->val_;
    list_element<T>* tmp = head_;
    head_ = head_->next_;
    delete tmp;
    return ret;
  }
  bool is_empty() {
    return head_ == nullptr;
  }
};
template<typename T>
struct list_iterator {
  list_iterator(list_element<T>& le) : current_position(&le) {}
  list_element<T>* current_position;
  bool end() {
    return current_position == nullptr;
  }
  T operator*() const {
    return current_position->val_;
  }
  list_iterator& operator++() {
    current_position = current_position->next_;
    return *this;
  }
  bool operator==(const list_iterator& li) const {
    return current_position == li.current_position;
  }
};

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

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  using node_value_type = V;
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

  using edge_list = list<Edge>;
  using edge_list_element = list_element<Edge>;
  using edge_list_iterator = list_iterator<Edge>;
  std::vector<list<Edge>> incident_edges;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    num_edges_ = 0;
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
      return graph_->node_positions[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    size_type degree() const {
      return incident_edges.size;
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_->incident_edges[index_].head_);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_->incident_edges[index_].tail_->next_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph_ == n.graph_ && index_ == n.index_;
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
      //const Point& p0 = position();
      //const Point& p = n.position();
      //return p0.x + p0.y + p0.z < p.x + p.y + p.z;
      if (graph_ < n.graph_) {
        return true;
      } else if (graph_ > n.graph_) {
        return false;
      } else {
        return index_ < n.index_;
      }
    }

    node_value_type& value() {
      return graph_->node_values[index_];
    }

    const node_value_type& value() const {
      return graph_->node_values[index_];
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type index_;
    Node(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_positions.size();
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
  Node add_node(const Point& position, const node_value_type& v=node_value_type()) {
    // HW0: YOUR CODE HERE
    node_positions.push_back(position);
    node_values.push_back(v);
    Node n(this, node_positions.size() - 1);        
    incident_edges.push_back(list<Edge>());
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
      return inverse ? graph_->node(graph_->node2_index[index_]) : 
        graph_->node(graph_->node1_index[index_]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return inverse ? graph_->node(graph_->node1_index[index_]) : 
        graph_->node(graph_->node2_index[index_]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return node1() == e.node1() && node2() == e.node2() && graph_ == e.graph_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ < e.graph_) {
        return true;
      } else if (graph_ > e.graph_) {
        return false;
      }
      size_type n1 = node1().index();
      size_type n2 = node2().index();
      size_type en1 = e.node1().index();
      size_type en2 = e.node2().index();
      if (n1 < en1) {
        return true;
      } else if (n1 > en1) {
        return false;
      } else if (n2 < en2) {
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type index_;
    bool inverse;
    Edge(const Graph* graph, size_type index, bool inverse=false)
      : graph_(const_cast<Graph*>(graph)), index_(index), inverse(inverse) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return node1_index.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i, false);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    /*
    for (size_type i = 0; i < num_edges_; i++) {
      if ((node1_index[i] == a.index() && node2_index[i] == b.index()) ||
          (node1_index[i] == b.index() && node2_index[i] == a.index())) {
        return true;
      }
    }
    */
    IncidentIterator iit = a.edge_begin();
    IncidentIterator iit_end = a.edge_end();
    for (;iit != iit_end; ++iit) {
      if ((*iit).node2() == b) {
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    IncidentIterator iit = a.edge_begin();
    Edge e;
    for (;!iit.end(); ++iit) {
      e = *iit;
      if (e.node2() == b) {
        e = Edge(this, e.index_, false);
        if (e.node1() == a) {
          return e;
        } else {
          return Edge(this, e.index_, true);
        }
      }
    }

    node1_index.push_back(a.index());
    node2_index.push_back(b.index());
    Edge e1(this, node1_index.size() - 1, false);
    Edge e2(this, node1_index.size() - 1, true);
    //edge_list_element ele1(e1);
    //edge_list_element ele2(e2);
    incident_edges[a.index()].push_back(e1);
    incident_edges[b.index()].push_back(e2);
    num_edges_++;
    return e1;        
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
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

    NodeIterator(const Graph* graph, size_type i) : graph_(const_cast<Graph*>(graph)), i(i) {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() {
      return Node(graph_, i);
    }
    Node operator*() const {
      return Node(graph_, i);
    }
    NodeIterator& operator++() {
      i++;
      return *this;
    }
    bool operator==(const NodeIterator& ni) const {
      return graph_ == ni.graph_ && i == ni.i;
    }
    size_type operator-(const NodeIterator& ni) const {
      assert(graph_ == ni.graph_);
      return i - ni.i;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type i;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
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
    IncidentIterator(edge_list_iterator eli) : eli(eli) {}
    IncidentIterator(edge_list_element* ele) : eli(edge_list_iterator(*ele)) {
    }

    bool end() {
      return eli.end();
    }
    Edge operator*() const {
      //return eli.current_position->val_;
      return *eli;
    }
    incident_iterator& operator++() {
      ++eli;
      return *this;
    }
    bool operator==(const incident_iterator& iit) const {
      return eli == iit.eli;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    edge_list_iterator eli;
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
    EdgeIterator(const Graph* graph, size_type i) : graph_(const_cast<Graph*>(graph)), i(i) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    EdgeIterator& operator++() {
      i++;
      return *this;
    }
    Edge operator*() {
      return Edge(graph_, i);
    }
    bool operator==(const EdgeIterator& ei) const {
      return graph_ == ei.graph_ && i == ei.i;
    }
    size_type operator-(const EdgeIterator& ei) const {
      assert(graph_ == ei.graph_);
      return i - ei.i;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type i;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges_);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> node_positions;
  std::vector<node_value_type> node_values;
  size_type num_edges_;

  /*
  struct internal_edge {
    unsigned node1_idx, node2_idx;
    internal_edge(id1, id2): node1_idx(id1), node2_idx(id2) {}
  };
  std::vector<internal_edge> edges;
  */
  std::vector<size_type> node1_index;
  std::vector<size_type> node2_index;
};

#endif // CME212_GRAPH_HPP
