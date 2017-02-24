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
  void empty() {
    head_ = nullptr;
    tail_ = nullptr;
    size = 0;
    /* need to destroy the elements */
  }
  bool is_empty() {
    return head_ == nullptr;
  }
  void erase(T& t) {
    auto p = head_;
    if (is_empty()) {
      return;
    }
    if (p->val_ == t) {
      head_ = p->next_;
      size--;
    } else {
      for(; p->next_ != nullptr; p = p->next_) {
        if (p->next_->val_ == t) {
          p->next_ = p->next_->next_;
          size--;
          break;
        }
      }
    }
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
template <typename V, typename E>
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
  using graph_type = Graph<V, E>;

  /** Predeclaration of Node type. */
  using node_value_type = V;
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  using edge_value_type = E;
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
    Point& position() const {
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
      return graph_->incident_edges[index_].size;
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

    size_type index() const {
      return index_;
    }

    double length() const {
      return norm(node1().position() - node2().position());
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

    edge_value_type& value() {
      return graph_->edge_values[index_];
    }

    const edge_value_type& value() const {
      return graph_->edge_values[index_];
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
    for (size_type i = 0; i < num_edges(); i++) {
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
  Edge add_edge(const Node& a, const Node& b, 
      const edge_value_type& v = edge_value_type()) {
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
    edge_values.push_back(v);
    Edge e1(this, node1_index.size() - 1, false);
    Edge e2(this, node1_index.size() - 1, true);
    incident_edges[a.index()].push_back(e1);
    incident_edges[b.index()].push_back(e2);
    return e1;        
  }

  /** Remove one element from a vector with O(1)
   * @param[v]    vector to be erased
   * @param[idx]  index of the element to be erased
   *
   * @post        v old[0:idx] = v new[0:idx]
   * @post        v old[idx+1:v old.size() - 1] = v new[idx+1:v new.size()]
   * @post        v old[v old.size() - 1] = v new[idx]
   *
   * Complexity O(1)
   */
  template <typename T>
  void erase(std::vector<T>& v, size_type idx) {
    assert(idx < v.size());
    v[idx] = v.back();
    v.pop_back();
  }
  /** Remove edge
   * @param[e0]   the edge to be removed
   *
   * @post        edge e0 will be removed from graph this
   *
   * Complexity O(d * 2) + O(1 * 3) = O(d), 
   * d = max(n.degree) where n is any nodes of this graph
   */
  size_type remove_edge(const Edge& e0) {
    Edge e;
    e = Edge(this, e0.index(), false);
    incident_edges[e.node1().index()].erase(e);
    e = Edge(this, e0.index(), true);
    incident_edges[e.node1().index()].erase(e);

    erase(node1_index, e.index());
    erase(node2_index, e.index());
    erase(edge_values, e.index());
    return num_edges();
  }
  size_type remove_edge(const Node& a, const Node& b) {
    Edge e0 = add_edge(a, b);
    return remove_edge(e0);
  }
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type idx = (*e_it).index();
    Edge e = Edge(this, idx, false);
    remove_edge(e);
    return EdgeIterator(this, idx);
  }
  /** Remove node
   * @param[n]   the node to be removed
   *
   * @post        edge node will be removed from graph this
   *
   * Complextity = O(d^2)
   * which is O(d)[for the number of edges] * O(d)[to remove each edge] + O(1) to remove node
   * d = max(n.degree) where n is any nodes of this graph
   */
  size_type remove_node(const Node& n) {
    Edge e, e_mov;
    size_type dg = n.degree(), i = 0, idx = n.index();
    for (incident_iterator iit = n.edge_begin(); i < dg; ++iit, ++i) {
      e = edge((*iit).index());
      if (e.node1() == n) {
        e.inverse = !e.inverse;
      }
      incident_edges[(*iit).node2().index()].erase(e);

      erase(node1_index, e.index());
      erase(node2_index, e.index());
      erase(edge_values, e.index());

      // change the corresponding edge stored in incident_edges
      Edge e_mov = edge(e.index());
      for (incident_iterator iit = e_mov.node1().edge_begin(); iit != e_mov.node1().edge_end(); ++iit) {
        if (*iit == e_mov) {
          iit.eli.current_position->val_.index_ = e_mov.index_;
          break;
        }
      }
      e_mov.inverse = !e_mov.inverse;
      for (incident_iterator iit = e_mov.node1().edge_begin(); iit != e_mov.node1().edge_end(); ++iit) {
        if (*iit == e_mov) {
          iit.eli.current_position->val_.index_ = e_mov.index_;
          break;
        }
      }
    }

    erase(node_values, idx);
    erase(node_positions, idx);

    incident_edges[idx].empty();
    erase(incident_edges, idx);
    for (incident_iterator iit = node(idx).edge_begin(); iit != node(idx).edge_end(); ++iit) {
      if ((*iit).inverse == false) {
        node1_index[(*iit).index()] = idx;
      } else {
        node2_index[(*iit).index()] = idx;
      }
    }
    return idx;
  }
  node_iterator remove_node(node_iterator n_it) {
    size_type idx = (*n_it).index();
    remove_node(*n_it);
    return NodeIterator(this, idx);
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
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> node_positions;
  std::vector<node_value_type> node_values;

  std::vector<size_type> node1_index;
  std::vector<size_type> node2_index;
  std::vector<edge_value_type> edge_values;
};

#endif // CME212_GRAPH_HPP
