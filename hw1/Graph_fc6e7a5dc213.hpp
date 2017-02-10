#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

/* HW0/HW1 */
#include <utility>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** HW1: Allow Node to support user-specified value */
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
      : nodes_(), adjacency_matrix_() {
    // HW0: YOUR CODE HERE
    /** No additional code needed */
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
      /** No additional code needed */
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // HW1: Modified
      return graph_->nodes_[uid_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** HW1: Return this node's value */
    node_value_type& value() {
        return graph_->nodes_[uid_].second;
    }

    /** HW1: Return this node's value */
    const node_value_type& value() const {
        return graph_->nodes_[uid_].second;
    }

    /** HW1: Return this node's degree (number of incident edges) */
    size_type degree() const {
        return graph_->adjacency_matrix_[uid_].size();
    }

    /** HW1 */
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_, uid_, 0);
    }

    /** HW1 */
    incident_iterator edge_end() const {
        return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ( (graph_ == n.graph_) && (uid_ == n.uid_) );
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
      return ( (graph_ < n.graph_) ||
             ( (graph_ == n.graph_) && (uid_ < n.uid_) ) );
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    /** HW0: Pointer back to the Graph */
    graph_type* graph_;

    /** HW0: Node's unique ID number */
    size_type uid_;

    /** HW0: Private constructor */
    /* @pre @a graph != nullptr */
    /* @pre @a uid <= graph_.size() */
    Node(const graph_type* graph, size_type uid) :
         graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();
  }

  /** HW0: Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   *
   * HW1: Modified to take a node_value_type as an input */
  node_type add_node(const Point& position, const node_value_type& value = node_value_type()) {

      /** Create pair with position & value, and add to nodes vector */
      nodes_.push_back(std::make_pair(position, value));

      /** Add node to first level of adjacency matrix */
      adjacency_matrix_.push_back({});
      return Node(this, size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ( (n.graph_ == this) && (n.uid_ < size()) ) ;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
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
      /** No additional code needed */
    }

    /** HW0: Return node 1 of this Edge */
    node_type node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, n1_uid_);
    }

    /** HW1: Return node 2 of this Edge */
    node_type node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, fetch_n2_uid());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ( (graph_ == e.graph_) &&
               ( std::min(n1_uid_, fetch_n2_uid()) ==
                 std::min(e.node1().uid_, e.node2().uid_) ) &&
               ( std::max(n1_uid_, fetch_n2_uid()) ==
                 std::max(e.node1().uid_, e.node2().uid_) ) );
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        size_type min_uid = std::min(n1_uid_, fetch_n2_uid());
        size_type max_uid = std::max(n1_uid_, fetch_n2_uid());
        size_type e_min_uid = std::min(e.n1_uid_, e.fetch_n2_uid());
        size_type e_max_uid = std::max(e.n1_uid_, e.fetch_n2_uid());
        return ( (graph_ < e.graph_) ||
               ( (graph_ == e.graph_) && ( (min_uid < e_min_uid) ||
               ( (min_uid == e_min_uid) && (max_uid < e_max_uid) ) ) ) );
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    /** Pointer back to the Graph */
    graph_type* graph_;

    /** Index of node 1 of edge */
    size_type n1_uid_;

    /** Index of node 2's uid in the vector corresponding to n1 in the adjacency matrix */
    size_type n2_aux_id_;

    /** HW1: Private constructor */
    /* @pre @a graph != nullptr */
    /* @pre @a n1_uid <= graph_.num_nodes() */
    /* @pre @a n2_aux_id <= graph_.adjacency_matrix_[n1_uid].size() */
    Edge(const graph_type* graph, size_type n1_uid, size_type n2_aux_id) :
         graph_(const_cast<graph_type*>(graph)),
          n1_uid_(n1_uid), n2_aux_id_(n2_aux_id) {
    }

    /** Get node 2 uid */
    size_type fetch_n2_uid() const {
        return graph_->adjacency_matrix_[n1_uid_][n2_aux_id_];
    }


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return std::distance(edge_begin(), edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0/HW1
    for (size_type ei = 0; ei < adjacency_matrix_[a.uid_].size(); ++ei) {
        if (adjacency_matrix_[a.uid_][ei] == b.uid_)
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
  edge_type add_edge(const Node& a, const Node& b) {
    // HW0/HW1
    for (size_type ei = 0; ei < adjacency_matrix_[a.uid_].size(); ++ei) {
        if (adjacency_matrix_[a.uid_][ei] == b.uid_) {
            return Edge(this, a.uid_, ei);
        }
    }

    /** Modify adjacency matrix to include edge (a, b) and edge (b,a) */
    adjacency_matrix_[a.uid_].push_back(b.uid_);
    adjacency_matrix_[b.uid_].push_back(a.uid_);

    return Edge(this, a.uid_, adjacency_matrix_[a.uid_].size() - 1);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    adjacency_matrix_.clear();
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
        /** No additional code needed */
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** HW1: Node dereference operator */
    value_type operator*() const {
        return Node(graph_, idx_);
    }

    /** HW1: Node increment operator*/
    node_iterator& operator++() {
        ++idx_;
        return *this;
    }

    /** HW1: Node iterator equality operator */
    bool operator==(const node_iterator& ni) const{
       return ( (graph_ == ni.graph_) && (idx_ == ni.idx_) );
    }

   private:
    friend class Graph;

    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;

    size_type idx_;

    /** HW1: Private constructor */
    /* @pre @a graph != nullptr */
    /* @pre @a idx <= graph_.num_nodes() */
    NodeIterator(const graph_type* graph, size_type idx) :
      graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** HW1 */
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }

  /** HW1 */
  node_iterator node_end() const {
     return NodeIterator(this, num_nodes());
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
        /** No additional code needed */
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** HW1: Incident iterator dereference operator */
    value_type operator*() const {
        return Edge(graph_, idx_, idx_aux_);
    }

    /** HW1: Incident iterator increment operator */
    incident_iterator& operator++() {
        ++idx_aux_;
        return *this;
    }

    /** HW1: Incident iterator equality operator */
    bool operator==(const incident_iterator& ii) const {
        return ( (graph_ == ii.graph_) &&
                 (idx_ == ii.idx_) && (idx_aux_ == ii.idx_aux_) );
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    graph_type* graph_;

    size_type idx_;

    size_type idx_aux_;

    /** HW1: Private constructor */
    /* @pre @a graph != nullptr */
    /* @pre @a idx <= graph_.num_nodes() */
    /* @pre @a idx_aux <= graph_.adjacency_matrix[idx].size() */
    IncidentIterator(const graph_type* graph, size_type idx, size_type idx_aux) :
      graph_(const_cast<graph_type*>(graph)), idx_(idx), idx_aux_(idx_aux) {
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
        /** No additional code needed */
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** HW1: Edge iterator dereference operator */
    value_type operator*() const {
       return Edge(graph_, idx_, idx_aux_);
    }

    /** HW1: Edge iterator increment operator */
    edge_iterator& operator++() {
        ++idx_aux_;
        find_next_edge();
        return *this;
    }

    /** HW1: Edge iterator equality operator */
    bool operator==(const edge_iterator& eit) const{
       return ( (graph_ == eit.graph_) &&
                (idx_ == eit.idx_) && (idx_aux_ == eit.idx_aux_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    graph_type* graph_;

    size_type idx_;

    size_type idx_aux_;

    /** HW1: Private constructor */
    /* @pre @a graph != nullptr */
    /* @pre @a idx <= graph_.num_nodes() */
    /* @pre @a idx_aux <= graph_.adjacency_matrix_[idx].size() */
    EdgeIterator(const graph_type* graph, size_type idx, size_type idx_aux) :
      graph_(const_cast<graph_type*>(graph)), idx_(idx), idx_aux_(idx_aux) {
        find_next_edge();
    }

    /** HW1: Helper function to find next edge to iterate over */
    void find_next_edge() {
        while ( true ) {
            if ( idx_ == graph_->num_nodes() ) {
                /** No more edges to iterate over */
                break;
            }
            else if ( idx_aux_ == graph_->adjacency_matrix_[idx_].size() ) {
                /** No more edges to iterate over for this node; move to next one */
                idx_aux_ = 0;
                ++idx_;
            }
            else if ( idx_ > graph_->adjacency_matrix_[idx_][idx_aux_] ) {
                /** Only iterate over edges (a, b) where a < b in order to
                avoid iterating over the same edge twice */
                ++idx_aux_;
            }
            else {
                break;
            }
        }
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** HW1 */
  edge_iterator edge_begin() const {
      return EdgeIterator(this, 0, 0);
  }

  /** HW1 */
  edge_iterator edge_end() const {
      return EdgeIterator(this, num_nodes(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** HW1: Vector of pairs to store node position & value */
  std::vector<std::pair<Point, node_value_type>> nodes_;

  /** HW1: Vector of vectors to store indices of adjacent nodes */
  std::vector<std::vector<size_type>> adjacency_matrix_;

};

#endif // CME212_GRAPH_HPP
