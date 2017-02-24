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
template <typename V, typename E>
class Graph {

 private:

  // HW0: YOUR CODE HERE

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** HW1: Allow Node to support user-specified data */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** HW2: Allow Edge to support user-specified data */
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

  /** Type of node info struct, which stores the node's data. */
  struct node_info;
  using node_info_type = node_info;
  /** Type of edge info struct, which stores the edge's data. */
  struct edge_info;
  using edge_info_type = edge_info;

  /** Type of adjacency list, which stores the uids of the nodes
      that is adjacent to node */
  using adj_list_type = std::vector<edge_info_type>;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
      : nodes_(), idx2uid_() {
    // HW0: YOUR CODE HERE
    /** No additional code needed */
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODE INFO
  //

  /** @struct Graph::node_info *
   * @brief Struct representing the node's data.
   *
   * Node info objects are used to store information about the node.
   */
  struct node_info {
      size_type idx_;
      Point p_;
      node_value_type v_;
      adj_list_type adj_;

     /** Constructor */
     node_info(size_type idx, Point p, node_value_type v, adj_list_type adj) :
       idx_(idx), p_(p), v_(v), adj_(adj) {
     }
  };

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

    /** HW2: Return this node's (modifiable) position */
    Point& position() {
        return graph_->nodes_[uid_].p_;
    }

    /** Return this node's position. */
    const Point& position() const {
        // HW0: YOUR CODE HERE
        // HW1: Modified
        return graph_->nodes_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].idx_;
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
        return graph_->nodes_[uid_].v_;
    }

    /** HW1: Return this node's value */
    const node_value_type& value() const {
        return graph_->nodes_[uid_].v_;
    }

    /** HW1: Return this node's degree (number of incident edges) */
    size_type degree() const {
        return graph_->nodes_[uid_].adj_.size();
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
    return idx2uid_.size();
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
  node_type add_node(const Point& position,
                     const node_value_type& value = node_value_type()) {

      /** Create struct to store node's data */
      adj_list_type adj;
      node_info_type ni(num_nodes(), position, value, adj);

      /** Add to nodes vector */
      nodes_.push_back(ni);
      /** Add to vector of uids */
      idx2uid_.push_back( nodes_.size() - 1 );

      return Node(this, nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if ( n.uid_ >= nodes_.size() ) {
        return false;
    }
    return ( (n.graph_ == this) && ( (int)nodes_[n.uid_].idx_ != -1 ) &&
             (idx2uid_[ nodes_[n.uid_].idx_ ] == n.uid_ ) );
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, idx2uid_[i]);
  }

  /** Remove a node from the graph, or return 0 to indicate
   *  node @a n is not graph @a g.
   * @pre      Graph g != null
   * @return   0 if has_node(@a n) == false; 1 if has_node(@a n) == true
   * @post     has_node(@a n) == false
   * @post     If old graph !has_node(@a n), new num_nodes() == old num_nodes().
   *           Else,                         new num_nodes() == old num_nodes() - 1.
   * @post     For i == n.index(), new index(i) == old index( num_nodes() - 1 )
               For i != n.index(), 0 <= i < num_nodes() - 1, new index(i) = old index(i)
   * @post     All edges adjacent to @a n are removed:
               If old has_edge(@a n.uid_, b.uid_) == true or
                  old has_edge(b.uid_, @a n.uid_) == true,
                  new has_edge(@a n.uid_, b.uid_) == false and
                  new has_edge(b.uid_, @a n.uid_) == false.
   * Complexity of O( max degree()^2 )
   */
  size_type remove_node(const Node& n) {

      /* First check that the node exists */
      if ( !has_node(n) ) {
          return 0;
      }

      /* Remove any edges associated with this node */
      if ( n.degree() > 0 ) {
          incident_iterator start = n.edge_begin();
          while( start != n.edge_end() ) {
              start = remove_edge(start);
          }

          /*for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
              remove_edge( n, Node(this, (*ii).n2_uid_) );
          }*/
      }

      /* Get index of node to remove */
      size_type idx = n.index();

      if ( idx != (num_nodes() - 1) ) {
          /* Copy last node to the spot where the node to
             be removed is stored */
          idx2uid_[idx] = idx2uid_.back();
          /* Update the replacement node's index */
          nodes_[ idx2uid_.back() ].idx_ = idx;
      }

      /* Remove last node, which is either the node to be removed itself,
         or a duplicate of the node that was moved */
      idx2uid_.pop_back();
      /* Chnge removed node's index to -1 */
      nodes_[ n.uid_ ].idx_ = -1;

      return 1;
  }

  /** Remove a node from the graph.
   * @pre      Graph g != null
   * @pre      g.node_begin() <= n_it < g.node_end()
   * @return   n_it
   * @post     has_node(@a n) == false
   * @post     If old graph !has_node(@a n), new num_nodes() == old num_nodes() and
                                             new n_it = old n_it.
   *           Else,                         new num_nodes() == old num_nodes() - 1 and
                                             new n_it = @a *n.
   * @post     For i == n.index(), new index(i) == old index( num_nodes() - 1 )
               For i != n.index(), 0 <= i < num_nodes() - 1, new index(i) = old index(i)
   * @post     All edges adjacent to @a n are removed:
               If old has_edge(@a n.uid_, b.uid_) == true or
                  old has_edge(b.uid_, @a n.uid_) == true,
                  new has_edge(@a n.uid_, b.uid_) == false and
                  new has_edge(b.uid_, @a n.uid_) == false.
   * Complexity of O( max degree()^2 )
   */
  node_iterator remove_node(node_iterator n_it) {
      remove_node(*n_it);
      return n_it;
  }


  //
  // EDGE INFO
  //

  /** @struct Graph::edge_info *
   * @brief Struct representing the edge's data.
   *
   * Edge info objects are used to store information about the edge.
   */

  struct edge_info {
      size_type n2_uid_;
      edge_value_type value_;

      /** Constructor */
      edge_info(size_type n2_uid, edge_value_type value) :
        n2_uid_(n2_uid), value_(value) {
      }
  };

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
      return Node(graph_, n2_uid_);
    }

    /** Set initial edge length */
    double length() const {
        return norm( node1().position() - node2().position() );
    }

    /** HW2: Return this edge's data */
    edge_value_type& value() {
        return graph_->nodes_[n1_uid_].adj_[fetch_n2_aux_id()].value_;
    }

    /** HW1: Return this edge's value */
    const edge_value_type& value() const {
        return graph_->nodes_[n1_uid_].adj_[fetch_n2_aux_id()].value_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ( (graph_ == e.graph_) &&
               ( std::min(n1_uid_, n2_uid_) ==
                 std::min(e.node1().uid_, e.node2().uid_) ) &&
               ( std::max(n1_uid_, n2_uid_) ==
                 std::max(e.node1().uid_, e.node2().uid_) ) );
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        size_type min_uid = std::min(n1_uid_, n2_uid_);
        size_type max_uid = std::max(n1_uid_, n2_uid_);
        size_type e_min_uid = std::min(e.n1_uid_, e.n2_uid_);
        size_type e_max_uid = std::max(e.n1_uid_, e.n2_uid_);
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

    /** Pointer back to the graph */
    graph_type* graph_;

    /** Index of node 1 */
    size_type n1_uid_;

    /** Index of node 2 */
    size_type n2_uid_;

    /** HW1: Private constructor */
    /* @pre @a graph != nullptr */
    /* @pre @a n1_uid <= graph_.num_nodes() */
    /* @pre @a n2_aux_id <= graph_.adjacency_matrix_[n1_uid].size() */
    Edge(const graph_type* graph, size_type n1_uid, size_type n2_uid) :
         graph_(const_cast<graph_type*>(graph)), n1_uid_(n1_uid), n2_uid_(n2_uid) {
    }

    /** Get index of node2 in node1's adjacency list */
    size_type fetch_n2_aux_id() const {
        for (size_type i = 0; i < graph_->nodes_[n1_uid_].adj_.size(); ++i) {
            if (graph_->nodes_[n1_uid_].adj_[i].n2_uid_ == n2_uid_)
                return i;
        }
        return graph_->nodes_[n1_uid_].adj_.size();
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

    /* First check if nodes a and b are in graph */
    if ( !has_node(a) || !has_node(b) ) {
        return false;
    }

    for (size_type ei = 0; ei < a.degree(); ++ei) {
        if ( b.uid_ == nodes_[a.uid_].adj_[ei].n2_uid_ ) {
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
  edge_type add_edge(const Node& a, const Node& b,
                     const edge_value_type& data = edge_value_type()) {
    // HW0/HW1
    if ( has_edge(a, b) ) {
        return Edge(this, a.uid_, b.uid_);
    }

    /** Create edge (a,b) and edge (b,a) along with associated data */
    edge_info_type ei_a(b.uid_, data);
    edge_info_type ei_b(a.uid_, data);

    /** Add edges to node a's and node b's adjacency lists */
    nodes_[a.uid_].adj_.push_back(ei_a);
    nodes_[b.uid_].adj_.push_back(ei_b);

    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove an edge (@a a, @a b) from the graph, or return 0 to indicate
   *  edge (@a a, @a b) is not graph g.
   * @pre      Graph g != null
   * @return   0 if has_edge(@a a, @a b) == false; 1 if has_edge(@a a, @a b) == true
   * @post     has_edge(@a a, @a b) == false and has_edge(@a b, @a a) == false
   * @post     If old graph !has_edge(@a a, @a b), new num_edges() == old num_edges().
   *           Else,                               new num_edges() == old num_edges() - 1.
   * @post     Symmetric edges are removed:
   *           If old has_edge(@a a, @a b) == true or old has_edge(@a b, @a a) == true,
   *              new has_edge(@a a, @a b) == false and new has_edge(@a b, @a a) == false.
   * Complexity of O( max( a.degree(), b.degree() ) )
   */
  size_type remove_edge(const Node& a, const Node& b) {

      /* Check that the edge exists */
      if ( !has_edge(a,b) ) {
          return 0;
      }

      /* Find where in node a's adjacency list edge (a, b) is stored */
      for (size_type i = 0; i < a.degree(); ++i) {
          if ( b.uid_ == nodes_[a.uid_].adj_[i].n2_uid_ ) {
              /* Replace edge to be deleted with last edge stored
                 in node a's adjacency list */
              nodes_[a.uid_].adj_[i] = nodes_[a.uid_].adj_.back();
          }
      }
      /* Remove last edge stored in node a's adjacency list,
         which is now a duplicate */
      nodes_[a.uid_].adj_.pop_back();

     /* Do the same thing for edge (b, a) */
      for (size_type i = 0; i < b.degree(); ++i) {
          if ( a.uid_ == nodes_[b.uid_].adj_[i].n2_uid_ ) {
              nodes_[b.uid_].adj_[i] = nodes_[b.uid_].adj_.back();
          }
      }
      /* Remove last edge stored in node b's adjacency list,
         which is now a duplicate */
      nodes_[b.uid_].adj_.pop_back();

     return 1;

  }


  /** Remove an edge @a e from the graph, or return 0 to indicate
   *  edge @a e is not graph g.
   * @pre      Graph g != null
   * @return   0 if has_edge(@a e.node1(), @a e.node2()) == false;
   *           1 if has_edge(@a e.node1(), @a e.node2()) == true
   * @post     has_edge(@a e.node1(), @a e.node2()) == false and
   *           has_edge(@a e.node2(), @a e.node1()) == false
   * @post     If old graph !has_edge(@a e.node1(), @a e.node2()),
   *                 new num_edges() == old num_edges().
   *           Else, new num_edges() == old num_edges() - 1.
   * @post     Symmetric edges are removed:
   *           If old has_edge(@a e.node1(), @a e.node2()) == true or
   *              old has_edge(@a e.node2(), @a e.node1()) == true,
   *              new has_edge(@a e.node1(), @a e.node2()) == false and
   *              new has_edge(@a e.node2(), @a e.node1()) == false.
   * Complexity of O( max( e.node1().degree(), e.node2().degree() ) )
   */
  size_type remove_edge(const Edge& e) {
      return remove_edge( Node(this, e.n1_uid_), Node(this, e.n2_uid_) );
  }

  /** Remove an edge from the graph, or return 0 to indicate
   *  the edge is not graph g.
   * @pre      Graph g != null
   * @pre      g.edge_begin() <= e_it < g.edge_end()
   * @return   e_it
   * @post     has_edge(@a (*e_it).node1(), @a (*e_it).node2()) == false and
               has_edge(@a (*e_it).node2(), @a (*e_it).node1()) == false
   * @post     If old graph !has_edge(@a (*e_it).node1(), @a (*e_it).node2()),
                     new num_edges() == old num_edges().
   *           Else, new num_edges() == old num_edges() - 1.
   * @post     Symmetric edges are removed:
               If old has_edge(@a (*e_it).node1(), @a (*e_it).node2()) == true or 
                  old has_edge(@a (*e_it).node2(), @a (*e_it).node1()) == true,
                  new has_edge(@a (*e_it).node1(), @a (*e_it).node2()) == false and
                  new has_edge(@a (*e_it).node2(), @a (*e_it).node1()) == false.
   * Complexity of O( max( (*e_it).node1().degree(), (*e_it).node2().degree() ) )
   */
  edge_iterator remove_edge(edge_iterator e_it) {
      remove_edge(*e_it);
      return e_it;
  }

  /** Remove an edge incident on node n from the graph, or return 0 to indicate
   *  the edge is not graph g.
   * @pre      Graph g != null
   * @pre      n.edge_begin() <= i_it < n.edge_end()
   * @return   i_it
   * @post     has_edge(@a (*i_it).node1(), @a (*i_it).node2()) == false and
               has_edge(@a (*i_it).node2(), @a (*i_it).node1()) == false
   * @post     If old graph !has_edge(@a (*i_it).node1(), @a (*i_it).node2()),
                     new num_edges() == old num_edges().
   *           Else, new num_edges() == old num_edges() - 1.
   * @post     Symmetric edges are removed:
               If old has_edge(@a (*i_it).node1(), @a (*i_it).node2()) == true or 
                  old has_edge(@a (*i_it).node2(), @a (*i_it).node1()) == true,
                  new has_edge(@a (*i_it).node1(), @a (*i_it).node2()) == false and
                  new has_edge(@a (*i_it).node2(), @a (*i_it).node1()) == false.
   * Complexity of O( max( (*i_it).node1().degree(), (*i_it).node2().degree() ) )
   */
  incident_iterator remove_edge(incident_iterator i_it) {
      remove_edge(*i_it);
      return i_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    idx2uid_.clear();
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
        /* Get node's uid */
        size_type uid = graph_->idx2uid_[idx_];
        return Node(graph_, uid);
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
        /* Get node1's uid */
        size_type n1_uid = graph_->idx2uid_[idx_];
        size_type n2_uid = graph_->nodes_[n1_uid].adj_[idx_aux_].n2_uid_;
        return Edge(graph_, n1_uid, n2_uid);
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
        /* Get node1's uid */
        size_type n1_uid = graph_->idx2uid_[idx_];
        size_type n2_uid = graph_->nodes_[ graph_->idx2uid_[idx_] ].adj_[idx_aux_].n2_uid_;
        return Edge(graph_, n1_uid, n2_uid);
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
            if ( idx_ >= graph_->num_nodes() ) {
                /** No more edges to iterate over */
                break;
            }
            else if ( idx_aux_ == graph_->nodes_[ graph_->idx2uid_[idx_] ].adj_.size() ) {
                /** No more edges to iterate over for this node; move to next one */
                idx_aux_ = 0;
                ++idx_;
            }
            else if ( graph_->idx2uid_[idx_] >
                      graph_->nodes_[ graph_->idx2uid_[idx_] ].adj_[idx_aux_].n2_uid_ ) {
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
  std::vector<node_info_type> nodes_;

  /** HW2: Vector of node uids */
  std::vector<size_type> idx2uid_;

};

#endif // CME212_GRAPH_HPP
