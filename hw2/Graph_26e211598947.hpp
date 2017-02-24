#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <tuple>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V,typename E> class Graph {
 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  //using graph_type = Graph<V>;
  using graph_type = Graph<V,E>;

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

  /** Synonym for node value type */
  using node_value_type = V;
  typedef E edge_value_type;
  /** Synonym for size type */
  using size_type = unsigned;
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() { }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODE ITERATOR
  //

  /** @class Graph::NodeIterator
  * @brief Class representing graph's node iterator
  */


  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */


  class Node: private totally_ordered<Node> {
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
     *
     * Representation Invariants:
     * RI(Node): g_ != nullptr
     * RI(Node): uid_ < g_ -> num_nodes()
     */


    Node() {
    }

    /** Return this node's value */
    node_value_type& value() {
      return graph_->nodeData[uid_].second;
    }

    /** Return this node's value - const */
    const node_value_type& value() const {
      return graph_->nodeData[uid_].second;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodeData[uid_].first;
    }

    /** Return this node's position (modifiable) */
    Point& position() {
      return graph_->nodeData[uid_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.uid_==uid_ && n.graph_==graph_);
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
      return (uid_ < n.uid_ && graph_ == n.graph_) || (graph_ < n.graph_);
    }

    size_type degree() const {
      return graph_->adjList[uid_].size();
    }

    incident_iterator edge_begin() const {
      return incident_iterator(graph_, uid_, 0);
    }

    incident_iterator edge_end() const {
      //return incident_iterator(graph_, uid_, this->degree() );
      return incident_iterator(graph_, uid_, graph_->adjList[uid_].size() );
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;
    size_type uid_;

    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        assert(graph_ != nullptr);
        assert(uid_ <= graph_->num_nodes());
      }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodeData.size();
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
  Node add_node(const Point& position, const node_value_type& passedValue = node_value_type()) {
    // add the Point to the vector holding the Graph data
    nodeData.push_back(std::make_pair(position,passedValue));
    std::vector<std::pair<unsigned,edge_value_type>> v;
    adjList.push_back( v );

    return Node(this, num_nodes()-1);
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

  /* *
  * remove_node
  * @brief Removes node @a n from graph object
  *
  * @post Invalidates node iterator for @a n
  * @post num_nodes() = old num_nodes() -1
  * Complexity: ~ O(N)
  *
  */

  size_type remove_node(const Node& n) {
    size_type returnVal = 0;

    if (n.index() < num_nodes()) {
      returnVal = 1;
    }

    /* Remove edges from neighbors in adjList */
    std::vector<Edge> v;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      v.push_back(*it);
    }

    for (size_type i = 0; i < v.size(); i++) {
      remove_edge(v[i]);
    }


    /* Remove node's inner vector from adjList */
    adjList.erase(adjList.begin() + n.index());

    /* Remove 'n' entry from nodeData */
    nodeData.erase(nodeData.begin() + n.index());

    /* Decrement all nodes above 'n' in adjList */
    for (size_type i = 0; i < adjList.size(); ++i) {
      for (size_type j = 0; j < adjList[i].size(); ++j) {
        if (adjList[i][j].first > n.index() ) {
          --adjList[i][j].first;
        }
      }
    }

    return returnVal;
  }

  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    remove_node(n);
    return node_begin();
  }


  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   *
   * Representation invariant:
   * RI(Edge): assert(n1_index_ != n2_index_);
   */
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_,n1_index_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,n2_index_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool dir1 = (e.n1_index_ == n1_index_ && e.n2_index_ == n2_index_);
      bool dir2 = (e.n1_index_ == n2_index_ && e.n2_index_ == n1_index_);
      return (e.graph_== graph_ && (dir1 == true || dir2 == true));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((n1_index_ < e.n1_index_ && graph_ == e.graph_) || graph_ < e.graph_);
    }

    double length() const {
      Point p1 = node1().position();
      Point p2 = node2().position();
      double length_ = norm(p1 - p2);
      return length_;
    }

    /* value() method
     * @brief Returns an edge's value
     * @pre Edge @a e must be valid edge of graph object
     * Complexity: O(d)
     */
    edge_value_type& value() {
      size_type n2_idx = node2().index();
      for (size_type i = 0; i < graph_->adjList[n1_index_].size(); i++) {
        if (graph_->adjList[n1_index_][i].first == n2_idx) {
          return graph_->adjList[n1_index_][i].second;
        }
      }
    }

    const edge_value_type& value() const {
      return (const edge_value_type&) value();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    size_type n1_index_;
    size_type n2_index_;

    Edge(const Graph* graph, size_type n1_index, size_type n2_index  )
      : graph_(const_cast<Graph*>(graph)), n1_index_(n1_index), n2_index_(n2_index) {
        assert(n1_index_ != n2_index_);
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type sum = 0;
    for (size_type i = 0; i < adjList.size(); ++i) {
      sum = sum + adjList[i].size();
    }
    return sum/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

   Edge edge(size_type i) const {
     size_type sum = 0;
     edge_iterator ei = edge_begin();
     while (sum < i) {
       ++ei;
       ++sum;
     }
     return (*ei);
   }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for(auto ni = a.edge_begin(); ni != a.edge_end(); ++ni) {
      if (b.index() == (*ni).node2().index()) {
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
    // check if edge exists
    for (size_type i = 0; i < adjList[a.uid_].size(); i++) {
      if (adjList[a.uid_][i].first == b.uid_) {
        return Edge(this, a.index(),b.index());
      }
   }

    adjList[a.uid_].push_back(std::make_pair(b.uid_,0));
    adjList[b.uid_].push_back(std::make_pair(a.uid_,0));

    return Edge(this, a.index(),b.index());
  }

  /* *
  * remove_edge
  * @brief Removes edge @a e from graph object
  *
  * @post Invalidates edge iterator for @a e
  * @post num_edges() = old num_edges() -1
  * Complexity: O(D) (where D is max degree in graph)
  *
  */

  size_type remove_edge(const Edge& e) {

    Node n1 = e.node1();
    Node n2 = e.node2();
    size_type returnVal = 0;

    /* Remove first instance of edge in adjList */
    for (size_type i = 0; i < adjList[n1.index()].size(); i++) {
      if (adjList[n1.index()][i].first == n2.index()) {
        adjList[n1.index()].erase(adjList[n1.index()].begin() + i);
        returnVal = 1;
      }
    }

    /* Remove second instance of edge in adjList */
    for (size_type i = 0; i < adjList[n2.index()].size(); i++) {
      if (adjList[n2.index()][i].first == n1.index()) {
        adjList[n2.index()].erase(adjList[n2.index()].begin() + i);
        returnVal = 1;
      }
    }
    return returnVal;
  }

  size_type remove_edge(const Node& n1 , const Node& n2) {
    Edge e = Edge(this, n1.index(), n2.index());
    size_type r = remove_edge(e);
    return r;
  }


  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    remove_edge(e);
    return edge_begin();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodeData.clear();
    adjList.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator.
   *
   * Representation invariant:
   * currNode <= graph_->num_nodes()
   */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
      : currNode(0), graph_(NULL) {
    }

    Node operator*() const {
      return Node(graph_, currNode);
    }

    node_iterator& operator++() {
      currNode++;
      return *this;
    }

    bool operator==(const node_iterator& n) const {
      return (currNode == n.currNode && graph_ == n.graph_);
    }


   private:
    friend class Graph;

    size_type currNode;
    Graph* graph_;

    NodeIterator(const Graph* graph, int passedInt) {
      graph_ = const_cast<Graph*>(graph);
      currNode = passedInt;
      assert(currNode <= graph_->num_nodes());
    }
  };

  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  node_iterator node_end() const {
    return NodeIterator(this, this->num_nodes() );
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
    /*IncidentIterator() : main_id_(0), curr_nid_(0), graph_(NULL)
      {
    }*/

    Edge operator*() const {
      size_type node_id = graph_->adjList[main_id_][curr_nid_].first;
      return Edge(graph_, main_id_, node_id);
    }

    incident_iterator& operator++() {
      curr_nid_++;
      return *this;
    }

    bool operator==(const incident_iterator& iit) const {
      return (iit.graph_==graph_ && iit.main_id_==main_id_ && iit.curr_nid_==curr_nid_);
    }

   private:
    friend class Graph;

    Graph* graph_;
    size_type main_id_;
    size_type curr_nid_;

    IncidentIterator(const Graph* graph, size_type main_id, size_type curr_nid) {
      graph_ = const_cast<Graph*>(graph);
      main_id_ = main_id;
      curr_nid_ = curr_nid;
      assert(main_id_ < graph_->adjList.size());
      assert(curr_nid_ <= graph_->adjList[main_id].size());
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. *
   * Representation invariant
   * RI(EdgeIterator): n1_id_ <= graph_->num_nodes()
   * RI(EdgeIterator): n2_id_ <= adjList[n1_id_].size()
   */
  class EdgeIterator : private totally_ordered <EdgeIterator> {
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

    /* Dereferences EdgeIterator @a ei
     * @return Edge
     */
    Edge operator*() const {
      size_type node2_id = graph_->adjList[n1_id_][n2_id_].first;
      return Edge(graph_,n1_id_,node2_id);
    }

    /* Increments EdgeIterator @a ei
     * @return reference to edge
     */
    EdgeIterator& operator++() {
      n2_id_++;

      // if reached end of inner vector
      if (n2_id_ >= graph_->adjList[n1_id_].size()) {
        ++n1_id_;
        n2_id_ = 0;
      }

      while (n1_id_ < graph_->adjList.size() && n1_id_ > graph_->adjList[n1_id_][n2_id_].first) {
        // if reached end of inner vector
        if (n2_id_>=(graph_->adjList[n1_id_].size()-1)) {
          ++n1_id_; // go to next entry of outer vector
          n2_id_ = 0; // reset to start of inner vector

        } else {
        // if still in inner vector
          n2_id_++;
        }
      }

      return *this;
    }

    /* Tests if two edge iterators are equivalent
     * @return boolean: true if equivalent, false otherwise
     */
    bool operator==(const EdgeIterator& e_it) const {
      return (e_it.graph_==graph_ && e_it.n1_id_==n1_id_ && e_it.n2_id_==n2_id_);
    }

   private:
    friend class Graph;

    Graph* graph_;
    size_type n1_id_;
    size_type n2_id_;

    EdgeIterator(const Graph* graph, size_type n1_id, size_type n2_id) {
      graph_ = const_cast<Graph*>(graph);
      n1_id_ = n1_id;
      n2_id_ = n2_id;
    }



  };


  edge_iterator edge_begin() const {
    return EdgeIterator(this,0,0);
  }

  edge_iterator edge_end() const {
    size_type noNodes = adjList.size();
    return EdgeIterator(this, noNodes, 0);
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector< std::pair<Point,node_value_type> > nodeData;
  std::vector< std::vector<std::pair<size_type,edge_value_type> > > adjList;

};

#endif // CME212_GRAPH_HPP
