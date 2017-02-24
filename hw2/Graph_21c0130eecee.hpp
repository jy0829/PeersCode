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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  using point_type = Point;

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
  
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  // better to use member initializer list
  // same result as assignment but faster
  Graph()
  :points_(), adjacency_(), nedges_(0) {
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
  class Node: totally_ordered<Node> {
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
      // just do nothing
    }
    

    /** Return this node's position. */
    const Point& position() const {
        return graph_->points_[uid_].first;

    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    
    /**
     * Return this node's value
     * @pre this is a valid node:
     *      graph_ is a valid pointer to graph
     *      uid_ is a key of graph_->points
     * @post graph_->points_[uid_].second = result
     */
    node_value_type& value(){
      return graph_->points_[uid_].second;
    }
    
    /**
     * Const version of value()
     */
    const node_value_type& value() const{
      return graph_->points_[uid_].second;
    }
    
    /**
     * Return this node's degree (number of adjacent nodes)
     * @pre this is a valid node
     *      graph_ is a valid pointer to graph
     *      uid_ is a key of graph_->adjacency_
     * @post graph_->adjacency_[uid_].size() = result
     */
    size_type degree() const{
      return graph_->adjacency_[uid_].size();
    }
    
    /**
    * Return incident_iterator begin
    *
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, uid_, 0);
    }

    /**
    * Return incident_iterator end
    */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_, degree());
    }

    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((n.graph_ == graph_) and (n.uid_ == uid_));
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * 
     * Ordering is the lexicographic order on (graph_, uid_)
     */
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE
      return (graph_ < n.graph_) or (graph_ == n.graph_ and uid_ < n.uid_);
    }
    
    
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer back to the Graph container
    graph_type* graph_;

    // This element's unique id
    size_type uid_;
    
    

    // Private constructor
    Node(const graph_type* graph, size_type uid)
        : graph_(const_cast<graph_type*>(graph)), uid_(uid){
    }
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return points_.size();
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
    // HW0: YOUR CODE HERE
    points_.push_back({position, value});
    adjacency_.push_back(std::vector<size_type>());
    return Node(this, size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ((n.graph_ == this) and (n.uid_ < size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert (i < num_nodes());
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
  class Edge: totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      // just do nothing
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, id1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, id2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((e.id1_ == id1_) and (e.id2_ == id2_)
              and (e.graph_ == graph_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * 
     * Ordering is lexicographic order on (graph_, id1_, id2_)
     */
    bool operator<(const Edge& e) const {
        // HW0 TODO global order
        // TODO can you compare pointers?
      return ((graph_ < e.graph_) or 
      (graph_ == e.graph_ and ((e.id1_ < id1_) or (e.id1_ == id1_ and e.id2_ < id2_))));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type id1_;
    size_type id2_;
    
    // Private constructor
    Edge(const graph_type* graph, size_type id1, size_type id2)
        : graph_(const_cast<graph_type*>(graph)), id1_(id1), id2_(id2) {
    }


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * We keep track of the number of edges while we add them
   * Thus, we don't need to iterate over the adjacency vector
   * 
   * @pre nedges_ = sum(e.size() for e in adjacency_) / 2
   * 
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE 
    // 1. non-efficient method
    // unsigned nedges = 0;
    // for (auto e: adjacency_){
    //   nedges += e.size();
    // }
    // assert(nedges_ == nedges/2);
    // return nedges/2;

    // 2. efficient method
    return nedges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
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
    // HW0: YOUR CODE HERE
    assert(this == a.graph_ and this == b.graph_);
    assert(a.uid_ < size() and b.uid_ < size());
    auto v = adjacency_[a.uid_];
    return std::find(v.begin(), v.end(), b.uid_) != v.end();
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
    // HW0: YOUR CODE HER
    if (!(has_edge(a, b))){
      adjacency_[a.uid_].push_back(b.uid_);
      adjacency_[b.uid_].push_back(a.uid_);
      ++nedges_;
      
    }
    assert(has_edge(a, b) and has_edge(b, a));
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
      points_.clear();
      adjacency_.clear();
      nedges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: totally_ordered<NodeIterator> {
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
    
    /**
    * Dereference op of NodeIterator, returns a Node
    * @pre n_.uid_ is an id in [0, n.graph_->size)
    * @pre n_ is a valid node
    */
    Node operator*() const{
      return n_;
    }
    
    /**
    * Increment op of NodeIterator, returns a reference to a NodeIterator
    * @post n_.uid_ = n_.uid_ + 1
    */
    NodeIterator& operator++(){
      ++n_.uid_;
      return *this;
    }
    
    /**
    * Equality op of NodeIterator, return bool
    * @param[in] it, other NodeIterator to compare with
    * @pre this and @a it are valid iterators
    * 
    */ 
    bool operator==(const NodeIterator& it) const{
      return (n_.graph_ == it.n_.graph_) and (n_.uid_ == it.n_.uid_);
    }
    

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    
    /**A node storing the current node
    * n_.graph_ is a valid pointer to the graph
    * n_.uid_ can actually be = n.graph_.size()
    */
    Node n_;
    
    /**
    * Constructor for NodeIter
    * @param n, Node object
    */
    NodeIterator(Node& n)
    : n_(n){
    }
    
    
  };

  // HW1 #2: YOUR CODE HERE
  /**
  * Return iterator begin
  * @post result.n is a Node s.t.
  *     - graph_ = this
  *     - uid_ = 0 (if there is no node in the graph, this is NOT a valid node)
  */
  node_iterator node_begin() const{
    Node n = Node(this, 0);
    return NodeIterator(n);
  }

  /**
  * Return iterator end
  * @post result.n_ is a Node s.t.
  *   - graph_ = this
  *   - uid_ = size() (thus it is NOT a valid node)
  */
  node_iterator node_end() const{
    Node n = Node(this, size());
    return NodeIterator(n);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: totally_ordered<IncidentIterator> {
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

    /**
    * Dereference operator for IncidentIterator
    * @pre graph_ is a valid graph
    * @pre centerid_ is a valid id for a node of the graph
    * @pre outerpos_ is a valid index in adjacency_[centerid_]
    *    - outerpos_ in [0, ajdacency_[centerid_].size())
    * @post result is a valid edge of the graph
    */
    Edge operator*() const{
      size_type outerid = (graph_->adjacency_[centerid_][outerpos_]);
      return Edge(graph_, centerid_, outerid);
    }
    

    /**
    * Increment operator 
    * @post outerpost is incremented
    */
    IncidentIterator& operator++(){
      ++outerpos_;
      return *this;
    }
    
    /**
    * Equality operator
    * @param[in] @a it IncidentIterator to compare with this
    * @pre @a it and this are valid iterators
    */
    bool operator==(const IncidentIterator& it) const{
      return ((graph_    == it.graph_) and
              (centerid_ == it.centerid_) and
              (outerpos_ == it.outerpos_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    
    // pointer to the graph
    graph_type* graph_;
    // id of the center node 
    size_type centerid_;
    // position of the outer node in the adjacency_[centerid_]
    size_type outerpos_;
    
    /**
    * Constructor
    * @param @a graph a valid pointer to graph
    * @param @a centerid a valid id for a node
    * @param @a outerpos valid index in the adjacency list of @a centerid_
    * @pre @a centerid_ in [0, @a graph_->ajacency_.size())
    * @pre @a outerpos in [0, @a graph_->ajdacency_[@a centerid_].size())
    */
    IncidentIterator(const graph_type* graph, size_type centerid, size_type outerpos)
    : graph_(const_cast<graph_type*>(graph)), centerid_(centerid), outerpos_(outerpos){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: totally_ordered<EdgeIterator> {
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
    /**
    * Dereference operator
    * @pre graph_ is a valid graph
    * @pre centeri_ is a valid nodeitetor on graph_
    * @pre centeri_ != graph_node_end()
    * @pre outi_ is a valid incident iterator of the *centeri
    * @pre outi_ != (*centeri_).edge_end()
    * @post result is a valid edge of the graph_
    */
    Edge operator*() const{
      return Edge(graph_, (*centeri_).index(), (*outi_).node2().index());
    }

    /**
    * Increment operator
    * @pre graph is a valid graph
    * @pre centeri_ is a valid nodeitetor on graph_
    * @pre outi_ is a valid incident iterator of the *centeri
    * @post result.graph_ = graph_
    * @post if there is i, j st 
    *               j is in graph_->adjacency[i]
    *               and j > i
    *               and i >= *centeri_.index()
    *               and if (i == *centeri_.index()), then 
    *                  j is the first j after outi_.node2.index in ajdacency[i] st. j > i
    *               and if i > *centeri_.index, then
    *                   j is the first index in ajdacency[i] st. j > i
    *       then result points to this edge
    *       else
    *           result is centeri_ = node_end() and outi_ = Incidentiterator(centeri, 0)
    *         
    * 
    */
    EdgeIterator& operator++() {
      ++outi_;
      while(centeri_ != graph_->node_end()){
        while(outi_ != (*centeri_).edge_end()){
          if ((*outi_).node2().index() > (*centeri_).index()) return *this;
          ++outi_;
        }
        ++centeri_;
        outi_ = (*centeri_).edge_begin();
      }

      return *this;
    }

    /**
    * Equality operator
    * @param ei other edge it to compare with
    */
    bool operator==(const EdgeIterator& ei) const {
      return (graph_ == ei.graph_) and (centeri_ == ei.centeri_) and (outi_ == ei.outi_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // pointer to the graph, not necessary (because already in centeri and outi)
    // but more practical
    graph_type* graph_;
    // iterator over nodes
    NodeIterator centeri_;
    // iterator over the incident edges of *centeri
    IncidentIterator outi_;

    /**
    * Constructor
    * @pre graph_ is a valid graph
    * @pre centeri_ is a valid nodeitetor on graph_
    * @pre outi_ is a valid incident iterator of the *centeri
    */
    EdgeIterator(const graph_type* graph, NodeIterator& centeri, IncidentIterator& outi)
    : graph_(const_cast<graph_type*>(graph)), centeri_(centeri), outi_(outi){
    }

  };

  // HW1 #5: YOUR CODE HERE
  /**
  * Return edge iterator begin
  * @pre centeri = node_begin
  * @pre outi = edge_begin of centeri
  */
  edge_iterator edge_begin() const {
    NodeIterator centeri = node_begin();
    IncidentIterator outi = (*centeri).edge_begin();
    return EdgeIterator(this, centeri, outi);
  }

  /**
  * Return edge iterator end
  * @pre centeri = node_end
  * @pre outi = Incidentoperator(size(), 0) (for cohrence with ++)
  */
  edge_iterator edge_end() const {
    NodeIterator centeri = node_end();
    IncidentIterator outi = IncidentIterator(this, size(), 0);
    return EdgeIterator(this, centeri, outi);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<std::pair<point_type, node_value_type>> points_;
  std::vector<std::vector<size_type>> adjacency_;
  size_type nedges_;  
};



#endif // CME212_GRAPH_HPP
