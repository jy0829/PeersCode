#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>

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

  /** Predeclare the internal struct */
  struct node_internal;
  struct edge_internal;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

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

  /** Type of lists storing internal graph objects: nodes, edges,
   *      adjacency list
   */
  using size_type = unsigned;
  std::vector<node_internal> node_list;
  std::vector<edge_internal> edge_list;
  std::vector<std::vector<std::pair<size_type, size_type>>> adj_list;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
      : node_list(), edge_list(), adj_list(){
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
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->node_list[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Return this node's value.
     *  @return For any internal node @a a, the object specified as the val
     *      attribute of @a a.
     */
    node_value_type& value(){
      return graph_->node_list[uid_].val;
    }

    /** Return this node's value if const.
     *  @return For the internal node @a a associated to this node, the object specified
     *      as the @a val attribute of @a a.
     */
    const node_value_type& value() const{
      return graph_->node_list[uid_].val;
    }

    /** Modify node value using specified integer input. Function takes in a value @a i
     *      and passes that to the @a val attribute of this node. Returns this value.
     *  @param[in] i   Integer input to be passed to node value.
     *  @return For the internal node @a a associated to this node, returns
     *      the object @a i now specified as the @a val attribute of @a a.
     */
    node_value_type& value(int i){
      graph_-> node_list[uid_].val = i;
      return graph_->node_list[uid_].val;
    }

    /** Returns this node's degree, the number of adjacent edges.
     *  @return  Returns the number of edges adjacent to this node.
     *  @pre     0 < size of the adjacency list.
     *  @post    Output is not less than zero and must be strictly less than
     *             total number of nodes in the graph.
     */
    size_type degree() const{
      return graph_->adj_list[uid_].size();
    }

    /** Return initial value for an incident-edge iterator.
     *  @return Returns the incident iterator pointing to the lowest-ordered
     *            edge in the graph.
     *  @pre    Graph must contain at least one edge to call edge with unique id 0.
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, this, 0);
    }

    /** Return end value for an incident-edge iterator.
     *  @return Returns the incident iterator pointing to the edge in the graph
     *            @a i such that, for all edges in the graph @a j, the unique
     *            identifier of edge @a j is less than that of edge @a i.
     *  @pre    Graph must contain at least one edge for a call to the adjacency list.
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, this, degree()-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (graph_ == n.graph_ && uid_ == n.uid_){
        return true;
      }
      return false;
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
      if (Node::uid_ < n.uid_ && Node::graph_ == n.graph_){
        return true;
      }
      if (Node::graph_ < n.graph_){
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid){
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    if (node_list.size() > 0){
      return node_list.size();
    }
    return 0;
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
    node_list.push_back(node_internal(position, val));
    adj_list.push_back(std::vector<std::pair<size_type, size_type>>());
    return Node(this, node_list.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this == n.graph_ && n.uid_ < node_list.size()){
      return true;
    }
    return false;
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Graph::Node(graph_, node_a_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Graph::Node(graph_, node_b_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((graph_ == e.graph_) &&
        ((node_a_id_ == e.node_a_id_ && node_b_id_ == e.node_b_id_) ||
        (node_a_id_ == e.node_b_id_ && node_b_id_ == e.node_a_id_))){
        return true;
      }
    return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_ && node_a_id_ < e.node_a_id_ && node_b_id_ < e.node_b_id_){
        return true;
      }
      if (graph_ < e.graph_){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    Graph* graph_;
    unsigned int node_a_id_;
    unsigned int node_b_id_;
    size_type edge_uid_;
    Edge(const Graph* graph, unsigned int node_a_id, unsigned int node_b_id, size_type edge_uid)
      : graph_(const_cast<Graph*>(graph)), node_a_id_(node_a_id), node_b_id_(node_b_id), edge_uid_(edge_uid){
    }
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
    return Edge(this, edge_list[i].node_ids[0], edge_list[i].node_ids[1], i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned int i = 0; i < adj_list[a.uid_].size(); i++){
      if(adj_list[a.uid_][i] == b.uid_){
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
    //check if already have edge and return if so
    for (unsigned int i = 0; i<adj_list[a.uid_].size(); i++){
      if (b.uid_ == adj_list[a.uid_][i].first){
        return Edge(this, a.uid_, b.uid_, adj_list[a.uid_][i].second);
      }
    }
    //otherwise, push to edge list and adjacency list
    adj_list[a.uid_].push_back(std::make_pair(b.uid_, edge_list.size()-1));
    adj_list[b.uid_].push_back(std::make_pair(a.uid_, edge_list.size()-1));
    if(a<b){
      edge_list.push_back(edge_internal(a, b));
    }
    else{
      edge_list.push_back(edge_internal(b, a));
    }
    //return the edge added
    return Edge(this, a.uid_, b.uid_, edge_list.size()-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_list.clear();
    edge_list.clear();
    adj_list.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {

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

    /** Call the dereferencing operator on NodeIterator.
     *  @return For this iteration, indicated by class variable @a p_, along graph,
     *              returns Node of Node type associated to unique identifier @a p_.
     *  @pre    Class variable @a p_ must be assigned in the set [0, number of nodes).
     *  @post   For output @a n, index of @a n = @a p_.
     */
    Node operator*() const{
      return Graph::Node(graph_, p_);
    }

    /** Call the increment operator on NodeIterator. Iterates through the nodes of the graph.
     *  @return The NodeIterator object that points to the "next" node of the graph,
     *              i.e. the one with unique identifier equal to the current node's unique
     *              identifier, plus 1. Returns valid object but without deferencability
     *              if the node is the node @a i such that for all nodes @a j in the graph,
     *              the index of @a j < the index of @a i.
     */
    NodeIterator& operator++(){
      p_ = p_+1;
      return *this;
    }

    /**  Calls equality operator on NodeIterator. Tests if two NodeIterators are the same.
     *   @return Boolean object that is true if the NodeIterators point to the same node,
     *               and false if they point to different nodes. It also returns true for
     *               two iterators with index outside [0, number of nodes), but if and only
     *               if this index @a p_ is the same.
     */
    bool operator==(const NodeIterator& x) const{
      return p_ == x.p_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    int p_;
    NodeIterator(const Graph* graph, int p)
      : graph_(const_cast<Graph*>(graph)), p_(p){
    }

  };

  /**  Initial value for the node_iterator of node_iterator type. The
   *       node iterator starts at this value in the graph when iterating over
   *       nodes in the graph.
   *   @return NodeIterator that dereferences to the node of index 0.
   */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /**  End value for node_iterator, of node_iterator type. The node iterator
   *       halts at this value when iterating over all nodes in the graph.
   *   @return NodeIterator that dereferences to the node i whose index
   *       is such that for any node j, the index of j is less than the
   *       index of i.
   *   @pre    Graph must contain at least one node.
   */
  node_iterator node_end() const{
    return NodeIterator(this, node_list.size()-1);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered <IncidentIterator> {
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

    /**  Dereference operator for the incident iterator. Returns the Edge
     *       that the iterator defines in the vector of nodes adjacent to a
     *       node specified by the class.
     *   @return An edge i such that, for class attributes home node @a n_ and
     *       iterating position @a q_, i is the edge whose first node is @a
     *       n_ and whose second node has index q_.
     *   @pre    Graph contains at least @a q_ edges incident to node @a n_.
     *   @post   Returned edge is a valid edge of the graph.
     */
    Edge operator*() const{
      Node home_node = *n_;
      return Edge(graph_, home_node.uid_, graph_->adj_list[home_node.uid_][q_].first,
          graph_->adj_list[home_node.uid_][q_].second);
    }

    /**  Increment operator for incident iterator. Iterates over the edges adjacent
     *       to a specified node.
     *   @return An IncidentIterator object that references an edge next to the edge
     *       currently referenced by the iterator, where "next" here means
     *       "subsequent in a vector of adjacent node indices." If the edge currently
     *       referenced is the edge i whose adjacency index is not less than the adjacency
     *       index of any other edge j in the set of nodes adjacent to the node referenced
     *       by class variable @a n_, then incrementing i will not give an iterator
     *       with a valid dereference.
     */
    IncidentIterator& operator++(){
      q_ = q_+1;
      return *this;
    }

    /**  Equality operator for the incident iterator. Tests if two iterators
     *       point to the same edge.
     *   @param[in] x  Incident iterator to compare this iterator to.
     *   @return       Returns true if @a x and the current iterator have the
     *                     same base node and the same adjacency index, indicating
     *                     that the iterators dereference to the same edge.
     *                     Otherwise returns false.
     *   @pre          Base node @a n_ and adjacency index @q q_ must be initialized
     *                     as attributes of this iterator and the iterator @a x.
     */
    bool operator==(const IncidentIterator& x) const{
      return (n_ == x.n_ && q_ ==x.q_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    Node* n_;
    int q_;
    IncidentIterator(const Graph* graph, const Node* n, int q)
      : graph_(const_cast<Graph*>(graph)), n_(const_cast<Node*>(n)), q_(q){
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /** Dereference an EdgeIterator object.
     *  @return For this iteration, indicated by class variable @a r_, along graph,
     *              returns edge of Edge type with unique index @a r_.
     *  @pre    Class variable @a r_ must be assigned in the set [0, number of edges).
     *  @post   For output @a e, edge index of @a e = @a r_.
     */
    Edge operator*() const{
      return graph_->edge(r_);
    }

    /** Call the increment operator on EdgeIterator. Iterates through the edges of the graph.
     *  @return The EdgeIterator object that points to the "next" edge of the graph,
     *              i.e. the one with index equal to the current edge's index plus 1.
     *              Returns valid object but without deferencability if the class variable
     *              @a r_ plus 1 is not in the set [0, number of edges), such as if the
     *              edge is the edge in the graph such that all other edges in the graph
     *              have lesser index.
     */
    EdgeIterator& operator++(){
      r_ = r_+1;
      return *this;
    }

    /**  Equality operator on EdgeIterator. Tests if two EdgeIterators are the same.
     *   @return Boolean object that is true if the EdgeIterators point to the same edge,
     *               and false if they point to different edges. It also returns true for
     *               two iterators with iterative index outside [0, number of edges),
     *               but if and only if this index (denoted by class variable @a r_ is
     *               the same.
     */
    bool operator==(const EdgeIterator& x) const{
      if (graph_ == x.graph_ && r_ == x.r_){
        return true;
      }
    return false;
    }

   private:
    friend class Graph;
    Graph* graph_;
    unsigned int r_;
    EdgeIterator(const Graph* graph, unsigned int r)
      : graph_(const_cast<Graph*>(graph)), r_(r){
    }
  };

  /**  Initial value for the edge iterator of edge_iterator type. The
   *       edge iterator starts at this reference to an edge when
   *       iterating over edge in the graph.
   *   @return EdgeIterator that dereferences to the edge of index 0.
   */
  edge_iterator edge_begin() const{
    return edge_iterator(this, 0);
  }

  /**  End value for edge_iterator, of edge_iterator type. The edge iterator
   *       halts at this value when iterating over all edges in the graph.
   *   @return EdgeIterator that dereferences to the edge i whose index
   *       is such that for any node j, the index of j is less than the
   *       index of i.
   */
  edge_iterator edge_end() const{
    return edge_iterator(this, edge_list.size());
  }

  private:
    //declare internal struct for nodes
    struct node_internal{
      node_value_type& value();
      const node_value_type& value() const;


      node_internal(Point p){
        position = p;
        node_value_type val;
      }
      node_internal(Point p, node_value_type v){
        position = p;
        val = v;
      }
      Point position;
      node_value_type val;
    };

    //declare internal struct for edges
    struct edge_internal{

      edge_internal(Node node_a, Node node_b){
        node_ids.push_back(node_a.uid_);
        node_ids.push_back(node_b.uid_);
      }

      std::vector<int> node_ids;
      size_type edge_uid;
    };

    node_internal* nodes_;
    edge_internal* edges_;
    size_type size_;

};

#endif // CME212_GRAPH_HPP
