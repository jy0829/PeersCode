#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 *  @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 *  @brief A template for 3D undirected graphs.
 *
 *  Users can add and retrieve nodes and edges. Edges are unique (there is at
 *  most one edge between any pair of distinct nodes).
 */

template <typename V, typename E>

class Graph {

 private:

  //
  //  PRIVATE TYPE DEFINITIONS
  //
  
  // Internal structure for node.
  struct node_info;
  
  // Internal structure for edge.
  struct edge_info;

 public:

  //
  //  PUBLIC TYPE DEFINITIONS
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

  /** The V argument of the Graph template. */
  using node_value_type = V; 

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** The E argument of the Graph template. */
  using edge_value_type = E;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of sizes. */
  using size_type = unsigned;
  
  /** Type of unique id's of nodes. */    
  using uid_type = unsigned;

  /** Type of indices in vector i2u. */
  using idx_type = int;

  //
  //  CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), i2u_(){
  }

  /** Default destructor */
  ~Graph() = default;

  //
  //  NODES
  //

  /** @class Graph::Node
   *  @brief Class representing the graph's nodes.
   *
   *  Node objects are used to access information about the Graph's nodes.
   */
  class Node :private totally_ordered<Node> {
   
   public:
    /** Construct an invalid node. */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_-> nodes_[uid_].position_;
    }

    /** Return this node's unique id. */
    uid_type uid() const {
      return uid_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    idx_type index() const {
      return graph_ -> nodes_[uid_].idx_;
    }

    /** Return this node's position. */
    Point& position(){
      return graph_-> nodes_[uid_].position_;
    }
    
    /** Return this node's value. */
    node_value_type& value(){
      return graph_-> nodes_[uid_].value_;
    }

    /** Return this node's value (overriding for correctness).*/
    const node_value_type& value() const{
      return graph_-> nodes_[uid_].value_;
    }
 
    /** Return this node's degree, a number in the range [0, @num_edges()).
     *
     *  Recall that the degree() of a node is the number of edges incident
     *  to it.
     *
     *  @result: Number in [0, @a num_edges()) 
     */
    size_type degree() const{
      return graph_ -> nodes_[uid_].adj_.size();
    }

    /** Returns first position of data structure containing the adjacency 
     *  information of node @a uid_
     *
     *  @result: Position of the first element of the vector 
     *           nodes_[uid_]_adj_.
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, uid_, 0);
    }

    /** Returns last position of data structure containing the adjacency      
     *  information of node @a uid_.
     *
     *  @result: One position past the abstract range of the vector 
     *           @a nodes_[@a uid_]_adj_.
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *  Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (uid_ == n.uid() && graph_ == n.graph_);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     *  This ordering function is useful for STL containers such as
     *  std::map<>. It need not have any geometric meaning.
     *
     *  The node ordering relation must obey trichotomy: For any two nodes x
     *  and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if(graph_ == n.graph_){
        return (uid_ < n.uid());
      }
      return graph_ < n.graph_;
    }

   private:
    
    Graph* graph_; // Pointer back to the Graph container
    size_type uid_; // This node's unique id

    /** Private Constructor 
     * 
     *  The job of the private constructor of node is to establish and 
     *  defend the representation invariance of Graph::Node. If the Graph 
     *  ever tries to construct an invalid node, it will be detected by this 
     *  constructor.
     *
     *  RI(Node) = graph_ -> i2u_[g -> nodes_[uid_].idx_] == uid_.
     *
     */
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid){
      assert(graph_ -> i2u_[graph_ -> nodes_[uid_].idx_] == uid_);
    }

    // Allow Graph to access Node's private member data and functions.*/
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   *  Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   *  @param[in] position The new node's position
   *  @param[in] value The new node's value
   *  @post new num_nodes() == old num_nodes() + 1
   *  @post result_node.index() == old num_nodes()
   *
   *  Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    // Add uid_ of new node to vector of valid nodes.
    i2u_.push_back(nodes_.size());

    // Create and populate struct node_info of node being added.
    node_info n;
    n.idx_ = i2u_.size() -1;
    n.position_ = position;
    n.value_ = value;

    // Add adjacency information vector for new node
    std::vector<edge_info> adj;
    n.adj_ = adj;

    // Add new node to nodes_ vector
    nodes_.push_back(n); 

    // Return new node
    return Node(this,nodes_.size()-1);   
  }

  /** Determine if a Node belongs to this Graph
   *  @return True if @a n is currently a Node of this Graph
   *
   *  Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_ == this && n.index() < (int)size() && n.index() != -1);
  }

  /** Return the node with index == i.
   *  @pre 0 <= i < size()
   *  @post result_node.idx_ == i
   *
   *  Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i2u_[i]);        
  }

   /** Return the node with uid == i.
   *  @post result_node.uid_ == i
   *
   *  Complexity: O(1).
   */
  Node node_uid(size_type i) const {
    return Node(this,i);        
  }

  /** Remove the Node n and all Edges that are incident to it.
   *  @param[in] n  Node to be removed
   *
   *  @post If node was present on the graph:
   * 
   *        1) new size() = old size() - degree
   *
   *        2) All edges incident to Node @a n are removed; that is,
   *           if Node b was incident to Node @a n: 
   *           new b.degree() = old b.degree() - 1, and
   *           new num_edges() = old num_edges() - old n.degree().
   *
   *        3) Node n is invalidated:
   *           new n.index() = -1.
   *
   *        4) All outstanding nodes and edges remain valid; that is,
   *           for all previously valid Nodes m different to Node n,
   *           i2u_[nodes_[m.uid_].idx_] == m.uid_ after the function call.
   *
   *        5) Invalidates NodeIterator node_end()-1 and all NodeIterators 
   *           that previously pointed to @a a.
   *           Invalidates all EdgeIterators and IncidentIterators pointing 
   *           to an Edge in which @a n is incident.
   *           For all Nodes b adjacent to @a, invalidates the Incident Iterator
   *           b.edge_end() - 1.
   *           All other outstanding iterators remain valid.
   *       
   * @result = 1 if node was present in the graph and 0 otherwise.
   *
   *  Complexity: No more than O(size())
   */

  size_type remove_node(const Node& n) {
    if(has_node(n)){
      // Remove all edges incident to n
      auto ii = n.edge_begin();
      while (ii != n.edge_end()) {
        ii = remove_edge(ii);
      }
      // Change the index of the vector that being swapped
      nodes_[i2u_[size()-1]].idx_ = nodes_[n.uid_].idx_; 
      // Delete node from the vector of valid nodes
      i2u_[nodes_[n.uid_].idx_] = i2u_[size()-1];
      i2u_.pop_back();
      // Invalidate the deleted node
      nodes_[n.uid_].idx_ = -1;
    }
    return 0;
  }

  /** Remove the Node n that @a n_it points to and all Edges that are incident 
    * to it.
    * @param[in] @a n_it.
    * @result = n_it,
    * 
    * @post See size_type remove_node(Node& n);
    *   
    * Complexity: No more than O(size())
    *
    */
  node_iterator remove_node(node_iterator n_it){
    remove_node(*n_it);
    return n_it;
  }

  //
  //  EDGES
  //

  /** @class Graph::Edge
   *  @brief Class representing the graph's edges.
   *
   *  Edges are order-insensitive pairs of nodes. Two Edges with the same
   *  nodes are considered equal if they connect the same nodes, in either
   *  order.
   */
  class Edge :private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return the unique id of a node of this Edge */
    Node node1() const {
      return graph_-> node_uid(uid_a_);
    }

    /** Return the unique id of the other node of this Edge */
    Node node2() const {
      return graph_-> node_uid(uid_b_);
    }

    /** Return this edge's value. 
     * @pre graph_-> nodes_[uid_a_].adj_[uid_b_].value_ ==
     *      graph_-> nodes_[uid_b_].adj_[uid_a_].value_;  
     *
     * @pre uid_b_ is an element of nodes_[uid_a_].adj_.
     *
     */
    edge_value_type& value(){
      size_type n = graph_ -> nodes_[uid_a_].adj_.size();
      for(size_type i = 0; i < n; ++i){
        if(graph_ -> nodes_[uid_a_].adj_[i].uid_other_ == uid_b_){
          return graph_-> nodes_[uid_a_].adj_[i].value_ ;
        }
      }
      // To silence compiler warning (this line will never be reached if
      // the preconditions are met).
      return graph_-> nodes_[uid_a_].adj_[0].value_ ;
    }

    /** Return this edge's value (overriding for correctness).*/
    const edge_value_type& value() const{
      size_type n = graph_ -> nodes_[uid_a_].adj_.size();
      for(unsigned int i = 0; i < n; ++i){
        if(graph_ -> nodes_[uid_a_].adj_[i].uid_other_ == uid_b_){
          return graph_-> nodes_[uid_a_].adj_[i].value_ ;
        }
      }
    }
    
    /** Return the length of this Edge; that is, the Euclidian distance 
     *  between the positions of it's nodes. 
     */
    double length() const {
      return norm(graph_ -> nodes_[uid_a_].position_ - 
                  graph_ -> nodes_[uid_b_].position_);
    }

    /** Test whether this edge and @a e are equal.
     *
     *  Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_ &&
	      std::min(uid_a_,uid_b_) == std::min(e.uid_a_, e.uid_b_) &&
              std::max(uid_a_,uid_b_) == std::max(e.uid_a_, e.uid_b_)   );
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     *  This ordering function is useful for STL containers such as
     *  std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(graph_ == e.graph_){
        return ((std::min(uid_a_,uid_b_)  < std::min(e.uid_a_,e.uid_b_)  ||
                (std::min(uid_a_,uid_b_) == std::min(e.uid_a_,e.uid_b_) &&
                 std::max(uid_a_,uid_b_)  < std::max(e.uid_a_,e.uid_b_))));
      }
      return graph_ < e.graph_;
    }

   private: 
    Graph* graph_;    //  Pointer back to the Graph container
    size_type uid_a_; //  This edge's node1() unique id
    size_type uid_b_; //  This edge's node2() unique id

    /** Private Constructor 
     * 
     *  The job of the private constructor of edge is to establish and 
     *  defend the representation invariance of Graph::Edge. If the Graph 
     *  ever tries to construct an invalid edge, it will be detected by this 
     *  constructor.
     *
     * RI(Edge) = uid_a_ != uid_b
     *         && graph_ -> i2u_[graph_ -> nodes_[uid_a_].idx_] == uid_a_
     *         && graph_ -> i2u_[graph_ -> nodes_[uid_b_].idx_] == uid_b_.
     *
     */

    Edge(const Graph* graph, uid_type uid_a, uid_type uid_b)
      : graph_(const_cast<Graph*>(graph)), uid_a_(uid_a), uid_b_(uid_b){
      assert(uid_a_ != uid_b_);
      assert(graph_ -> i2u_[graph_ -> nodes_[uid_a_].idx_] == uid_a_);
      assert(graph_ -> i2u_[graph_ -> nodes_[uid_b_].idx_] == uid_b_);
    }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
  };

  /** Return the total number of edges in the graph.
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return count_edges;
  }

  /** Return the edge with index @a i.
   *  @pre 0 <= @a i < num_edges()
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { 
    return *(std::next(edge_begin(),i));
  }

  /** Test whether two nodes are connected by an edge.
   *  @pre @a a and @a b are valid nodes of this graph
   *  @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check if Node b is adjacent to node a
    for(auto ii = a.edge_begin(); ii != a.edge_end(); ++ii){
      edge_type e = *ii;
      if(e.node2() == b){
        return true;
      }
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already
   *  exists.
   *
   *  @pre @a a and @a b are distinct valid nodes of this graph
   *  @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   *  @post has_edge(@a a, @a b) == true
   *  @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *        Else,                        new num_edges() == old num_edges()+1.
   *
   *  Can invalidate edge indexes -- in other words, old edge(@a i) might not
   *  equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, 
                const Node& b, 
                const edge_value_type& value = edge_value_type()) {
    if(not has_edge(a,b)){
      // Add adjacency information to nodes incident to new edge.
      edge_info e_a, e_b;
      e_a.uid_other_ = b.uid();
      e_b.uid_other_ = a.uid();
      e_a.value_ = value;
      e_b.value_ = value;
      // Push information to both nodes' adjacency vectors
      nodes_[a.uid()].adj_.push_back(e_a);
      nodes_[b.uid()].adj_.push_back(e_b);
      // Increase the edge count by one.
      count_edges += 1;
    }
    return Edge(this,a.uid(),b.uid());
  }

  /* Remove Edge e that connects Node @a a and Node @a b from this graph
   * @param[in] a  First node incident on the edge to be removed,
   *            b  Second node incident on the edge to be removed.
   *
   * @pre @a a and @a b are distinct valid nodes of the graph; that is,
   *      a.uid_ != b.uid_
   *      0 <= graph_ -> i2u_[graph_ -> nodes_[a.uid()].idx_] == a.uid
   *      0 <= graph_ -> i2u_[graph_ -> nodes_[b.uid()].idx_] == b.uid.
   *
   * @post new num_edges() = old num_edges() - 1
   *
   *       @a e is eliminated from nodes_[a.uid()].adj_ & nodes_[b.uid()].adj_
   *       new a.degree() = old a.degree() - 1
   *       new b.degree() = old b.degree() - 1 
   *
   *       Invalidates Edge(this, @a a, @a b), all the oustanding edges remain
   *       valid.
   *
   *       Invalidates all EdgeIterators and IncidentIterators pointing to 
   *       Edge(this, @a a, @a b) or to an edge incident to @a a or @a b. All 
   *       other outstanding iterators remain valid.
   * 
   * 
   * @result = 1 if the edge was succesfully removed and 0 otherwise.
   *
   * Complexity: No more than O(size() + num_edges())
   */

  size_type remove_edge(const Node& a, const Node& b){
    // If the graph has edge e
    if(has_edge(a,b)){
      // Remove edge from nodes_[a.uid()].adj_
      for(size_type i = 0; i < a.degree(); ++i){
        if(nodes_[a.uid_].adj_[i].uid_other_ == b.uid()){
           nodes_[a.uid_].adj_[i] = nodes_[a.uid_].adj_[a.degree()-1];
           nodes_[a.uid_].adj_.pop_back();
        }
      } 
      // Remove edge from nodes_[b.uid()].adj_
      for(size_type i = 0; i < b.degree(); ++i){
        if(nodes_[b.uid_].adj_[i].uid_other_ == a.uid()){
           nodes_[b.uid_].adj_[i] = nodes_[b.uid_].adj_[b.degree()-1];
           nodes_[b.uid_].adj_.pop_back();
        }
      }
      count_edges-=1;
      return 1;
    }
    return 0;
  }

   /* Remove Edge @a e from this graph
   * @param[in] e Edge that will be removed.
   *
   * @pre @a e.node1() and @a e.node2() are distinct valid nodes of the graph; 
   *    e.node1().uid_ != e.node2().uid_
   *    graph_ -> i2u_[graph_ -> nodes_[e.node1().uid()].idx_] == e.node1().uid_
   *    graph_ -> i2u_[graph_ -> nodes_[e.node2().uid()].idx_] == e.node2().uid.
   *
   * @post new num_edges() = old num_edges() - 1
   *
   *       @a e is eliminated from 
   *       nodes_[e.node2().uid()].adj_ & nodes_[e.node2().uid()].adj_
   *       new e.node1().degree() = old e.node1().degree() - 1
   *       new e.node2().degree() = old e.node2().degree() - 1. 
   *
   *       Invalidates the Edge @a e, all the oustanding edges remain valid.
   *
   *       Invalidates all EdgeIterators and IncidentIterators pointing to @a e
   *       or to an edge incident to @a e.node1() or @a e.node2(). All other 
   *       outstanding iterators remain valid.
   * 
   * @result = 1 if the edge was succesfully removed and 0 otherwise.
   *
   * Complexity: No more than O(size() + num_edges())
   */
  
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }
  
   /* Remove Edge e that @a e_it points at from this graph
   * @param[in] e_it Pointer to Edge that will be removed.
   *
   * @pre Edge (*e_it) has to meet all the conditions specified in 
   *      size_type remove_edge(const Edge& e).
   * 
   * @post Same post-conditions than in size_type remove_edge(const Edge& e).
   *
   * @result = e_it;
   *
   * Complexity: No more than O(size() + num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return e_it;
  }
     /* Remove Edge e that @a i_it points at from this graph
   * @param[in] i_it Pointer to Edge that will be removed.
   *
   * @pre Edge (*i_it) has to meet all the conditions specified in 
          size_type remove_edge(const Edge& e).
   * 
   * @post Same post-conditions than in size_type remove_edge(const Edge& e).
   *
   * @result = i_it;
   *
   * Complexity: No more than O(size() + num_edges())
   */
  incident_iterator remove_edge(incident_iterator i_it){
    remove_edge(*i_it);
    return i_it;
  }

  /** Remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *
   *  Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    i2u_.clear();
    count_edges = 0;
  }

  //
  //  Node Iterator
  //

  /** @class Graph::NodeIterator
   *  @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private equality_comparable<NodeIterator>{
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

    /** Dereference operator:
     *  Access node at the current position of the i2u_ vector.
     */ 
    Node operator*() const{
      return Node(graph_,graph_ -> i2u_[curr_]);
    }

    /** Increment iterator:
     *  Increment the current position by one.
     */
    NodeIterator& operator++(){
      ++curr_;
      return *this;
    }
    
    /** Equality operator:
     *  Check if two iterators are equal; two iterators are equal if they have
     *  the same current position in the same graph.
     */
    bool operator==(const NodeIterator& ni) const{
      return (curr_ == ni.curr_ && graph_ == ni.graph_);
    }

   private:
    friend class Graph;
    Graph* graph_; // Pointer to a graph.
    size_type curr_; // Current position in @a i2u_ vector. 

    /** Private Constructor 
     * 
     *  The job of the private constructor of NodeIterator is to establish and 
     *  defend the representation invariance of Graph::NodeIterator.
     *  If the Graph ever tries to construct an invalid node iterator, it will 
     *  be detected by this constructor.
     *
     */
    NodeIterator(const Graph* graph, size_type curr)
      : graph_(const_cast<Graph*>(graph)), curr_(curr){
      assert(graph_ != nullptr);
      assert(curr_ <= graph_ -> size());
    }
  };

  // Returns first position of @a i2u_ vector
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  // Returns last position of  @a i2u_ vector
  node_iterator node_end() const{
    return NodeIterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   *  @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private equality_comparable<IncidentIterator> {
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

    /** Dereference operator:
     *  Access edge that connects node(@a uid_) with the node at the current
     *  position of @a nodes_[@a node_id_].adj_ vector.
     *  This vector contains the adjacency information of node(@a uid_).
     */ 
    Edge operator*() const{
      return Edge(graph_, 
                  node_id_, 
                  graph_ -> nodes_[node_id_].adj_[curr_].uid_other_);
    }
    /** Increment operator:
     *  Increment the current position at @a nodes_[@a node_id_].adj_ by one.
     */
    IncidentIterator& operator++(){
      ++curr_;
      return *this;
    }
    /** Equality operator:
     *  Check if two incident iterators are equal. Two incident iterators are 
     *  equal if they iterate through the same node in the same graph, and if
     *  they share the same current position.
     */
    bool operator==(const IncidentIterator& ii) const{
      return (graph_   == ii.graph_   &&
              node_id_ == ii.node_id_ && 
              curr_    == ii.curr_);
    }

   private:
    friend class Graph;
    Graph* graph_; // Pointer to graph
    uid_type node_id_; // Node you're iterating through
    size_type curr_; // Current position in @a nodes_[@a node_id_].adj_

    /** Private Constructor */
    IncidentIterator(const Graph* graph, uid_type node_id, size_type curr)
        : graph_(const_cast<Graph*>(graph)), node_id_(node_id), curr_(curr){
      assert(graph-> i2u_[graph_ -> nodes_[node_id_].idx_] == node_id_);
      assert(curr_ <= graph_ -> nodes_[node_id_].adj_.size());
    }
  };

  //
  //  Edge Iterator
  //

  /** @class Graph::EdgeIterator
   *  @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private equality_comparable<EdgeIterator>{
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
    /** Dereference operator:
     *  Access edge that connects the nodes stored at the current position of  
     *  the @a edges_ vector.
     */ 
    Edge operator*() const{
      return Edge(graph_, 
                  graph_ -> i2u_[curr_n_], 
                  graph_ -> 
                  nodes_[graph_ -> i2u_[curr_n_]].adj_[curr_e_].uid_other_);
    }
    /** Increment operator:
     *  Increment the current position at @a edges_ by one.
     */
    EdgeIterator& operator++(){
      ++curr_e_;
      fix();
      return *this;
    }
    /** Equality operator:
     *  Check if two edge iterators are equal. Two edge iterators are equal if 
     *  they have the same current position at the (@a graph).(@a edges)_ vector
     *  in the same graph.
     */
    bool operator==(const EdgeIterator& ei) const{
      return (graph_  == ei.graph_  && 
              curr_n_ == ei.curr_n_ && 
              curr_e_ == ei.curr_e_);
    }

   private:
    friend class Graph;

    Graph* graph_; // Pointer to graph
    size_type curr_n_; // Current node
    size_type curr_e_; // Current edge

    void fix(){
      bool end = false;
      while (!end){
        // If you've reached the end of i2u_, set curr_e_ to 0 and stop.
        if (curr_n_ == graph_ -> size()){
          curr_e_ = 0;
          return;
        }
        // If you've reached the end of the adj_ vector of curr_n_, 
        // move to the first position of the next node.
        else if (curr_e_ == 
                 graph_ -> nodes_[graph_ -> i2u_[curr_n_]].adj_.size()){
          ++curr_n_;
          curr_e_ = 0;
        }
        // Guarantee specific node ordering to prevent duplicate edges
        else if (graph_ -> i2u_[curr_n_] > 
                 graph_ -> 
                 nodes_[graph_ -> i2u_[curr_n_]].adj_[curr_e_].uid_other_){
          ++curr_e_; // Skip this edge
        }
        else{
          end = true; // you're done fixing!
        }
      }
    }

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type curr_n, size_type curr_e)
    : graph_(const_cast<Graph*>(graph)), curr_n_(curr_n), curr_e_(curr_e){
      fix();
    }
  };
  // Returns first position at the adj_ vector of the first node
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0, 0);
  }

  // Returns one position past the range of the abstract end.
  edge_iterator edge_end() const{
    return EdgeIterator(this, size(), 0);
  }

 private:
  
  /* Internal type for nodes. */
  struct node_info{
    // Number in [0,size()) if node is valid or -1 if node is invalid. 
    idx_type idx_; 
    // Cartesian location. 
    Point position_; 
    // Value of the node. 
    node_value_type value_; 
    // Adjacency vector containing edge value and uid_ of other incident node.
    std::vector<edge_info> adj_;  
  };

  /* Internal type for edge. */
  struct edge_info{
    // Unique id of other node incident on this edge.
    uid_type uid_other_;
    // Value of type edge_value_type of this edge.
    edge_value_type value_;
  };

  // count of the total number of edges
  size_type count_edges = 0;

  // STD container that holds the nodes' data:
  // (unique id, index, position, value, and adjacency vector
  std::vector<node_info> nodes_;

  // Mapping from index to unique uid
  std::vector<uid_type> i2u_;

};

#endif // CME212_GRAPH_HPP
