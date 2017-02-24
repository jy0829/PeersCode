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
template<typename V, typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;
 

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

  using edge_value_type = E;
  //
  // CONSTRUCTORS AND DESTRUCTORc
  //

  /** Construct an empty graph. */
  Graph(): nodes_(), valid_nodes_(), num_edges_(0), edges_() {
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
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[node_uid_].point_;
    }

    Point& position() {
      return graph_->nodes_[node_uid_].point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[node_uid_].index_;
    }
    /** Return this node's value, the value can be changed **/
    node_value_type& value(){
      return graph_->nodes_[node_uid_].value_;
    }
    /** Return this node's value, the value cannot be changeds **/
    const node_value_type& value() const{
      return graph_->nodes_[node_uid_].value_;
    }


    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // const node_value_type& value() const;

    // return the number of the edges that incident to the node
    size_type degree() const{
      return graph_->connect_edges_[node_uid_].size();
    }

    // return the iterator that points to the start of the node's incident edge
     incident_iterator edge_begin() const{
      return incident_iterator(graph_, node_uid_, 0);

     }

    // return the iterator that points to the end of the node's incident edge
    incident_iterator edge_end() const{
      return incident_iterator(graph_, node_uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //(void) n;          // Quiet compiler warning
      return (graph_ == n.graph_ && node_uid_ == n.node_uid_);
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
      //(void) n; 
      if (graph_ < n.graph_) {
        return true;
      } 
      else if (graph_ == n.graph_ && node_uid_ < n.node_uid_){
        return true;
      } 
      //else{
      //  return false;
      //}       
      return false;
    }



   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    /* pointer to the graph */
    graph_type *graph_; 
    /* unique indentification number for the node */
    size_type node_uid_; 
    /* private constructor */
    Node(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), node_uid_(uid){

      }

    /** Helper method to return the appropriate element 
     */
    internal_node& fetch() const {
        return graph_->nodes_[node_uid_];
        assert(false);
    }


  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return valid_nodes_.size();
    //return 0;
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
   * Complexity: O(1) amortized operations.c
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    internal_node new_internal_node;
    new_internal_node.point_ = position;
    new_internal_node.index_ = valid_nodes_.size();
    new_internal_node.value_ = value;
    nodes_.push_back(new_internal_node);
    valid_nodes_.push_back(nodes_.size()-1);
    connect_edges_.resize(nodes_.size());
    connect_edges_values_.resize(nodes_.size());
    //(void) position;      // Quiet compiler warning
    return Node(this, nodes_.size()-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    return (this == n.graph_);
    //(void) n;            // Quiet compiler warning
    //return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //assert(i >= 0 && i < num_nodes());

    //(void) i;             // Quiet compiler warning
    return Node(this, valid_nodes_[i]);        // Invalid node
  }

  // remove node
  size_type remove_node(const Node &n){
    if (!has_node(n)){
      return 0;
    }
    size_type temp = 0;
    // find the node
    for (size_type i = 0; i < valid_nodes_.size(); i++){
      if (n.node_uid_ == valid_nodes_[i]){
        temp = i;
        break;
      }
    }
    valid_nodes_.erase(valid_nodes_.begin()+temp);
    for (size_type i = temp; i< valid_nodes_.size(); i++){
      internal_node& in = nodes_[valid_nodes_[i]];
      in.index_ = i;
    }
    // remove all incident edges
    auto in_edge = connect_edges_[n.node_uid_];
    for (auto ie : in_edge){
      remove_edge(node(nodes_[ie].index_), n);
    }
    return temp;
  }

  node_iterator remove_node(node_iterator n_it){
    if (n_it != node_end()){
      remove_node(*n_it);
    }
    return n_it;
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
      return Node(graph_, node_uid1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node_uid2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //return false;
      return (e.graph_ == graph_ && node1() == e.node1() && node2() == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //return false;
      if (node_uid1_ < e.node_uid1_){
        return true;
      } 
      else if (node_uid1_ == e.node_uid1_ && node_uid2_ < e.node_uid2_){
        return true;
      } 
      
      return false;
    }

    // calculate the length of the edge
    double length() const{
      return norm(node1().position() - node2().position());

    }

    edge_value_type &value(){
      size_type node_min = std::min(node_uid1_, node_uid2_);
      size_type node_max = std::max(node_uid1_, node_uid2_);
      Edge e = Edge(graph_, node_min, node_max);
      for (auto it = e.node1().edge_begin(); it != e.node1().edge_end(); ++it){
        if ((*it).node_uid2_ == e.node_uid2_){
          return graph_->connect_edges_values_[e.node_uid1_][it.node2_list_id_];
        }
      }
    }

    const edge_value_type &value() const{
      size_type node_min = std::min(node_uid1_, node_uid2_);
      size_type node_max = std::max(node_uid1_, node_uid2_);
      Edge e = Edge(graph_, node_min, node_max);
      for (auto it = e.node1().edge_begin(); it != e.node1().edge_end(); ++it){
        if ((*it).node_uid2_ == e.node_uid2_){
          return graph_->connect_edges_values_[e.node_uid1_][it.node2_list_id_];
        }
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type node_uid1_;  // id for the first node
    size_type node_uid2_;  // id for the second node
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const graph_type* graph, size_type node_uid1, size_type node_uid2)
      : graph_(const_cast<graph_type*> (graph)), node_uid1_(node_uid1), node_uid2_(node_uid2){

      } 

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;  
    assert(i >= 0 && i < this->num_edges());           // Quiet compiler warning
    size_type node_index1 = edges_[i][0];
    size_type node_index2 = edges_[i][1];
    return Edge(this, node_index1, node_index2);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //(void) a; (void) b;   // Quiet compiler warning
    for (size_type i = 0; i < connect_edges_[a.node_uid_].size(); i++){
      if (connect_edges_[a.node_uid_][i] == b.node_uid_){
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
    //(void) a, (void) b;   // Quiet compiler warning
    
    assert(a.node_uid_ != b.node_uid_);
    if (!has_edge(a, b)){
      connect_edges_[a.node_uid_].push_back(b.node_uid_);
      connect_edges_[b.node_uid_].push_back(a.node_uid_);
      connect_edges_values_[a.node_uid_].push_back(edge_value_type());
      connect_edges_values_[b.node_uid_].push_back(edge_value_type());
      std::vector<size_type> edge_temp;
      edge_temp.push_back(a.node_uid_);
      edge_temp.push_back(b.node_uid_);
      edges_.push_back(edge_temp);
      num_edges_++;
    }
    return Edge(this, a.node_uid_, b.node_uid_);        // Invalid Edge
    
    //return Edge();
  }

  size_type remove_edge(const Node& a, const Node& b){
    if (!has_edge(a, b)){
      return a.degree();
    }
    size_type res = a.degree();
    for (size_type i = 0; i < connect_edges_[a.node_uid_].size(); i++){
      if (connect_edges_[a.node_uid_][i] == b.node_uid_ ){
        connect_edges_[a.node_uid_].erase(connect_edges_[a.node_uid_].begin()+i);
        connect_edges_values_[a.node_uid_].erase(connect_edges_values_[a.node_uid_].begin()+i);
        res = i;
        break;
      }
    }

    for (size_type i = 0; i < connect_edges_[a.node_uid_].size(); i++){
      if (connect_edges_[b.node_uid_][i] == a.node_uid_ ){
        connect_edges_[b.node_uid_].erase(connect_edges_[b.node_uid_].begin()+i);
        connect_edges_values_[b.node_uid_].erase(connect_edges_values_[b.node_uid_].begin()+i);
        break;
      }
    }
    num_edges_--;
    return res;


  }

  size_type remove_edge(const Edge &e){
    return remove_edge(e.node1(), e.node2());
  }

  edge_iterator remove_edge(edge_iterator e_it){
    Edge e = *e_it;
    size_type idx = remove_edge(e);
    return edge_iterator(this, e.node1().node_uid_, idx);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    valid_nodes_.clear();
    connect_edges_.clear();
    edges_.clear();
    num_edges_ = 0;

  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
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
    // Deference NodeIterator and return current node 
    // @pre @a idx_ should be valid(>0) 
    // @post return an instance of node corresponding to the current iterator 
    Node operator*() const{
      return graph_->node(idx_);
    }

    // increment and return the next NodeIterator
    // @pre @a idx_ < graph_->size() && @a idx_>=0
    // @post return the next valid iterator 
    NodeIterator& operator++(){
      idx_++;
      return (*this);

    }

    // check whether two NodeIterators are equivalent
    // @param @a nodeIterator NodeIterator
    // @post return true when @a nodeIterator and current NodeIterator 
    // points to the same graph and the same node, otherwise return false */
    bool operator==(const NodeIterator& node_iterator) const{
      return ( (graph_ == node_iterator.graph_) && (idx_ == node_iterator.idx_) );
    }

   private:
    friend class Graph;
    graph_type *graph_; 
    size_type idx_;
    NodeIterator(const graph_type* graph, size_type idx):
      graph_(const_cast<graph_type*>(graph)), idx_(idx){

      }
    // HW1 #2: YOUR CODE HERE
  };


  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // return the node_iterator that points to the first node of the graph
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  // return the node_iterator that points to the last node of the graph
  node_iterator node_end() const{
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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

    // deference Incident Iterator and return current edge information
    Edge operator*() const{
      return Edge(graph_, node1_id_, graph_->connect_edges_[node1_id_][node2_list_id_]);
    }

    // increment current iterator and return the next position
    // @pre @a node2_list_id_  < connect_edges[node1_id_].size()-1
    // return incident iterator points to the next
    IncidentIterator& operator++(){
      node2_list_id_++;
      return  (*this);
    }

    // check wheter two IncidentIterator are equivalent 
    // @return true the they have are in the same graph and same incident edge 
    // @return otherwise return false
    bool operator==(const IncidentIterator& x) const{
      return ((graph_ == x.graph_) && (node1_id_ == x.node1_id_) && (node2_list_id_ == x.node2_list_id_));

    }

   private:
    friend class Graph;

    Graph* graph_;
    // node id for the first node
    size_type node1_id_;

    // the index in the vector that correspond to the second node
    size_type node2_list_id_;
    // HW1 #3: YOUR CODE HERE

    IncidentIterator(const Graph* graph, size_type node1_id, size_type node2_list_id)
      : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_list_id_(node2_list_id){

      }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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

    // deference edge iterator and return current inforamtion
    Edge operator*() const{
      assert(graph_ != nullptr && node1_id_ < graph_->connect_edges_.size());
      return Edge(graph_, node1_id_, graph_->connect_edges_[node1_id_][node2_list_id_]);
    }

    // increment edge iterator and return next edge information
    EdgeIterator& operator++(){
      if (graph_->connect_edges_[node1_id_].size() == 0){
        node1_id_++;
        node2_list_id_ = 0;
      }
      else if (node2_list_id_ >= graph_->connect_edges_[node1_id_].size()-1 ){
        node1_id_++;
        node2_list_id_ = 0;
      }
      else{
        node2_list_id_++;
      }
      return (*this);
    }
    // check wheter two edgeiterators are equivalent
    bool operator==(const EdgeIterator& x) const{
      return (graph_ == x.graph_ && node1_id_ == x.node1_id_ && node2_list_id_ == x.node2_list_id_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    // node id for the first node
    size_type node1_id_;

    // the index in the vector that correspond to the second node
    size_type node2_list_id_;

    EdgeIterator(const Graph* graph, size_type node1_id, size_type node2_list_id)
      : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_list_id_(node2_list_id){

      }

    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // return the edge_iterator points to the start of the edge
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0, 0);
  }

  // return the edge_iterator points to the end of the edge
  edge_iterator edge_end() const{
    return EdgeIterator(this, this->connect_edges_.size(), 0 );
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node{
      Point point_;
      size_type index_;
      node_value_type value_;
  };


  // Container stores all the nodes
  std::vector<internal_node> nodes_;

  // Container stores all the valid node's id
  std::vector<size_type> valid_nodes_;

  // number of edges
  size_type num_edges_;

  /* Container stores all connected edge 
   * The index of the outer vector represents the @a node_uid_ of the first node
   * The internal std::vector stores the @a node_uid_ of the second node
   */

  std::vector<std::vector<size_type>> connect_edges_;

  /* Container stores all connected edge 
   * The index of the outer vector represents the @a node_uid_ of the first node
   * The internal std::vector stores the @a edge_value of the second node
   */

  std::vector<std::vector<edge_value_type>> connect_edges_values_;
  /* Container stores all the edge information
   * The index of the outer vector represents the index of the edge
   * The internal std::vector stores only two elements, 
   * the first element: the @a node_uid_ of the first node
   * the second elemetn: the @a node_uid_ of the second node
   */ 
  std::vector<std::vector<size_type>> edges_;


};

#endif // CME212_GRAPH_HPP
