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

template <typename V>
class Graph {
 private:

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
  //  CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), edges_(), A_() {
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
      return graph_-> nodes_[index_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /** Return this node's value, of type node_value_type. */
    node_value_type& value(){
      return graph_-> nodes_[index_].second;
    }

    /** Return this node's value (overriding for correctness).*/
    const node_value_type& value() const{
      return graph_-> nodes_[index_].second;
    }
 
    /** Return this node's degree, a number in the range [0, graph_size-1).
     *  Recall that the degree() of a node is the number of edges incident
     *  to it.
     *
     *  @result: Number in [0, @a num_edges()) 
     */
    size_type degree() const{
      return graph_ -> A_[index_].size();
    }

    /** Returns first position of data structure containing the adjacency 
     *  information of node @a index_
     *
     *  @result: Position of the first element of the vector A_[index_].
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, index_, 0);
    }

    /** Returns last position of data structure containing the adjacency      
     *  information of node @a index_.
     *
     *  @result: One position past the abstract range of the vector 
     *           @a A_[@a index_].
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *  Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (index_ == n.index() && graph_ == n.graph_);
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
      return (index_ < n.index() && graph_ == n.graph_);
    }

   private:
    
    Graph* graph_; // Pointer back to the Graph container
    size_type index_; // This node's index

    /** Private Constructor 
     * 
     *  The job of the private constructor of node is to establish and 
     *  defend the representation invariance of Graph::Node. If the Graph 
     *  ever tries to construct an invalid node, it will be detected by this 
     *  constructor.
     *
     */
    Node(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index){
      assert(graph_ != nullptr);
      assert(index_ < graph_ -> num_nodes());
    }

    // Allow Graph to access Node's private member data and functions.*/
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   *  Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
    std::pair<Point, node_value_type> n(position,value);
    nodes_.push_back(n); // Add new node to nodes_ vector
    std::vector<size_type> a;
    A_.push_back(a); // Add adjacency information vector for new node
    return Node(this,size()-1);   
  }

  /** Determine if a Node belongs to this Graph
   *  @return True if @a n is currently a Node of this Graph
   *
   *  Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return  (n.graph_ == this && n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   *  @pre 0 <= @a i < num_nodes()
   *  @post result_node.index() == i
   *
   *  Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);        
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

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(id_a_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(id_b_);
    }

    /** Test whether this edge and @a e are equal.
     *
     *  Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_ &&
	      std::min(id_a_,id_b_) == std::min(e.id_a_, e.id_b_) &&
              std::max(id_a_,id_b_) == std::max(e.id_a_, e.id_b_)   );
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     *  This ordering function is useful for STL containers such as
     *  std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (graph_ == e.graph_ &&
              (std::min(id_a_,id_b_) < std::min(e.id_a_,e.id_b_)  ||
              (std::min(id_a_,id_b_) == std::min(e.id_a_,e.id_b_) &&
               std::max(id_a_,id_b_)  < std::max(e.id_a_,e.id_b_))));
    }

   private: 
    Graph* graph_;//  Pointer back to the Graph container
    size_type id_a_; //  This edge's index a
    size_type id_b_; //  This edge's index b

    /** Private Constructor 
     * 
     *  The job of the private constructor of edge is to establish and 
     *  defend the representation invariance of Graph::Edge. If the Graph 
     *  ever tries to construct an invalid edge, it will be detected by this 
     *  constructor.
     *
     */
    Edge(const Graph* graph, size_type id_a, size_type id_b)
      : graph_(const_cast<Graph*>(graph)), id_a_(id_a), id_b_(id_b){
      assert(graph != nullptr);
      assert(id_a_ != id_b_);
      assert(id_a_ < graph -> num_nodes());
      assert(id_b_ < graph -> num_nodes());
    }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
  };

  /** Return the total number of edges in the graph.
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   *  @pre 0 <= @a i < num_edges()
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { 
    return Edge(this,edges_[i].first,edges_[i].second);  
  }

  /** Test whether two nodes are connected by an edge.
   *  @pre @a a and @a b are valid nodes of this graph
   *  @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for(auto ei = a.edge_begin(); ei != a.edge_end(); ++ei){
      edge_type e = *ei;
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
  Edge add_edge(const Node& a, const Node& b) {
    if(not has_edge(a,b)){
      std::pair<size_type, size_type> e (a.index(),b.index());
      edges_.push_back(e); // Add new edge to edges_ container
      // Add adjacency information to nodes incident to new edge.
      A_[a.index()].push_back(b.index()); 
      A_[b.index()].push_back(a.index());
    }
    return Edge(this,a.index(),b.index());
  }

  /** Remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *
   *  Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    A_.clear();
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
     *  Access node at the current position of the nodes_ vector.
     */ 
    Node operator*() const{
      return Node(graph_,curr_);
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
    size_type curr_; // Current position in @a nodes_ vector. 

    /** Private Constructor 
     * 
     *  The job of the private constructor of NodeIterator is to establish and 
     *  defend the representation invariance of Graph::NodeIterator. If the Graph 
     *  ever tries to construct an invalid node iterator, it will be detected by
     *  this constructor.
     *
     */
    NodeIterator(const Graph* graph, size_type curr)
      : graph_(const_cast<Graph*>(graph)), curr_(curr){
      assert(graph_ != nullptr);
      assert(curr_ <= graph_ -> num_nodes());
    }
  };

  // Returns first position of @a nodes_ vector
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  // Returns last position of @a nodes_ vector
  node_iterator node_end() const{
    return NodeIterator(this, num_nodes());
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
     *  Access edge that connects node(@a index_) with the node at the current
     *  position of @a A_[@a node_id_] vector.
     *  This vector contains the adjacency information of node(@a index_).
     */ 
    Edge operator*() const{
      return Edge(graph_, node_id_, graph_->A_[node_id_][curr_]);
    }
    /** Increment operator:
     *  Increment the current position at @a A_[@a node_id_] by one.
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
    size_type node_id_; // Node you're iterating through
    size_type curr_; // Current position in @a A_[@a node_id_]

    /** Private Constructor */
    IncidentIterator(const Graph* graph, size_type node_id, size_type curr)
        : graph_(const_cast<Graph*>(graph)), node_id_(node_id), curr_(curr){
      assert(graph_ != nullptr);
      assert(node_id_ <= graph_ -> num_nodes());
      assert(curr_ <= graph_ -> A_[node_id_].size());
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
      return Edge(graph_, graph_ -> edges_[curr_].first, 
                          graph_ -> edges_[curr_].second);
    }
    /** Increment operator:
     *  Increment the current position at @a edges_ by one.
     */
    EdgeIterator& operator++(){
      ++curr_;
      return *this;
    }
    /** Equality operator:
     *  Check if two edge iterators are equal. Two edge iterators are equal if 
     *  they have the same current position at the (@a graph).(@a edges)_ vector
     *  in the same graph.
     */
    bool operator==(const EdgeIterator& ei) const{
      return (curr_ == ei.curr_ && graph_ == ei.graph_);
    }
   private:
    friend class Graph;
    Graph* graph_; // Pointer to graph
    size_type curr_; // current position at (@a graph).(@a edges_) vector.
    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type curr)
      : graph_(const_cast<Graph*>(graph)), curr_(curr){
    }
  };
  // Returns first position of @a edges_ vector
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }
  // Returns last position of @a edges_ vector (one past range of abstract end)
  edge_iterator edge_end() const{
    return EdgeIterator(this, num_edges());
  }

 private:

  // STD container that holds the nodes' data (unique ID and position)
  std::vector<std::pair<Point, node_value_type>> nodes_;

  // STD container that holds the edges' data (IDs of connected nodes)
  std::vector<std::pair<size_type, size_type>> edges_;
  
  // STD container that holds the adjacency info of the graph
  std::vector<std::vector<size_type>> A_;

};

#endif // CME212_GRAPH_HPP
