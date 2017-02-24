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
  //  typedef unsigned size_type;
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

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
    :nodes_(std::map<size_type,std::pair<Point,node_value_type> > ()),
    adjacent_nodes_(std::map<size_type, std::vector<size_type> > ()){
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
  class Node : private totally_ordered <Node> {
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
    Node() :graph_(0),
            node_id_(0) {
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(this->graph_ != NULL);
      return (graph_->nodes_.at(node_id_).first);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(this->graph_ != NULL);
      assert(graph_->size() != 0 and graph_->size() > node_id_);
      return node_id_;
    }

    /** Get the non-const value held by this node.
    */
    node_value_type& value(){
      assert(this->graph_ != NULL);
      return graph_->nodes_.at(node_id_).second;
    }
    
    /** Get the value held by this node.
    */
    const node_value_type& value() const{
      assert(this->graph_ != NULL);
      return graph_->nodes_.at(node_id_).second;
    }
    
    /** Get the number of nodes incident on current node.
    */
    size_type degree() const{
      return graph_->adjacent_nodes_.at(node_id_).size();
    }
    
    /** Get the first edge for iterator begin
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, node_id_, 
        graph_->adjacent_nodes_.at(node_id_).at(0));
    }

    /** Get the last edge for iterator end
    */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, node_id_, degree());

    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if(graph_ == n.graph_ && node_id_ == n.node_id_){
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
      assert(n.graph_ != NULL);
      if(graph_ == n.graph_ && node_id_ < n.node_id_) return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type node_id_;
    
    Node(const Graph* graph, size_type node_id)
         : graph_(const_cast<graph_type*>(graph)),
          node_id_(node_id) {
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1`
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    size_type index_p = size();
    //nodes_[index_p]= std::make_pair(const_cast<Point>(&position),
    nodes_[index_p]= std::make_pair(position,
      node_value_type());
    /** Update the attributes of graph. */
    return Node(this, index_p);
  }

  Node add_node ( const Point & position, const node_value_type & node_value){
    size_type index_p = size();
    nodes_[index_p]= std::make_pair(&position,node_value);
    adjacent_nodes_[index_p] = std::vector<size_type> ();
    /** Update the attributes of graph. */
    return Node(this, index_p);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    assert(n.graph_ != NULL);
    if(n.graph_ == this and n.node_id_ < size()){
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
    assert(i<size());
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
  class Edge : private totally_ordered <Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge(): graph_(0),
            node1_id_(0),
            node2_id_(0){
    }

    /** Return a node of this Edge */
    Node node1() const {
      assert(this->graph_ != NULL);
      return Node(graph_,node1_id_);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(this->graph_ != NULL);
      return Node(graph_,node2_id_);   
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(e.graph_ != NULL && this->graph_ != NULL);
      if(graph_ == e.graph_ &&
          ((node1_id_ == e.node1_id_ && node2_id_ == e.node2_id_) ||
          (node2_id_ == e.node1_id_ && node1_id_ == e.node2_id_))) 
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(e.graph_ != NULL && this->graph_ != NULL);
      size_type min_id = std::min(node1_id_,node2_id_);
      size_type max_id = std::max(node1_id_,node2_id_);
      size_type e_min_id = std::min(e.node1_id_,e.node2_id_);
      size_type e_max_id = std::max(e.node1_id_,e.node2_id_);
      if(graph_ == e.graph_ &&
          ((min_id < e_min_id) ||
          (min_id == e_min_id && max_id < e_max_id))) 
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type node1_id_;
    size_type node2_id_;

    Edge(const Graph* graph, size_type node1_id, size_type node2_id) 
         : graph_(const_cast<graph_type*>(graph)),
          node1_id_(node1_id),
          node2_id_ (node2_id){
    }

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity requested: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type num_edges = 0;
    for(auto& i : adjacent_nodes_){
      num_edges += i.second.size();
    }
    num_edges /= 2;
    return num_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity requested: No more than O(num_nodes() + num_edges()), hopefully less
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    edge_iterator iter = edge_begin();
    size_type count =0;
    while(count!=i || iter!= edge_end()) {
      ++iter;
      ++count;
    }
    if(iter!= edge_end())
      return *iter;
    else
      return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(has_node(a) && has_node(b) && adjacent_nodes_.size() > a.index())
      if( (adjacent_nodes_.find(a.index()) != adjacent_nodes_.end())&&
       std::find(adjacent_nodes_.at(a.index()).begin(),
                    adjacent_nodes_.at(a.index()).end(),
                    b.index()) != adjacent_nodes_.at(a.index()).end()){
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
  Edge add_edge(const Node& a, const Node& b) {
    size_type node1_id = a.index();
    size_type node2_id = b.index();
    if(!has_edge(a,b)){
      if(!has_node(a))
        adjacent_nodes_.at(node1_id) = std::vector<size_type> ();
      if(!has_node(b))
        adjacent_nodes_.at(node2_id) = std::vector<size_type> ();
      adjacent_nodes_[node1_id].push_back(node2_id);
      adjacent_nodes_[node2_id].push_back(node1_id);
    }
    return Edge(this,node1_id,node2_id);        
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adjacent_nodes_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered <NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_iter_(0),
      node_iter_id_(0) {
    }

    Node operator*() const{
      return graph_iter_->node(node_iter_id_);
    }

    NodeIterator& operator++(){
      size_type next_id = ++node_iter_id_;
      if(next_id>graph_iter_->size())
        next_id = graph_iter_->size();
      node_iter_id_ = next_id;
      return *this;
    }

    bool operator==(const NodeIterator& node_iter) const{
      assert(graph_iter_!=NULL);
      return node_iter.node_iter_id_ == node_iter_id_;
    }

   private:
    friend class Graph;
    graph_type * graph_iter_;
    size_type node_iter_id_;

    NodeIterator(const graph_type* graph, size_type node_id)
      : graph_iter_(const_cast<graph_type*>(graph)),
        node_iter_id_(node_id){
    }

  };

  node_iterator node_begin() const{
    return NodeIterator(this,(this->nodes_).begin()->first);
  }

  node_iterator node_end() const{
    return NodeIterator(this,size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered <IncidentIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_(0),
        node_(0),
        adj_node_idx_(0){
    }
    
    Edge operator*() const{
      return Edge(graph_, node_, 
        graph_->adjacent_nodes_.at(node_).at(adj_node_idx_));
    }
    
    IncidentIterator& operator++(){
      ++adj_node_idx_;
      return *this;
    }
    
    bool operator==(const IncidentIterator& input) const{
      return (node_ == input.node_ && 
        adj_node_idx_ == input.adj_node_idx_);
    }

   private:
    friend class Graph;

    /** node1 of edge == node_ and node2 == adj_node_
    */
    graph_type* graph_;
    size_type node_;
    unsigned adj_node_idx_;

    IncidentIterator(const graph_type * graph, size_type node1, unsigned node2_idx)
      : graph_(const_cast<graph_type*>(graph)),
        node_(node1),
        adj_node_idx_(node2_idx){
        }
  };

  incident_iterator incident_edge_begin(Node& a) const {
    size_type idx = a.index();
    if((this->adjacent_nodes_).at(idx).size()>0)
      return IncidentIterator(this, idx, 0);
    else
      return IncidentIterator();
  }

  incident_iterator incident_edge_end(Node& a) const{
    size_type idx = a.index();
    if((this->adjacent_nodes_).at(idx).size()>0)
      return IncidentIterator(this, idx, (this->adjacent_nodes_).at(idx).size());
    else
      return IncidentIterator();
   }


  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered <EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(0),
        node_(0),
        adj_node_idx_(0){
    }

    //returns edge corresponsing to the iterator
    Edge operator*() const{
      return Edge(graph_, node_, 
        graph_->adjacent_nodes_.at(node_).at(adj_node_idx_));
    }
    
    /** Increment node2 index till node1 has corresponding node2
    * then increment node1 and look for valid edge
    * if nothing can be found then return edge_iterator end values. 
    */
    EdgeIterator& operator++(){
      ++adj_node_idx_;
      check();
      return *this;
    }

     /* check not to access same edge twice by noting that if 
      * node2 id is less than node1 id then it is already accessed.
      * and increment to next position. 
      */
    void check(){
      while( node_< graph_->size() &&
             adj_node_idx_>=graph_->adjacent_nodes_.at(node_).size() ){
        adj_node_idx_=0;
        ++node_;
        if(node_ == graph_->size())
          break;
        while(node_<graph_->size() && 
          graph_->adjacent_nodes_.at(node_).size()==0  )
            ++node_;
        while(adj_node_idx_<graph_->adjacent_nodes_.at(node_).size() && 
              graph_->adjacent_nodes_.at(node_).at(adj_node_idx_) < node_)
          ++adj_node_idx_;
      }
    }
    
    bool operator==(const EdgeIterator& input) const{
      return (node_ == input.node_ && adj_node_idx_ == input.adj_node_idx_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type node_;
    size_type adj_node_idx_;

    EdgeIterator(const graph_type * graph, size_type node1, size_type node2_idx)
      : graph_(const_cast<graph_type*>(graph)),
        node_(node1),
        adj_node_idx_(node2_idx){}
  };

  //find the first valid edge
  edge_iterator edge_begin() const {
    auto iter = (this->nodes_).begin();
    size_type idx = iter->first;
    while(iter!=(this->nodes_).end() &&
      (this->adjacent_nodes_).at(idx).size()==0) 
      ++iter;
    return EdgeIterator(this, iter->first, 0);
  }

  //keep last edge as one with node1 out of bound.
  edge_iterator edge_end() const{
    return EdgeIterator(this, size(),0);
  }

 private:

    /**when Point* is used in pair the SFML_viewer add_nodes needs to be edited
    * so stuck with Point, node_value_type. TODO: Any benefit? 
    */
    std::map<size_type, std::pair<Point, node_value_type> > nodes_;
    std::map<size_type, std::vector<size_type> > adjacent_nodes_;

};

#endif // CME212_GRAPH_HPP

