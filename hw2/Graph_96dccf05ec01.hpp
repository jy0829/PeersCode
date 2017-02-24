#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <tuple>
#include <vector>

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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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

  /** Types for edge and node values*/
  using node_value_type = V;
  using edge_value_type = E;

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
  Graph() : nodes(),i2u_(){
  }

  /** Default destructor */
  ~Graph() = default;

  /** Class for the node information which will be stored in the nodes array of the graph*/
  struct node_info{
    Point position_;  //< Position of the node
    node_value_type value_;  //< Value of the node
    std::vector<std::pair<size_type,edge_value_type>> adjacency_;  //< vector with adjacent
                                                                   //nodes and edges
    size_type idx_;  //< Node index
    node_info(Point pos, node_value_type value, std::vector<std::pair<size_type,
              edge_value_type>> adjacency, size_type idx):position_(pos),value_(value),
              adjacency_(adjacency), idx_(idx){}
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
     */

    Node() {
    }

    /** Return this node's position. */
    Point& position() const {
      return graph_->nodes[uid_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes[uid_].idx_;
    }

    /** Return a reference to the value of a node.
    * @post result == reference to the value of the node.
    */
    node_value_type& value(){
      return graph_->nodes[uid_].value_;
    }
    /** Return the value of a node.
    * @post result == value of the node.
    */
    const node_value_type& value() const{
      return graph_->nodes[uid_].value_;
    }

    /** Return the degree of a node.
    * @post result == n<N. Where n is the number of adjacent nodes to the current
             node, and N the total number of nodes in the graph.
    */
    size_type degree() const{
    return graph_->nodes[uid_].adjacency_.size();
    }

    /** Return the first incident iterator.
    * @post result <= a. Where a any other valid incident iterator for the node.
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,uid_,0);
    }

    /** Return the last incident iterator.
    * @post result >= a. Where a is any other valid incident iterator for the node.
    */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_,uid_,this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (graph_ == n.graph_)
        if (uid_ ==  n.uid_)
          return true;
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
      if (graph_ == n.graph_){
        if (uid_ <  n.uid_)
          return true;
      }
      else{
        if (graph_ < n.graph_)
          return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph_pointer, size_type uid)
        : graph_(const_cast<Graph*>(graph_pointer)), uid_(uid) {
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
    return i2u_.size();
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
  Node add_node(const Point& position, const node_value_type& value= node_value_type()) {
    std::vector<std::pair<size_type, edge_value_type>> adjacency_vector;
    nodes.push_back(node_info(position,value,adjacency_vector,i2u_.size()));
    i2u_.push_back(nodes.size()-1);
    return Node(this,nodes.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this)
      if (n.graph_->nodes[n.uid_].idx_ < i2u_.size())
        return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u_[i]);
  }

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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_,node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,node2_id);
    }

    edge_value_type& value(){
      unsigned index = 0;
      for (unsigned i=0; i < graph_->nodes[node1_id].adjacency_.size(); ++i){
        if (graph_->nodes[node1_id].adjacency_[i].first == node2_id){
          index = i;
          break;
        }
      }
      return graph_->nodes[node1_id].adjacency_[index].second;
    }

    const edge_value_type& value() const{
      unsigned index = 0;
      for (unsigned i=0; i < graph_->nodes[node1_id].adjacency_.size(); ++i){
        if (graph_->nodes[node1_id].adjacency_[i].first == node2_id){
           index = i;
           break;
        }
      }
      return graph_->nodes[node1_id].adjacency_[index].second;
    }

    double length() const{
      Node n1 = node1();
      Node n2 = node2();
      return norm_2(n1.position()-n2.position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      auto e1 = std::make_tuple(graph_,std::min(node1_id,node2_id),std::max(node1_id,node2_id));
      auto e2 = std::make_tuple(e.graph_,std::min(e.node1_id,e.node2_id),std::max(e.node1_id, e.node2_id));
      return e1==e2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      auto e1 = std::make_tuple(graph_,std::min(node1_id,node2_id),std::max(node1_id,node2_id));
      auto e2 = std::make_tuple(e.graph_,std::min(e.node1_id,e.node2_id),std::max(e.node1_id, e.node2_id));
      return e1<e2;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type node1_id;
    size_type node2_id;
    Edge(const Graph* graph_pointer, size_type node1, size_type node2)
        :graph_(const_cast<Graph*>(graph_pointer)), node1_id(node1), node2_id(node2){
    }
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return std::distance(this->edge_begin(),this->edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    size_type count = 0;
    for (auto it = this->edge_begin(); it != this->edge_end(); ++it){
      if (count==i)
        return *it;
      ++count;
    }
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned int i=0; i<nodes[a.uid_].adjacency_.size(); i++)
       if (b.uid_ == nodes[a.uid_].adjacency_[i].first)
         return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value= edge_value_type()) {
    if (!has_edge(a,b)){
      nodes[a.uid_].adjacency_.push_back(std::make_pair(b.uid_,value));
      nodes[b.uid_].adjacency_.push_back(std::make_pair(a.uid_,value));
    }
    return Edge(this,a.uid_,b.uid_);
  }

  // REMOVE FUNCTIONS

  /** Remove a node and the incident edges.
   * @param[in] n  Node to be removed
   * @pre 0<=n.index()<g.num_nodes()
   * @return n.index().
   * @post new num_nodes() == old num_nodes()-1
   * @post If node @a n has incident edges, new num_edges() < old num_edges()
   * The node is invalidated by setting its index to -1. The incident edges
   * are removed and all their information.
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_node(const Node& n){
    // First delete edges in the adjacency matrix
    auto ei = n.edge_begin();
    while (ei!=n.edge_end()){
      Edge e = *ei;
      Node n2 = e.node2();
      size_type removed=remove_edge(n,n2);
      (void) removed;
    }
    // Delete node from i2u_ and invalidate node in nodes (idx_ = -1)
    unsigned index_to_delete = nodes[n.uid_].idx_;
    i2u_.erase(i2u_.begin()+index_to_delete);
    nodes[n.uid_].idx_ = -1;
    // Decrease the index by one of nodes with higher index
    for (auto ni = this->node_begin(); ni!= this->node_end(); ++ni){
      Node n2 = *ni;
      if (n2.uid_ > n.uid_)
        nodes[n2.uid_].idx_ -= 1;
    }
    return index_to_delete;
  }

  /** Remove a node and the incident edges.
   * @param[in] n_it  Iterator pointing to node to be removed
   * @pre 0<=*n_it.index()<g.num_nodes()
   * @return n_it.
   * @post n_it points to the next node. If there is
   *       no next node, n_it==edge_end().
   * @post new num_nodes() == old num_nodes()-1
   * @post If node *n_it has incident edges, new num_edges() < old num_edges()
   * The node is invalidated by setting its index to -1. The incident edges
   * are removed and all their information. Iterators [n_it, last) are invalidated
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  node_iterator remove_node(node_iterator n_it){
    Node n = *n_it;
    size_type node_deleted = remove_node(n);
    (void) node_deleted;
    return n_it;
  }

   /** Remove an edge.
   * @param[in] n1  First node of edge to be removed
   * @param[in] n2  Second node of edge to be removed
   * @pre has_edge(@a n1,@a n2)==true
   * @return 0
   * @post new num_edges() == old num_edges()-1.
   * The edge is removed and all its information
   * Complexity: No more than O(num_edges())
   */
  size_type remove_edge(const Node& n1, const Node& n2){
    for (unsigned i=0; i < nodes[n1.uid_].adjacency_.size(); ++i){
      if (nodes[n1.uid_].adjacency_[i].first == n2.uid_){
        nodes[n1.uid_].adjacency_.erase(nodes[n1.uid_].adjacency_.begin()+i);
        break;
      }
    }
    for (unsigned i=0; i < this->nodes[n2.uid_].adjacency_.size(); ++i){
      if (nodes[n2.uid_].adjacency_[i].first == n1.uid_){
        nodes[n2.uid_].adjacency_.erase(nodes[n2.uid_].adjacency_.begin()+i);
        break;
      }
    }
    return 0;
  }

  /** Remove an edge.
   * @param[in] e  Edge to be removed
   * @pre has_edge(n1,n2)==true. With n1=e.node1() and n2=e.node2(),
   *      or n2=e.node1() and n1=e.node2().
   * @return 0
   * @post new num_edges() == old num_edges()-1.
   * The edge is removed and all its information
   * Complexity: No more than O(num_edges())
   */
  size_type remove_edge(const Edge& e){
    Node n1 = e.node1();
    Node n2 = e.node2();
    return remove_edge(n1,n2);
  }

  /** Remove an edge.
   * @param[in] e_it  Edge iterator pointing the edge to be removed
   * @pre has_edge(n1,n2)==true. With n1=*e_it.node1(), n2=*e_it.node2(),
   *      or n2=*e_it.node1() and n1=*e_it.node2().
   * @return edge_iterator e_it. e_it points to the next edge. If no next
   *         edge exists, e_it==edge_end().
   * @post new num_edges() == old num_edges()-1.
   * The edge is removed and all its information. Iterators [e_it, last) are invalidated
   * Complexity: No more than O(num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it){
    Edge e = *e_it;
    size_type removed = remove_edge(e);
    (void) removed;

    /** Check if the iterator e_it points the next edge. If it doesn't because
     it is invalid, search for the iterator to next edge and return it */
    size_type node_id = e_it.graph_->i2u_[e_it.n_idx_];
    if (e_it.index_ == e_it.graph_->nodes[node_id].adjacency_.size()){
      e_it.index_=0;
      ++e_it.n_idx;
      while (e_it.n_idx_< e_it.graph_->num_nodes()){
        size_type node_id = e_it.graph_->i2u_[e_it.n_idx_];
        while (e_it.index_ < e_it.graph_->nodes[node_id].adjacency_.size()){
          if (e_it.graph_->nodes[node_id].adjacency_[e_it.index_].first > node_id){
             return e_it;//return next edge in the adjacency matrix
          }
          ++e_it.index_;
        }
        e_it.index_= 0;
        ++e_it.n_idx_;
      }
    }
    return e_it; //return edge_end()
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    i2u_.clear();
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

    /** Return the value of the iterator.
    * @post result is a node n with n.uid_==index_ and n.graph_==graph_
    */
    Node operator*() const{
      return graph_->node(index_);
    }

    /** Return the next iterator.
    * @post result is a NodeIterator with new_index_ == old_index_ + 1
    */
    NodeIterator& operator++(){
      ++index_;
      return *this;
    }

    /** Test whether this NodeIterator and @a n are equal.
    *
    * Equal NodeIterator have the same graph and the same index.
    */
    bool operator==(const NodeIterator& n) const{
      return index_==n.index_ and graph_==n.graph_;
    }

   private:
    friend class Graph;
     const graph_type* graph_;
     size_type index_;
     NodeIterator(const Graph* graph_pointer, size_type index)
       : graph_(graph_pointer), index_(index){
     }
  };

  /** Return the first node iterator.
  * @post result == node_iterator iter such that *iter.index() == 0
           and *iter.graph_ == this.
  */
  node_iterator node_begin() const{
     return NodeIterator(this,0);
  }

  /** Return the last node iterator.
    * @post result == node_iterator iter such that *iter.index() == n and
             *iter.graph_ == this. Where n is the total number of nodes.
    */
  node_iterator node_end() const{
     return NodeIterator(this, this->num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    /** Return the value of the iterator.
    * @post result is an edge e with e.uid_==index_, e.graph_==graph_,
             e.node1()==node_id_ and e.node2()==n. Such that edge e connects
             node n and node node_id_.
    */
    Edge operator*() const{
      return Edge(graph_,node_id_,graph_->nodes[node_id_].adjacency_[index_].first);
    }

    /** Return the next iterator.
    * @post result is a IncidentIterator with new_index_ == old_index_ + 1
    */
    IncidentIterator& operator++(){
      ++index_;
      return *this;
    }

    /** Test whether this IncidentIterator and @a e are equal.
    *
    * Equal IncidentIterators have the same graph, the same index
       and the same node_id_.
    */
    bool operator==(const IncidentIterator& e) const{
      return index_==e.index_ and graph_==e.graph_ and node_id_==e.node_id_;
    }

   private:
    friend class Graph;
     const graph_type* graph_;
     size_type node_id_;
     size_type index_;
     IncidentIterator(const Graph* graph_pointer,size_type node_id, size_type index)
       : graph_(graph_pointer), node_id_(node_id), index_(index){
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

    /** Return the value of the iterator.
    * @post result is an edge e
    */
    Edge operator*() const{
      return Edge(graph_,graph_->i2u_[n_idx_],graph_->nodes[graph_->i2u_[n_idx_]].adjacency_[index_].first);
    }

    /** Return the next iterator.
    * @post result is a EdgeIterator ei_new such that *ei_new==edge(i+1). Where i satisfies
            *ei_old==edge(i);
    */
    EdgeIterator& operator++(){
      ++index_;
      while (n_idx_< graph_->num_nodes()){
        size_type node_id = graph_->i2u_[n_idx_];
        while (index_ < graph_->nodes[node_id].adjacency_.size()){
          if (graph_->nodes[node_id].adjacency_[index_].first > node_id){
             return *this;
          }
          ++index_;
        }
        index_= 0;
        ++n_idx_;
      }
      return *this;
    }

    /** Test whether this EdgeIterator and @a e are equal.
    *
    * Equal EdgeIterators have the same graph and the same indeces.
    */
    bool operator==(const EdgeIterator& e) const{
      return index_==e.index_ and graph_==e.graph_ and n_idx_==e.n_idx_;
    }

   private:
    friend class Graph;
     const graph_type* graph_;
     size_type index_;
     size_type n_idx_;
     EdgeIterator(const Graph* graph_pointer, size_type index, size_type n_idx)
       : graph_(graph_pointer), index_(index), n_idx_(n_idx){
     }
  };

  /** Return the first edge iterator.
  * @post result == edge_iterator iter pointing the first
          edge in the adjacency matrix.
  */
  edge_iterator edge_begin() const{
    auto ei = EdgeIterator(this,0,0);
    while (ei.n_idx_< ei.graph_->num_nodes()){
      size_type node_id = ei.graph_->i2u_[ei.n_idx_];
      while (ei.index_ < ei.graph_->nodes[node_id].adjacency_.size()){
        if (ei.graph_->nodes[node_id].adjacency_[ei.index_].first > node_id){
           return ei;//return first edge in the adjacency matrix
        }
        ++ei.index_;
      }
      ei.index_= 0;
      ++ei.n_idx_;
    }
    return ei;//return edge_last() if the adjacency matrix is empty
  }
  /** Return the last edge iterator.
  * @post result == edge_iterator iter such that *iter.index_ == num_nodes(),
           *iter.n_idx=i2u_.size(),and *iter.graph_ == this. Where e is the
           total number of edges in the graph.
  */
  edge_iterator edge_end() const{
    return EdgeIterator(this,0,this->num_nodes());
  }
 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  // * @nodes is a list with the information of nodes contained
  //          in node_info objects.
  std::vector<node_info> nodes;
  // * @i2u_ maps indeces to node's uids;
  std::vector<size_type> i2u_;
};

#endif // CME212_GRAPH_HPP
