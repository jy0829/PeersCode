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
template <typename V, typename E>
class Graph {
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V,E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  using node_value_type = V;
  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes_(), i2u() {}

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
    /** Construct an invalid node.*/
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      //assert(valid());
      return graph_->nodes_[uid_].p_;
    }

    /** Return this node's position. */
    Point& position(){
      //assert(valid());
      return graph_->nodes_[uid_].p_;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //assert(valid());
      return graph_->nodes_[uid_].index_;
    }

    /** Return the uid of the node*/
    size_type uid() const {
      //assert(valid());
      return uid_;
    }

    /** Return this nodes value. */
    node_value_type& value(){
      //assert(valid());
      return graph_->nodes_[uid_].value_;
    }

    /** Return this nodes value. */
    const node_value_type& value() const{
      //assert(valid());
      return graph_->nodes_[uid_].value_;
    }
    
    /** Return the degree of this node*/
    size_type degree() const{
      //assert(valid());
      return graph_->nodes_[uid_].adj_.size();
    }

    /** Return the first Incident Iterator of this node*/
    incident_iterator edge_begin() const{
      return incident_iterator(graph_, index(), 0);
    }

    /**Return this first invalid Incident Iterator of this node*/
    incident_iterator edge_end() const{
      return incident_iterator(graph_, index(), degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      //assert(valid());
      return (uid_== n.uid_ && graph_==n.graph_);
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
      //assert(valid());
      return (graph_< n.graph_ || (graph_==n.graph_ && uid_<n.uid_));
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type uid_;
    Node(const graph_type* graph, size_type uid)
       : graph_(const_cast<graph_type*>(graph)), uid_(uid){
    }

    //Testing for validity of a node
    bool valid() const{
      return (uid_>=0 && uid_ < graph_->nodes_.size() &&  
      graph_->nodes_[uid_].index_ < (int)graph_->i2u.size() &&
      graph_->i2u[graph_->nodes_[uid_].index_] == uid_);
    }

 };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u.size();
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
    node_type new_node=Node(this,nodes_.size());

    i2u.push_back(nodes_.size());
    std::vector<edge_info> adj;
    node_info info(position, value, adj, num_nodes()-1);
    nodes_.push_back(info);


    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (nodes_[n.uid_].index_<(long)num_nodes() && nodes_[n.uid_].index_ != -1 && n.graph_==this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u[i]);
  }


/** Remove a node from the graph
 * @param[in] a   The node to be removed
 * @return        1 if not is removed successfully, 0 otherwise 
 * @post          new num_nodes() = old num_nodes() - 1 if 1 is returned
 *                new num_edges() = old num_edges() - a.degree() if 1 is returned
 *                If the last node is not removed, the node with position num_nodes()-1 will now have index()
 *                of the previously removed node. If the last node is removed then no indices will change.
 *           
 * @notes         All Node objects besides @a a remain valid
 *                All Edge objects that are not incident to @a a remain valid
 *                Invalidates node iterators node_end()-1 and the iterator pointing to @a, 
 *                does not affect other node iterators
 *                Invalidates all incident and edge iterators that point to edges that are incident to @a a or
 *                are incident to a node that is adjacent to @a a . All other edge and incident iterators
 *                remain unaffected.
 *
 * Complexity:  O(1), assuming that the max degree of each node is O(1)(Or O(d^2) if d is the max degree of a node)
 */  
  size_type remove_node(const Node& a){
    if(has_node(a)){

      //Remove each edge that is incident to a
      auto it=a.edge_begin();
      while(it!=a.edge_end()){
        it=remove_edge(it);
      }
   
      //Replace the node to be deleted in i2u by its last element and then remove the last one
      //Update the index of the node that was moved
      int swap=i2u.back();
      i2u[nodes_[a.uid_].index_]=swap;
      nodes_[swap].index_=nodes_[a.uid_].index_;
      i2u.pop_back();
     

      //Set the index of the deleted node to -1
      nodes_[a.uid_].index_=-1;

      return 1;
    }
    return 0;
  }

/** Remove a node from the graph
 * @param[in][out]  n_it An iterator to the node that will  be removed
 * @return          If removing the last element, return the iterator node_end() 
 *                  Otherwise, return an iterator to the object that had previously been at node_end()-1, 
 * @post            new num_nodes() = old num_nodes() - 1 if @a n_it pointed to a valid node in the graph
 *                  new num_edges() = old num_edges() - a.degree() if @a n_it pointed to a valid node in the graph
 *                  If the last node is not removed, the node with position num_nodes()-1 will now have index()
 *                  of the previously removed node. If the last node is removed then no indices will change.
 *
 * @notes           All Node objects besides @a *n_it remain valid
 *                  All Edge objects that are not incident to @a *n_it remain valid
 *                  Invalidates node iterators node_end()-1 and n_it, 
 *                  does not affect other node iterators
 *                  Invalidates all incident and edge iterators that point to edges that are incident to @a a or
 *                  are incident to a node that is adjacent to @a a . All other edge and incident iterators
 *                  remain unaffected.
 *
 * Complexity:  O(1), assuming that the max degree of each node is O(1)(Or O(d^2) if d is the max degree of a node)
 */ 
  node_iterator remove_node(node_iterator n_it){
    remove_node(*n_it);
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, uid2_);
    }

    /**Returns the lenght of this Edge*/
    double length() const{
      return norm(node1().position()-node2().position());
    }

    /**Returns the value associated to this Edge*/
    const edge_value_type& value() const{
      return graph_->nodes_[uid1_].adj_[adj_index(uid1_,uid2_)].value_;
    }
    /**Returns the value associated to this Edge*/
    edge_value_type& value(){
      return graph_->nodes_[uid1_].adj_[adj_index(uid1_,uid2_)].value_;
     }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_==e.graph_ && ((uid1_ ==e.uid1_ && uid2_ == e.uid2_)||
             (uid1_ == e.uid2_ && uid2_ == e.uid1_)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_<e.graph_){
        return true;}
      else if(graph_==e.graph_ && std::min(uid1_,uid2_)<std::min(e.uid1_,e.uid2_)){
        return true;}
      else if(graph_==e.graph_ && std::min(uid1_,uid2_)==std::min(e.uid1_,e.uid2_)
             && std::max(uid1_,uid2_)<std::max(e.uid1_,e.uid2_)){
        return true;}
      return false;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type uid1_;
    size_type uid2_;
    Edge (const Graph* graph, size_type uid1, size_type uid2) :
      graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {}

    //Helper method to find the index where node with uid b is stored
    //in the adjacency vector associated to the node with uid a
    size_type adj_index(size_type a, size_type b) const{
      for (size_type i=0; i < graph_->nodes_[a].adj_.size(); ++i) {
        if (graph_->nodes_[a].adj_[i].uid_other_ == b) { 
          return i;
        }
      }
      return graph_->nodes_[a].adj_.size();
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return std::distance(edge_begin(),edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(),i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto i=a.edge_begin(); i !=a.edge_end(); ++i) {
      if ((*i).node2().uid_ == b.uid_) { 
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {

    if(has_edge(a,b)){
      return Edge(this, a.uid_, b.uid_);
    }

    edge_info A=edge_info(a.uid_, value);
    edge_info B=edge_info(b.uid_, value);

    nodes_[a.uid_].adj_.push_back(B);
    nodes_[b.uid_].adj_.push_back(A);

    return Edge(this, a.uid_,b.uid_);
  }
  
  /** Remove an edge from the graph
 * @param[in] a,b The nodes that are incident to the edge that is to be removed
 * @return        1 if not is removed successfully, 0 otherwise 
 * @post          new num_edges() = old num_edges() - 1 if 1 is returned
 *                new a.degree()  = old a.degree() -1 if 1 is returned
 *                new b.degree()  = old b.degree() -1 if 1 is returned
 * 
 * @notes         All Node objects remain valid.
 *                All Edge objects besides the one being deleted remain valid.
 *                All node iterators are unaffected.
 *                Invalidates all incident and edge iterators that are incident to @a a or @a b.
 *                All other incident and edge iterators remain valid.
 *
 * Complexity:  O(1), assuming that the max degree of each node is O(1)(Or O(d) if d is the max degree of a node)
 */  
  size_type remove_edge(const Node& a, const Node& b){
    return remove_edge(Edge(this, a.uid_, b.uid_));
  }

/** Remove an edge from the graph
 * @param[in] e   The edge that is to be removed
 * @return        1 if not is removed successfully, 0 otherwise 
 * @post          new num_edges() = old num_edges() - 1 if 1 is returned
 *                new node1().degree() = old node1().degree() -1 if 1 is returned
 *                new node2().degree() = old node2().degree() -1 if 1 is returned 
 *
 *           
 * @notes         All Node objects remain valid.
 *                All Edge objects besides the one being deleted remain valid.
 *                All node iterators are unaffected.
 *                Invalidates all incident and edge iterators that are incident to either nodes of @a e.
 *                All other incident and edge iterators remain valid.
 *
 * Complexity:  O(1), assuming that the max degree of each node is O(1)(Or O(d) if d is the max degree of a node)
 */
  size_type remove_edge(const Edge& e){
    if(has_edge(e.node1(), e.node2())){
      
      //Finding where the edge is stored
      size_type adj1=e.adj_index(e.uid1_,e.uid2_);
      size_type adj2=e.adj_index(e.uid2_,e.uid1_);

      //Replacing the edge to be removed with the last edge and then removing the last element
      nodes_[e.uid1_].adj_[adj1]=nodes_[e.uid1_].adj_.back();
      nodes_[e.uid1_].adj_.pop_back();
      nodes_[e.uid2_].adj_[adj2]=nodes_[e.uid2_].adj_.back();
      nodes_[e.uid2_].adj_.pop_back();
      return 1;
    }
    return 0;
  }

/** Remove an edge from the graph
 * @param[in] e_it An edge iterator to the edge that is to be removed
 * @return         If e_it is edge_end() -1, return the end iterator. Else if e_it points to the same edge as 
 *                 node1().edge_end()-1 or the edge iterator pointing to node1().edge_end() is not valid,
 *                 return the next valid edge iterator @a e_it. Otherwise return the edge 
 *                 iterator that was previously pointing to the same edge as node1().edge_end().
 *
 * @post           new num_edges() = old num_edges() - 1 if 1 is returned
 *                 new node1().degree() = old node1().degree() -1 if edge is succesfully removed
 *                 new node2().degree() = old node2().degree() -1 if edge is succesfully removed 
 *
 * @notes          All Node objects remain valid.
 *                 All Edge objects besides the one being deleted remain valid.
 *                 All node iterators are unaffected.
 *                 Invalidates all incident and edge iterators that are incident to either of the nodes of @a *e_it
 *                 All other incident and edge iterators remain valid.
 *
 * Complexity:  O(1), assuming that the max degree of each node is O(1)(Or O(d) if d is the max degree of a node)
 */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    e_it.fix();
    return e_it;
  }

/** Remove an edge from the graph
 * @param[in] e_it An incident iterator to the edge that is to be removed
 * @return         If e_it points to the last element of the adjacency vector of node1(), return
 *                 the end iterator for node1(). Otherwise, return an incident iterator that points
 *                 to the previously last element of the adjacency vector of node1().
 *
 * @post           new num_edges() = old num_edges() - 1 if 1 is returned
 *                 new node1().degree() = old node1().degree() -1 if edge is succesfully removed
 *                 new node2().degree() = old node2().degree() -1 if edge is succesfully removed 
 *
 * @notes          All Node objects remain valid.
 *                 All Edge objects besides the one being deleted remain valid.
 *                 All node iterators are unaffected.
 *                 Invalidates all incident and edge iterators that are incident to either of the nodes of @a *e_it
 *                 All other incident and edge iterators remain valid.
 *
 * Complexity:  O(1), assuming that the max degree of each node is O(1)(Or O(d) if d is the max degree of a node)
 */

  incident_iterator remove_edge(incident_iterator e_it){
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    i2u.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator> {
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

    /** Dereferences the Node Iterator and returns its corresponding Node*/
    Node operator*() const{
      return Node(graph_, graph_->i2u[index_]); 
    }

    /**Increments the Node Iterator and returns the next Node Iterator*/
    node_iterator& operator++(){
      index_++;
      return *this;
    }

    /**Tests whether the current Node Iterator and @a x are equal*/
    bool operator== (const NodeIterator& x) const{
      return (index_ == x.index_ && graph_ == x.graph_) ;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type index_;
    NodeIterator (const graph_type* graph, size_type index) :
                  graph_(const_cast<graph_type*>(graph)), index_(index) {}
  };

  /** Returns a Node Iterator that corresponds to the first Node*/
  node_iterator node_begin() const{
    return NodeIterator(this,0); 
  }
  
  /** Returns a Node Iterator that corresponds to the first invalid Node*/
  node_iterator node_end() const{
    return NodeIterator(this, num_nodes());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private equality_comparable<IncidentIterator> {
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

    /** Deferences the Incident Iterator and returns the corresponding Edge*/
    Edge operator*() const{
      return Edge(graph_,graph_->i2u[index_],graph_->nodes_[graph_->i2u[index_]].adj_[aux_index_].uid_other_);
    }

    /**Increments the Incident Iterator and returns the next Incident Iterator*/
    IncidentIterator& operator++(){
      ++aux_index_;
      return *this;
    }
    /** Test whether the current Incident Iterator and @a x are equal*/
    bool operator == (const IncidentIterator& x) const{
      return (graph_ == x.graph_ && index_ == x.index_ && aux_index_ == x.aux_index_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type index_;
    //This is the index of where the second node is stored 
    //in the adjacency matrix of @a index_
    size_type aux_index_;
    IncidentIterator (const graph_type* graph, size_type index, size_type aux_index) :
                     graph_(const_cast<graph_type*>(graph)), index_(index), aux_index_(aux_index){}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private equality_comparable<EdgeIterator>{
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

    /** Deferences the Edge Iterator and returns the corresponding Edge*/
    Edge operator*() const{
      return Edge(graph_,graph_->i2u[index_], graph_->nodes_[graph_->i2u[index_]].adj_[aux_index_].uid_other_);
    }
   /**Increments the Edge Iterator and calls fix() which makes sure the 
     *increment was legitimate and returns the next valid Edge Iterator*/
    EdgeIterator& operator++(){
      ++aux_index_;
      fix();
      return *this;
    }

    /** Test whether the current Edge Iterator and @a x are equal*/
    bool operator==(const EdgeIterator& x) const {
      return (graph_ == graph_ && index_ == x.index_ && aux_index_ == x.aux_index_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type index_;
    //This is the index of where the second node is stored
    //in the adjacency matrix of @a index_
    size_type aux_index_;
    EdgeIterator (const graph_type* graph, size_type index, size_type aux_index) :
                     graph_(const_cast<graph_type*>(graph)), index_(index), aux_index_(aux_index){
      fix();
    }
    // Function that increments the edge iterator until a valid Edge Iterator is reached 
    void fix(){
      while (true){
        //There's no more valid nodes 
        if (index_ >= graph_->num_nodes()){
          break;
        }
        //We've gone too far on this node
        else if(aux_index_ >= graph_->nodes_[graph_->i2u[index_]].adj_.size()){
          index_++;
          aux_index_=0;
        }
        //Only consider when the first node is smaller so we don't double count
        else if(graph_->i2u[index_]>graph_->nodes_[graph_->i2u[index_]].adj_[aux_index_].uid_other_){
          aux_index_++;
        }
        //We've got a good Edge Iterator
        else{
          break;
        }
      }
    }
  };

  /** Returns the first Edge Iterator of the graph*/
  edge_iterator edge_begin() const{
    return edge_iterator(this,0,0);
  }
  /** Returns the first invalid Edge Iterator of the graph*/
  edge_iterator edge_end() const{
    return edge_iterator(this, num_nodes(), 0);
  }

  /** Returns the index of the node with the inputted uid*/
  size_type u2i(size_type uid) const{
    return nodes_[uid].index_;
  }


 private:
  struct edge_info{
    size_type uid_other_;
    edge_value_type value_;
    edge_info(size_type uid_other, edge_value_type value) : uid_other_(uid_other), value_(value) {}
  };

  struct node_info{
    Point p_;
    node_value_type value_;
    std::vector<edge_info> adj_;
    long index_;
    node_info(Point p, node_value_type value, std::vector<edge_info> adj, size_type index) :
             p_(p), value_(value), adj_(adj), index_(index) {}
  };

  std::vector<node_info> nodes_;
  std::vector<size_type> i2u;
};



#endif // CME212_GRAPH_HPP
