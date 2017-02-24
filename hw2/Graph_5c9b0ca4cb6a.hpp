
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
 private:
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct node_info;
  struct edge_info;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_value_type = V;
  using edge_value_type =E;

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

  using idx_type = size_type;
  using uid_type = size_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes(), i2u(){
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
    }
/* Return if the node is valid or not
 * For the node to be valid it must have a positive uid 
 * and the index and uid must coincide between the nodes and the i2u vectors
 */

    bool valid() const{
      return uid_ >= 0 && uid_ < graph_ -> nodes.size() 
           && graph_ -> nodes[uid_].index < graph_->i2u.size()
           && graph_ -> i2u[graph_ -> nodes[uid_].index] == uid_;

    }

   

    /** Return this node's position. */
    const Point& position() const {
//      assert(valid());
      return graph_ -> nodes[uid_].position;
    }

    Point& position() {
//      assert(valid());
      return graph_ -> nodes[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_ -> nodes[uid_].index;
    }

    size_type& index() {
      return graph_ -> nodes[uid_].index;
    }
    /** Return the value associated to the node. */
    node_value_type& value(){
      return graph_ -> nodes[uid_].val;
    }
    const node_value_type& value() const{
      return graph_ -> nodes[uid_].val;
    }
    /** The degree of the node is the size of the vector of adjacent nodes. */
    size_type degree() const{
      return graph_ -> nodes[uid_].adj.size();
    }
    /** The incident iterator begins at the first index of the vector of adjacent nodes.
     * @post @a id = 0 
     */ 

    IncidentIterator edge_begin() const{
      return IncidentIterator(graph_, uid_, 0);
    }

     /** The incident iterator ends when there is no more 
      *  nodes adjacent to it
      *  @post @a id = degree()
      */ 
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_, degree());
    }    

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph_ == graph_ and uid_ == n.uid_);           
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
      if(n.graph_ == graph_){
        if(uid_ < n.uid_)
          return true;
      }
      else{
        if(graph_ < n.graph_)
          return true;
      }
      return false;
    }



   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_; //Pointer back to the Graph container
    uid_type uid_;// This element's unique identification number 
   

    Node(const Graph* graph, size_type node_uid)
        : graph_(const_cast<Graph*>(graph)), uid_(node_uid){
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
    return i2u.size(); // The valid nodes are in i2u
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
    node_info n_new = node_info(position,val);
    n_new.index = size(); // The new index available
    nodes.push_back(node_info(position, val));
    i2u.push_back(nodes.size()-1);
    return Node(this, nodes.size() -1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      return (n.graph_ == this) and (n.index() < num_nodes()) and (n.index() > -1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i2u[i]); //Return the node from the graph that has index i (see constructor) 
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_uid);  
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_uid);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if( graph_ == e.graph_){
        return (node1_uid == e.node1_uid and node2_uid == e.node2_uid) or (node2_uid == e.node1_uid and node1_uid == e.node2_uid);

      }
      return false;
    }
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(e.graph_ == graph_){
        if(node1_uid == e.node1_uid){
          return node2_uid < e.node2_uid;
        }
        else {
         return node1_uid < e.node1_uid;
        }
      }
      else
        return graph_ < e.graph_;
    }

    double length() const {
      return norm(node1().position() - node2().position());
    }
 
    edge_value_type& value(){
      for(size_type j = 0; j < graph_ -> nodes[node1_uid].adj.size(); j++){
        if(graph_ -> nodes[node1_uid].adj[j].uid_other == node2_uid)
          return graph_ -> nodes[node1_uid].adj[j].val; 
      }
    }

    const edge_value_type& value() const {
      for(size_type j = 0; j < graph_ -> nodes[node1_uid].adj.size(); j++){
        if(graph_ -> nodes[node1_uid].adj[j].uid_other == node2_uid)
          return graph_ -> nodes[node1_uid].adj[j].val; 
      }
    }



   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type node1_uid; //Uid of the "first" (no distinction) node
    size_type node2_uid; //Uid node2

    Edge(const Graph* graph, size_type node1, size_type node2)
      : graph_(const_cast<Graph*>(graph)), node1_uid(node1), node2_uid(node2){
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    int cmpt = 0;
    for(unsigned int i = 0; i < nodes.size();++i){
      cmpt += nodes[i].adj.size();
    }
    return cmpt/2; // The edges are counted twice
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
 //   if(a.valid() and b.valid()){
    for(size_type j = 0; j < nodes[a.uid_].adj.size(); j++){
      if(nodes[a.uid_].adj[j].uid_other == b.uid_)
        return true; 
    }
    for(size_type j = 0; j < nodes[b.uid_].adj.size(); j++){
      if(nodes[b.uid_].adj[j].uid_other == a.uid_)
        return true; 
    }
    
   // }
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
    if(!(has_edge(a,b))){
      nodes[a.uid_].adj.push_back(edge_info(b.uid_));
      nodes[b.uid_].adj.push_back(edge_info(a.uid_));
    }
    return Edge(this, a.uid_, b.uid_);
    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    i2u.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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
    /** Return the node corresponding to the node_id_
     *  @pre @a node_id_ < num_nodes()
     */
    Node operator*() const{
      return graph_ ->node(node_id_); 
    }
    /** We implement the operator ++ by incrementing the node id
     * @pre node_id_ < num_nodes()
     * @post node_id_ <= num_nodes()
     */
    NodeIterator& operator++() {
      ++node_id_;
      return *this  ;
    }
    /** Check the equalitu of the graph and of the id of the node
    */
    bool operator ==(const NodeIterator& x) const{
      return (graph_ == x.graph_ and node_id_ == x.node_id_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_id_;//Index of the node 
    
    NodeIterator(const Graph* graph, size_type node_id) 
      : graph_(const_cast<Graph*>(graph)), node_id_(node_id){
       assert(node_id_ <= graph_ -> size());
    }    
  };
  /** Return the node iterator starting at the first node
   */
  NodeIterator node_begin() const{

    return NodeIterator(this,0);     
   }
 /** Return the node iterator that can not be incremented or dereferenced
  * That is the last node of the graph 
  */

  NodeIterator node_end() const{
     return NodeIterator(this, this->size());
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
    IncidentIterator() {
    }
    /** Return the edge that links the nodes with index 
     * @a node_1 to the one with index @a node2
     * @pre this edge exists
     */
    Edge operator*() const{
     size_type node2_uid = graph_->nodes[main_node_uid_].adj[idx_adj_].uid_other;
      return Edge(graph_, main_node_uid_, node2_uid);
    }
    /** Return the incremented incidentiterator 
     * @pre the incident iterator hasn't reached the end yet
     * @post @a new idx_adj_ = @a old idx_adj_
     */
    IncidentIterator& operator++(){
      ++idx_adj_;
      return *this;
    }
    /** Test the equality of two incident iterator 
     * @param the incident iterator that we want to compare with 
     * @result true if all components are the same 
     */
    bool operator==(const IncidentIterator& x) const{
      return (x.graph_ == graph_ and x.main_node_uid_ == main_node_uid_ and x.idx_adj_ == idx_adj_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type main_node_uid_; // This will not change: we will iterate to its edges
    size_type idx_adj_;
    IncidentIterator(const Graph* graph, size_type node1, size_type node2) 
      : graph_(const_cast<Graph*>(graph)), main_node_uid_(node1), idx_adj_(node2){
      assert(idx_adj_ <= graph_-> nodes[main_node_uid_].adj.size());
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
    Edge operator*() const{
      size_type node2_ = graph_->nodes[graph_ -> i2u[node1_]].adj[idx_].uid_other;
      return Edge(graph_, graph_ ->i2u[node1_], node2_);
    }
    /** The operator to increment our edge iterator
     *  We return only edges such as the id of node1 is < id of node2
     *  We start by @a node1_ then we go through the idx of the adjency vector 
     *  For each idx if the node id related to it, is < @a node1_ then it is the incremented edge iterator
     *  @pre node1_ < graph_ -> adjacency[node1_][idx_]
     *  @post node2_ < graph_ -> adjacency[node1_][idx_]
     *  @pre and @post  idx_ < graph_->adjacency[node1_].size()-1
     */
    EdgeIterator& operator++(){
      while(node1_ < graph_->size()){

        while(idx_ < graph_->nodes[graph_->i2u[node1_]].adj.size()-1){
          ++idx_;
          if(graph_ -> i2u[node1_] < graph_->nodes[graph_->i2u[node1_]].adj[idx_].uid_other){
            return *this;
          }
        }
        ++node1_;
        idx_ = 0;
      }
      return *this;
    }
 
    bool operator==(const EdgeIterator& x) const{
      return (graph_ == x.graph_ and node1_ == x.node1_ and idx_ == x.idx_);

    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node1_; //index of node
    size_type idx_; // position of the other node in the vector adj of the previous node
    EdgeIterator(const Graph* graph, size_type n1, size_type id) 
      : graph_(const_cast<Graph*>(graph)), node1_(n1), idx_(id){
        fix();
    } 

    void fix(){
      while(true){
        if(node1_ == graph_-> num_nodes()){
          break;
        }
        else if(idx_ == graph_->nodes[graph_->i2u[node1_]].adj.size()){
          idx_ = 0;
          ++node1_;  

        }
        else if(graph_->i2u[node1_] > graph_->nodes[graph_->i2u[node1_]].adj[idx_].uid_other){
          ++idx_;
        }
        else{ break;}
      }

     


    }
  };
 /** First iterator, given our ++ operator, we start at @a node1_ = 0 and @a idx_ = 0
  * We are sure that 0 < adjacency[node1_][idx_]
 */
  edge_iterator edge_begin() const{
    return edge_iterator(this, 0, 0);
  }
 /** The last iterator is when node1_ can not be incremented */
  edge_iterator edge_end() const{
    return edge_iterator(this, num_nodes(),0);
  }

/** A fonction to remove a node from the graph 
 * @param @a n node to be removed 
 * @result 1 if we effectively removed the node, 0 if there was no node to remove
 * @post all edges previously adjacent to this node are removed
 * @post the index of @a n does not correspond to the uid of @a n in i2u 
 * @post all other nodes are not affected and stay valid
 * We need to re-define the index of the nodes to have the concordance n
 * necessary for nodes to be valid

 */

  size_type remove_node(const Node& n){
    if(!n.valid()){return 0;}
    if(num_nodes() ==0){return 0;}
    //We need to remove all edges that are incident to n
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){
      remove_edge(*ei);
    }
    size_type idx = n.index();
    //First, we remove this index from i2u
    i2u.erase(i2u.begin() + idx); // erase the element 
    //Then we need to re-do the concordance with the index in nodes: 
    for(unsigned int i = idx; i < i2u.size(); ++i){
      nodes[i2u[i]].index = i;
    }
   
    return 1;
  }
/*Remove the node that @a n_it points to
 *@param n_it the node iterator pointing to the node that needs to me removed
 *@result a node_iterator that does not iterate through @a old *n_it
 */
  node_iterator remove_node(node_iterator n_it){
    Node n = *n_it;
    remove_node(n);
    return n_it;

  }

/** A fonction to remove an edge from the graph 
 * @param @a a and @a b the two nodes of the edge to be removed 
 * @result 1 if we effectively removed the node, 0 if there was no node to remo$
 * @post the edge from @a a to @a b and the edge from @a b to @a a are removed
 * @post the corresponding index of these nodes are erased from the adjacent vectors
 * @post all other edges and so index of the adjacent vector are not affected

 * I think this function is the problem but I can't see what is going wrong...

 */


  size_type remove_edge(const Node& a, const Node& b){
//    if(!a.valid() or !b.valid()){return 0;}
    if(has_edge(a,b)){
    size_type a_in_b = 0;
    size_type b_in_a = 0;
    for(size_type j = 0; j < nodes[a.uid_].adj.size(); j++){
      if(nodes[a.uid_].adj[j].uid_other == b.uid_){
        b_in_a = j;
       // break;
      }
    }
    for(size_type j = 0; j < nodes[b.uid_].adj.size(); ++j){
      if(nodes[b.uid_].adj[j].uid_other == a.uid_){
        a_in_b = j;
      //  break;
      }
    }
    nodes[a.uid_].adj.erase(nodes[a.uid_].adj.begin() + b_in_a); 
    nodes[b.uid_].adj.erase(nodes[b.uid_].adj.begin() + a_in_b); 
    return 1;
    }
    else{ return 0;}
  } 

/*Remove the edge @a e 
 *@param the edge to be removed 
 *@result 1 if success 
 *@post the edge between the nodes linked by edge @a e is removed 
 *@post the other edges are not affected

 */

  size_type remove_edge(const Edge& e){
    size_type a_in_b = 0;
    size_type b_in_a = 0;
    for(size_type j = 0; j < nodes[e.node1_uid].adj.size(); j++){
      if(nodes[e.node1_uid].adj[j].uid_other == e.node2_uid){
        b_in_a = j;
        break;
      }
    }
    for(size_type j = 0; j < nodes[e.node2_uid].adj.size(); ++j){
      if(nodes[e.node2_uid].adj[j].uid_other == e.node1_uid){
        a_in_b = j;
        break;
      }
    }
    nodes[e.node1_uid].adj.erase(nodes[e.node1_uid].adj.begin() + b_in_a); 
    nodes[e.node2_uid].adj.erase(nodes[e.node2_uid].adj.begin() + a_in_b); 
    return 1;
  }


/*Remove the edge that @a e_it points to
 *@param e_it the edge iterator pointing to the node that needs to be removed
 *@result a edge_iterator that does not iterate through @a old *e_it
 */

  edge_iterator remove_edge(edge_iterator e_it){
    if(e_it == edge_end()){return edge_end();}
    remove_edge(*e_it);
    if(num_edges() == 0){ return edge_end();}
    return e_it;

  }

 private:

 
  /** New structure to take into account the value of the node */
  struct node_info{
    Point position; 
    node_value_type val;
    std::vector<edge_info> adj; //Vector of adjacent edges
    idx_type index; //Index of the node
    node_info(Point p, node_value_type val_) 
      : position(p), val(val_), adj(std::vector<edge_info>()), index(){
    }

  };

  struct edge_info{
    uid_type uid_other; 
    E val;
    edge_info(size_type uid) : uid_other(uid), val(){
    }
    edge_info(size_type uid, E val_) : uid_other(uid), val(val_){
    }
  };



  std::vector<node_info> nodes;
  std::vector<uid_type> i2u;

};

#endif // CME212_GRAPH_HPP
