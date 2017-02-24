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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
  using edge_value_type = E;

  // typedef V node_value_type;
  // typedef E edge_value_type;
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  std::vector<std::vector<Edge>> connections;
  std::vector<V> values;
  std::vector<E> values_e;
  std::vector<size_type> e_i2u_; 
  std::vector<size_type> n_i2u_; 
private: 
  std::vector<Point> master_points;
  std::vector<Node> master_nodes;
  std::vector<Edge> master_edges;
public:
  /** Construct an empty graph. */  
  Graph() {
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
 
    size_type uid_;
    size_type index_;
    Graph* graph_;

    Node() {
      // HW0: YOUR CODE HERE
     this->graph_=NULL;
     this->uid_=0;
     this->index_=0;
    }

    /** Return this node's position. */
    Point& position() {
	  return this->graph_->master_points[this->uid_];
    }
    const Point& position() const {
      // HW0: YOUR CODE HERE

	  return this->graph_->master_points[this->uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // returns value of node
    node_value_type& value(){
       return this->graph_->values[this->uid_];
    }
    const node_value_type& value() const{
       return this->graph_->values[this->uid_];
    }

    //returns the number of edges that include this node as one of the nodes
    size_type degree() const{
       return this->graph_->connections[this->uid_].size();
    }
    // returns an incident iterator pointing at the first edge of this node
    incident_iterator edge_begin() const{
       return incident_iterator(this,0); 
    }
    // returns an incident iterator pointing at one past the last edge of this node
    incident_iterator edge_end() const{
       return incident_iterator(this,this->graph_->connections[this->uid_].size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE 
      if (this->graph_==n.graph_){ 
         if (this->index_ == n.index_){
            return true;
         }
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
      // HW0: YOUR CODE HERE
      if (!(n.graph_<this->graph_)){
         if (this->uid_<n.uid_){
            return true;
         }
      }
      return false;
    }  
     
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(Graph* graph, size_type uid){
    this->uid_=uid;
    this->index_=graph->n_i2u_.size(); 
    this->graph_=graph;
    }

 };  //end of node


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->n_i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }


  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value (optional)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    this->master_points.push_back(position);
    Node new_node(this,this->master_nodes.size());
    this->n_i2u_.push_back(this->master_nodes.size());
    this->master_nodes.push_back(new_node);
    std::vector<Edge> empty;
    this->connections.push_back(empty);
    this->values.push_back(value);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_==this){
       for (size_type i=0; i<this->num_nodes();++i){
           if (this->master_nodes[n_i2u_[i]]==n){
              return true;
           }
       }
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre4 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return this->master_nodes[n_i2u_[i]];     
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
  class Edge :private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }
    const Node* node1_;
    const Node* node2_;
    size_type uid_;
    size_type index_;
    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return *node1_; 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return *node2_;
    }

    double length() const{
       return norm_2(node1_->position()-node2_->position());
    }
    edge_value_type& value(){
       return this->node1_->graph_->values_e[this->uid_];
    }
    
    const edge_value_type& value() const{
      return this->node1_->graph_->values_e[this->uid_];
    }
     
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->node1().graph_==e.node1().graph_){
         if ((*(this->node1_)==*(e.node1_))||(*(this->node1_)==*(e.node2_))){
            if ((*(this->node2_)==*(e.node1_))||(*(this->node2_)==*(e.node2_))){
               return true;
            } 
         }
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->node1().graph_<e.node1().graph_){
         if (this->uid_<e.uid_){
            return true;
         }
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge object
    Edge(const Node* node1,const Node* node2){
       this->node1_=node1;
       this->node2_=node2;
       this->uid_=node1->graph_->master_edges.size();
       this->index_=node1->graph_->e_i2u_.size();
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->e_i2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return this->master_edges[e_i2u_[i]];        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    Edge new_edge(&a,&b);
           for (size_type i=0;i< this->e_i2u_.size(); ++i){
              if(this->edge(i)==new_edge){
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
    // HW0: YOUR CODE HERE
    if (has_edge(a,b)){
       return Edge(&a,&b);
    }
    Edge new_edge(&(this->master_nodes[a.uid_]),&(this->master_nodes[b.uid_]));
    Edge new_edge2(&(this->master_nodes[b.uid_]),&(this->master_nodes[a.uid_]));
    this->e_i2u_.push_back(this->master_edges.size());
    this->master_edges.push_back(new_edge);
    this->connections[new_edge.node1_->uid_].push_back((new_edge));
    this->connections[new_edge.node2_->uid_].push_back((new_edge2));
    this->values_e.push_back(value);
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
    
  void clear() {
    // HW0: YOUR CODE HERE
    this->master_nodes.clear();
    this->master_points.clear();
    this->master_edges.clear();
    this->connections.clear();
    this->n_i2u_.clear();
    this->e_i2u_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator>{
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
    // return node that iterator points to
    Node& operator*() const{
      return (this->g_->master_nodes[this->g_->n_i2u_[this->index_]]);
    }
    // increment interator
    NodeIterator& operator++(){
      this->index_+=1;
      return *this;
    }
    // checks if two node iterators are equal
    bool operator==(const NodeIterator& node) const{
      return this->index_==node.index_;
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type index_;
    Graph* g_;  
    NodeIterator(const Graph* g,const size_type index) {
       this->index_=index;
       this->g_=const_cast<Graph*>(g);
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // returns node iterator starting at the first node
  node_iterator node_begin() const{
     NodeIterator begin(this,0);
     return begin;
  }
  //returns node iterator pointing 1 past the last node
  node_iterator node_end() const{
     return NodeIterator(this,this->size());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator :private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    
    IncidentIterator() {
      this->node_=NULL;
      this->index_=0;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // returns the edge that the incident iterator points to
    Edge& operator*() const{
       return (this->node_->graph_->connections[this->node_->uid_][this->index_]);
    }
    //increments the incident iterator
    IncidentIterator& operator++(){
       this->index_=this->index_+1;
       return *this;
    }
    //checks if two incident iterators are equal
    bool operator==(const IncidentIterator& II) const{
       if (*(this->node_)==*(II.node_)){
          if (this->index_==II.index_){
             return true;
          }
       }
       return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Node* node_;
    size_type index_;
    IncidentIterator(const Node* node, const size_type index) {
       this->node_=const_cast<Node*>(node);
       this->index_=index;
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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
     // returns the edge the edge iterator points to
     Edge& operator*() const{
         return (this->g_->master_edges[this->g_->e_i2u_[this->index_]]);
     }
     //increments the edge iterator
     EdgeIterator& operator++(){
         this->index_+=1;
         return  *this;
     }
     // checks if the edge iterator point to same edge
     bool operator==(const EdgeIterator& II) const{
         return((*this).index_==(II.index_));
     }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* g_;
    size_type index_;
    EdgeIterator(const Graph* g,const size_type index){
       this->index_=index;
       this->g_=const_cast<Graph*>(g);
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  //returns edge iterator starting at the first edge
  edge_iterator edge_begin() const{
       EdgeIterator begin(this,0);
       return begin;
   }
  //returns edge iterator pointing one past the last edge
   edge_iterator edge_end() const{
       EdgeIterator end(this,this->e_i2u_.size());
       return end;
   }
size_type remove_edge(const Edge& e){
                   std::cout<<e.index_<<std::endl;
      if(this->edge_begin()!=this->edge_end()){
         for (auto e_update=EdgeIterator(this,e.index_);e_update!=this->edge_end();++e_update){
//                   std::cout<<"here"<<(*e_update).index_<<std::endl;
               (*e_update).index_=((*e_update).index_-1);
                   
         }
         
         if(e.index_!=this->num_edges()){
            this->e_i2u_.erase(e_i2u_.begin()+e.index_);
         }else{
            this->e_i2u_.pop_back();
         }
         return 1;
      }
      return 0;
   }

   size_type remove_edge(const Node& a, const Node& b){
      size_type i=0;
      for (auto here=a.edge_begin();here!=a.edge_end();++here){
         if ((*here).node2()==b){
             size_type i2=0;
             for (auto here2=b.edge_begin();here2!=b.edge_end();++here2){
                if ((*here2).node2()==a){
                   this->connections[b.uid_].erase(this->connections[b.uid_].begin()+i2);
                   break;
                }
                ++i2;
             }
             remove_edge(*here);
             this->connections[a.uid_].erase(this->connections[a.uid_].begin()+i);
             return 1;
         }
         ++i;
      }
      return 0;
   }

  // edge_iterator remove_edge(edge_iterator e_it){

 //}


   size_type remove_node(const Node& n){
      for (auto n_update=NodeIterator(this,n.index()+1);n_update!=this->node_end();++n_update){
         (*n_update).index_=((*n_update).index_-1);
      }
      this->n_i2u_.erase(this->n_i2u_.begin()+n.index_);
      for (auto ii_it=n.edge_begin();ii_it!=n.edge_end();++ii_it){
          remove_edge(*ii_it);  
      } 
      return 1;      
   }
   
   node_iterator remove_node(node_iterator n_it){ 
      auto index_=(*n_it).index_;    
      for (auto n_update=n_it;n_update!=this->node_end();++n_update){
         (*n_update).index_=((*n_update).index_-1);
      }
      //(*n_it).index_=-1;
      this->n_i2u_.erase(this->n_i2u_.begin()+index_);
      return NodeIterator(this,index_);
   }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
