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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
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
 
    size_type index_;
    Graph* graph_;

    Node() {
      // HW0: YOUR CODE HERE
     this->graph_=NULL;
     this->index_=0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE

	  return this->graph_->master_points[this->index_];
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
       return this->graph_->values[this->index()];
    }
    const node_value_type& value() const{
       return this->graph_->values[this->index()];
    }

    //returns the number of edges that include this node as one of the nodes
    size_type degree() const{
       return this->graph_->connections[this->index()].size();
    }
    // returns an incident iterator pointing at the first edge of this node
    incident_iterator edge_begin() const{
       return incident_iterator(this,0); 
    }
    // returns an incident iterator pointing at one past the last edge of this node
    incident_iterator edge_end() const{
       return incident_iterator(this,this->graph_->connections[this->index()].size());
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
         if (this->index_<n.index()){
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
    Node(Graph* graph, size_type index){
    this->index_=index; 
    this->graph_=graph;
    }

 };  //end of node


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->master_nodes.size();
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
       for (size_type i=0; i<this->master_nodes.size();++i){
           if (this->master_nodes[i]==n){
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
    return this->master_nodes[i];     
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
         if (this->index_<e.index_){
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
    this->index_=node1->graph_->master_edges.size();
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->master_edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return this->master_edges[i];        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
   /*// old way
    Edge new_edge(&a,&b);
    for (size_type i=0; i<this->master_edges.size(); ++i){
        if (this->master_edges[i]==new_edge){
          // return true;
         }
    }
   // return false;
 // }
*/ //new way
    Edge new_edge(&a,&b);
           for (size_type i=0; i<this->connections[a.index()].size(); ++i){
              if((this->connections[a.index()][i])==new_edge){
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
    if (has_edge(a,b)){
       return Edge(&a,&b);
    }
    Edge new_edge(&a,&b);
    Edge new_edge2(&b,&a);
    this->master_edges.push_back(new_edge);
    this->connections[new_edge.node1_->index()].push_back((new_edge));
    this->connections[new_edge.node2_->index()].push_back((new_edge2));
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
      return (this->g_->master_nodes[this->index_]);
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
     return NodeIterator(this,this->master_nodes.size());
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
    Edge operator*() const{
       return (this->node_->graph_->connections[this->node_->index_][this->index_]);
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
    IncidentIterator(const Node* node, size_type index) {
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
     Edge operator*() const{
         return (this->g_->master_edges[this->index_]);
     }
     //increments the edge iterator
     EdgeIterator& operator++(){
         this->index_+=1;
         return  *this;
     }
     // checks if the edge iterator point to same edge
     bool operator==(const EdgeIterator& II) const{
         return((this->index_)==(II.index_));
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
       EdgeIterator end(this,this->master_edges.size());
       return end;
   }
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
};

#endif // CME212_GRAPH_HPP
