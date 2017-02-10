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
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

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
  Graph() : nodes_(),adjacency_matrix(std::vector<std::vector<size_type>>(0)){}

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
      return graph_->nodes_[index()].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /** Return this nodes value. */
    node_value_type& value(){
      return graph_->nodes_[index()].second; 
    }

    /** Return this nodes value. */
    const node_value_type& value() const{
      return graph_->nodes_[index()].second;
    }
    
    /** Return the degree of this node*/
    size_type degree() const{
      return graph_->adjacency_matrix[index()].size();
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
      return (index()== n.index_ && graph_==n.graph_);
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
      return (graph_< n.graph_ || (graph_==n.graph_ && index()<n.index_));
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type index_;
    Node(const graph_type* graph, size_type index)
       : graph_(const_cast<graph_type*>(graph)), index_(index){
    }

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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    node_type new_node=Node(this, size());
    nodes_.push_back(std::make_pair(position, value));
    std::vector<size_type> auxadj;
    adjacency_matrix.push_back(auxadj);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index_<size() && n.graph_==this);
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
      return Node(graph_, index_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, index_2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_==e.graph_ && ((index_1 ==e.index_1 && index_2 == e.index_2)||
             (index_1 == e.index_2 && index_2 == e.index_1)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_<e.graph_){
        return true;}
      else if(graph_==e.graph_ && std::min(index_1,index_2)<std::min(e.index_1,e.index_2)){
        return true;}
      else if(graph_==e.graph_ && std::min(index_1,index_2)==std::min(e.index_1,e.index_2)
             && std::max(index_1,index_2)<std::max(e.index_1,e.index_2)){
        return true;}
      return false;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type index_1;
    size_type index_2;
    Edge (const Graph* graph, size_type index_1, size_type index_2) :
      graph_(const_cast<Graph*>(graph)), index_1(index_1), index_2(index_2)
      {}

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
      if ( (*i).node2().index_ == b.index_) { 
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

    if(has_edge(a,b)) {
      return Edge(this, a.index_, b.index_);
    }
    adjacency_matrix[a.index_].push_back(b.index_);
    adjacency_matrix[b.index_].push_back(a.index_);
    return Edge(this, a.index_, b.index_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adjacency_matrix.clear();
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
      return Node(graph_, index_); 
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
      return Edge(graph_,index_,graph_->adjacency_matrix[index_][aux_index_]);
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
      return Edge(graph_,index_, graph_->adjacency_matrix[index_][aux_index_]);
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
        else if(aux_index_ >= graph_->adjacency_matrix[index_].size()){
          index_++;
          aux_index_=0;
        }
        //Only consider when the first node is smaller so we don't double count
        else if(index_>graph_->adjacency_matrix[index_][aux_index_]){
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
    return edge_iterator(this, 0,0);
  }
  /** Returns the first invalid Edge Iterator of the graph*/
  edge_iterator edge_end() const{
    return edge_iterator(this, num_nodes(), 0) ;
  }
 private:

  std::vector<std::pair<Point,node_value_type>> nodes_;
  std::vector<std::vector<size_type>> adjacency_matrix;

};

#endif // CME212_GRAPH_HPP
