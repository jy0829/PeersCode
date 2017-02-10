#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 *
 * Chloe: chlosmpsn
 * 456Potosi
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
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS

  /** Value types. */
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

  /** Construct an empty graph. */
  Graph(): nodes(), edges(), adjacent() {
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
      return graph_->nodes[index_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(index_ < graph_->size());  
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Returns the value of the node. */
    node_value_type& value(){
        return graph_->nodes[index_].node_value;
    }

    /** Overloading the previous function for constant types. */
    const node_value_type& value() const {
        return graph_->nodes[index_].node_value;
    }

    /** Returns the degree of the node. (number of edges starting from this node) */
    size_type degree() const{
        return graph_->adjacent[index_].size();
    }
    /** The incident iterator begins at the first index of the vector of adjacent nodes.
     */ 
    incident_iterator edge_begin() const{
        return incident_iterator(graph_, index_, 0);
    }
    /** The incident iterator ends after having seen all edges
     *  starting from the current node.
     */
    incident_iterator edge_end() const{
        return incident_iterator(graph_, index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_ && index_ == n.index_);
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
      if (&graph_ < &n.graph_){
          return true;
      }else if (&graph_ == &n.graph_){
          return index_ < n.index_;
      }else{ 
          return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer back to the Graph
    Graph* graph_;
    
    // This node's unique index
    size_type index_;

    /** Private Constructor */
    Node(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {
    }

    /** Helper method to return the appropriate node
     *  The node with index i is also placed in position i in the vector
     *  of nodes (in Graph).
     */

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes.size();
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
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type new_index = size();

    // create a new internal_node object
    internal_node new_node;
    new_node.position = position;
    new_node.node_value = node_value;
    nodes.push_back(new_node);

    adjacent.push_back(std::vector<size_type> ());

    return Node(this, new_index);      
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return(n.graph_ == this && n.index_ < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < num_nodes());
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(index1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(index2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_==e.graph_) {
          if((e.node1() == node1() && e.node2() == node2()) ||
              (e.node1() == node2() && e.node2() == node1()))
              return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if (&graph_ < &e.graph_){
            return true;
        }else if (&graph_ == &e.graph_){
            if (node1().index() < e.node1().index()){
                return true;
            }else if (node1().index() == node1().index()){
                return node2().index() < e.node2().index();
            }else{
                return false;
            }
        }else{ 
            return false;
        }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    
    Graph* graph_;
    size_type index1_; //ID of node 1
    size_type index2_; //ID of node 2
    
    Edge(const Graph* graph, size_type index1, size_type index2) :
        graph_(const_cast<Graph*>(graph)), index1_(index1), index2_(index2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 <= i && i < num_edges());
    internal_edge i_e = edges[i];
    return Edge(this, i_e.index1, i_e.index2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Go through all edges
    for(auto it = adjacent[a.index()].begin(); it != adjacent[a.index()].end(); ++it){
        if((*it) == b.index())
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
    // HW0: YOUR CODE HERE
    assert(has_node(a));
    assert(has_node(b));
    if(!has_edge(a,b)){
        internal_edge new_edge;
        new_edge.index1 = a.index();
        new_edge.index2 = b.index();
        edges.push_back(new_edge);

        adjacent[a.index()].push_back(b.index());
        adjacent[b.index()].push_back(a.index());
    }
    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear();
    edges.clear();
    adjacent.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Dereferencing operator to access the node of index @a index_.
    *  @pre @a index_ < @a graph_->size()
    */
    Node operator*() const {
        return graph_->node(index_);
    }

    /** Operator to get to the next node in the collection
     *  by incrementing @a index_ .
     *  @pre @a index_ < @a graph_->size()
     * */
    node_iterator& operator++(){
        ++index_;
        return (*this);
    }

    /** Test whether this node iterator and @a n are equal.
     */
    bool operator==(const node_iterator& n) const {
      return (graph_ == n.graph_ && index_ == n.index_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
 
    // Pointer back to the Graph
    Graph* graph_;
    
    // This node_iterator's index over the collection of nodes 
    size_type index_;

    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type index) :
        graph_(const_cast<Graph*>(graph)), index_(index) {
            assert(index<=graph->size());
    }   
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /* Create a node_iterator starting from the node of index 0.
   * @post @a index_ == 0
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  /* The end iterator has an index equal to the size of the graph.
   * @post @a index_ == graph_->size()
   */
  node_iterator node_end() const {
    return node_iterator(this, this->size());
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private  totally_ordered<IncidentIterator>{
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const


    /** Dereferencing operator to access the node of index @a index_.
    *  @pre and @post @a index2_ <= @a graph_->adjacent[@a index1_].size()
    */
    Edge operator*() const{
        return Edge(graph_, index1_, graph_->adjacent[index1_][index2_]);
    }

    /** Operator to get to the next edge in the collection
     *  by incrementing @a index2_ .
     *  @pre and @post @a index2_ <= @a graph_->adjacent[@a index1_].size()
     * */
    IncidentIterator& operator++(){
        ++index2_;
        return (*this);
    }

    /** Test whether this incident iterator and @a it are equal.
     */
    bool operator==(const IncidentIterator& it) const{
        return (graph_==it.graph_ and index1_==it.index1_ and index2_==it.index2_);
    } 

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type index1_; //ID of node 1 
    size_type index2_; //Position of the second node in the adjacency matrix
    IncidentIterator(const Graph* graph, size_type index1, size_type index2) :
        graph_(const_cast<Graph*>(graph)), index1_(index1), index2_(index2) {
            assert(index2_<=graph->adjacent[index1].size());
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const


    /** Dereferencing operator to access the edge of the two indices of the nodes.
    */
    Edge operator*() const {
        return Edge(graph_, index1_, graph_->adjacent[index1_][index2_]);
    }

    /** Operator to get to the next egde in the collection
     *  by incrementing first @a index1_, and then @a index2_ once @a index1_ is fixed.
     *
     *  To make sure we do not iterate over every edge twice, we return edges
     *  where the index of the first node is lower or equal than the index of the second node.
     *
     *  @pre @a node1.index() <= @a node2.index()
     *  @pre and @post @a index2_ < @a graph_->adjacent[@a index1_].size()
     * */
    EdgeIterator& operator++(){
        while(index1_ < graph_->adjacent.size()){ //go over index1 nodes
            while(index2_ < graph_->adjacent[index1_].size()-1){ //go over edges to neighbours from index1
                index2_++;
                if(index1_ <= graph_->adjacent[index1_][index2_]){ //only pick when node1<=node2
                    return *this;
                }
            }
            index1_++;
            index2_=0;
        }
        return *this;
    }

    /** Test whether this edge iterator and @a e are equal.
     */
    bool operator==(const EdgeIterator& e) const {
        return (graph_==e.graph_ and index1_==e.index1_ and index2_==e.index2_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type index1_;
    size_type index2_;
    EdgeIterator(const Graph* graph, size_type index1, size_type index2) :
        graph_(const_cast<Graph*>(graph)), index1_(index1), index2_(index2) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const


    /** The edge_iterator starts at the first edge starting from the first node
     *  
     *  @pre 0 <= @a graph_->adjacency[0][0]
     *  @post @a index1_ == 0
     *  @post @a index2_ == 0
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0, 0);
    }
    /** The end iterator is with @a index1_ beyond the size of the graph
     *  @post @a index1_ == @a graph_->size()
     */
    edge_iterator edge_end() const {
        return EdgeIterator(this, num_nodes(), 0);
    }



 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
      Point position;
      node_value_type node_value;
  };

  struct internal_edge {
      size_type index1;
      size_type index2;
  };

  std::vector<internal_node> nodes; 
  std::vector<internal_edge> edges;
  std::vector<std::vector<size_type>> adjacent;
};

#endif // CME212_GRAPH_HPP
