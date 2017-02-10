#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

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
  class Node : private totally_ordered<Node>{
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
      nid_ = -1; //Invalid node id
      g_ = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g_->nodes[nid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(nid_ < g_->size());
      return nid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /** Returns the value of this node
     *  @pre @a nid_ should be a valid node id i.e. 0<= @a nid_ < num_nodes()
     *  @pre @a g_ should not be NULL 
     *  @return value of this node 
     */
    node_value_type& value(){
       return g_->nodes[nid_].value;
    }
    /** Returns the value of this node
     * @pre @a nid_ should be a valid node id i.e. 0<= @a nid_ < num_nodes()
     * @pre @a g_ should not be NULL 
     * @return value of this node
     */
    const node_value_type& value() const{
      return g_->nodes[nid_].value;
    }
    /** Returns the degree of this node
     *  @pre @a nid_ should be a valid node id i.e. 0<= @a nid_ < num_nodes()
     *  @pre @a g_ should not be NULL 
     */
    size_type degree() const{
      return g_->adjacency[nid_].size();
    }
    /** Returns incident_iterator to the first edge of this node 
     *  @pre @a nid_ should be a valid node id i.e. 0<= @a nid_ < num_nodes()
     *  @pre @a g_ should not be NULL 
     */
    incident_iterator edge_begin() const{
    	return IncidentIterator(g_,0,nid_);
    }
    /** Returns incident_iterator to one past the last edge of this node
     *  @pre @a nid_ should be a valid node id i.e. 0<= @a nid_ < num_nodes()
     *  @pre @a g_ should not be NULL 
     */
    incident_iterator edge_end() const{
        return IncidentIterator(g_,degree(),nid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      if((this->g_ == n.g_) && (this->nid_==n.nid_))
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
      // HW0: YOUR CODE HERE
      (void) n;           // Quiet compiler warning
      if(this->g_<n.g_)
        return true;
      else if((this->g_==n.g_) && (this->nid_<n.nid_))
        return true;
      else
        return false;
    }
   
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    size_type nid_;     //node id
    Graph* g_;          //pointer to the graph
    
    //private constructor
    Node(const Graph* g, size_type nid): nid_(nid),g_(const_cast<Graph*>(g)){
    }

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
   * @param[in] value The new node's value. If not specified, uses the default value.
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * 
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value= node_value_type()){
    size_type old_num_nodes = num_nodes();

    internal_node new_node;
    new_node.position = position;
    new_node.value = value;
    new_node.index = old_num_nodes;

    nodes.push_back(new_node);
    
    //add empty vector to adjacency list
    adjacency.push_back(std::vector<std::pair<size_type,size_type>>());
    
    return Node(this, old_num_nodes);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    if(n.g_==this && n.index()<size())
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
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    assert(0<=i && i<num_nodes());
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
      eid_ = -1;
      g_ = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return g_->edges[eid_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return g_->edges[eid_].node2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if((this->g_ == e.g_) && (this->eid_==e.eid_))
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if(this->g_<e.g_)
        return true;
      else if((this->g_==e.g_) && (this->eid_<e.eid_))
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type eid_;
    Graph* g_;
    Edge(const Graph* g, size_type eid): eid_(eid),g_(const_cast<Graph*>(g)){
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
    (void) i;             // Quiet compiler warning
    assert(0<=i && i<num_edges());
    return Edge(this,i);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    for(size_type i=0; i<adjacency[a.index()].size(); i++)
    {
	if(adjacency[a.index()][i].first == b.index())
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
    (void) a, (void) b;   // Quiet compiler warning

    assert(has_node(a));
    assert(has_node(b));

    //Search if edge exists
    for(size_type i=0; i<adjacency[a.index()].size(); i++)
    {
	if(adjacency[a.index()][i].first == b.index())
	    return Edge(this,adjacency[a.index()][i].second);
    }
	
    //if edge not found, add edge
    size_type old_num_edges = num_edges();
    internal_edge new_edge;
    new_edge.index = old_num_edges;
    new_edge.node1 = a;
    new_edge.node2 = b;
    edges.push_back(new_edge);
        
    //Add edge to adjacency list
    adjacency[a.index()].push_back(std::make_pair(b.index(),old_num_edges));
    adjacency[b.index()].push_back(std::make_pair(a.index(),old_num_edges));

    return Edge(this, old_num_edges);
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
    adjacency.clear();
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
	uid_=-1;
 	g_=nullptr;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Returns the node which this iterator points to 
     *  @pre @a uid_ should be a valid node id, i.e. 0<= @a uid_ < num_nodes()
     *  @pre @a g_ should not be null
     */
    Node operator*() const{
        return g_->node(uid_);
    }
    /** Returns the next NodeIterator 
     * @pre @a uid_ should be a valid node id, i.e. 0<= @a uid_ < num_nodes()
     */
    NodeIterator& operator++(){
       uid_ = uid_+1;
       return *this;
    }
    /* Returns true if this NodeIterator and @a ni are equal */
    bool operator==(const NodeIterator& ni) const{
        bool value =  ((g_ == ni.g_) && (uid_ == ni.uid_));
	return value;
    }
    
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type uid_;     //node id 
    Graph* g_; 		//pointer to the graph

    //private constructor
    NodeIterator(const Graph* g,size_type uid): uid_(uid),g_(const_cast<Graph*>(g)){}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
  /** Returns node_iterator which points to the first node */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  /** Returns node_iterator which points to one past the last node */
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Returns the edge this IncidentIterator points to 
     *  @pre @a g_ should not be NULL
     *  @pre @a nid_ should be a valid node id, i.e. 0<= @a nid_ < num_nodes()
     *  @pre @a eid_ should be a valid index in the adjacency list of node with id @a nid_ , 
     *       i.e. 0<= @a eid_ < node(@a nid_).degree()
     */
    Edge operator*() const{
        return Edge(g_,g_->adjacency[nid_][eid_].second);
    }
    /* Returns the next IncidentIterator */
    IncidentIterator& operator++(){
	eid_++;
	return *this;
    }
    /* Returns true if this IncidentIterator and @a ii are equal */
    bool operator==(const IncidentIterator& ii) const{
	return ((g_==ii.g_) && (eid_==ii.eid_) && (nid_==ii.nid_));
    }
     
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type eid_; //index of edge in adjacency list of node nid_
    size_type nid_; //node id
    Graph* g_;      //pointer to the graph

    //private constructor
    IncidentIterator(const Graph* g, size_type eid, size_type nid): eid_(eid),nid_(nid),g_(const_cast<Graph*>(g)){}
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
	eid_=-1;
        g_=nullptr;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Returns the edge this EdgeIterator points to
     *  @pre @a eid_ should be a valid edge id, i.e. 0<= @a eid_ < num_edges() 
     *  @pre @a g_ should not be NULL
     */
    Edge operator*() const{
	return g_->edge(eid_);
    }
    /* Returns the next EdgeIterator */
    EdgeIterator& operator++(){
	eid_++;
	return *this;
    }
    /* Returns true if this and @a ei are equal */
    bool operator==(const EdgeIterator& ei) const{
 	return (g_==ei.g_) && (eid_==ei.eid_);
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type eid_; //edge id 
    Graph* g_;      //pointer to the graph

    //private constructor
    EdgeIterator(const Graph* g, size_type eid): eid_(eid),g_(const_cast<Graph*>(g)){}
    
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /* Returns edge_iterator which points to the first edge */
  edge_iterator edge_begin() const{
	return EdgeIterator(this, 0);
  }
  /* Returns edge_iterator which points to one past the last edge */
  edge_iterator edge_end() const{
	return EdgeIterator(this, num_edges());
  }
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  struct internal_node{
	size_type index; 
        node_value_type value;
  	Point position;
  };

  struct internal_edge{
	size_type index;
	Node node1;
	Node node2;
  };

  std::vector<internal_node> nodes;
  std::vector<internal_edge> edges;
  std::vector<std::vector<std::pair<size_type,size_type>>> adjacency;
};

#endif // CME212_GRAPH_HPP
