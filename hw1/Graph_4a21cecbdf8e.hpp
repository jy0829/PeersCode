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

  /** Construct an empty graph. */
  Graph(){   
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
      const Point& pos = graph_-> node_list[index_];
      return pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
      node_value_type& val = graph_ -> node_value_list[index_];
      return val;
    }
 
    const node_value_type& value() const{
      const node_value_type& val = graph_ -> node_value_list[index_];
    }

    size_type degree() const{
      return graph_ -> adj_list[index_].size();
    }

    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,index_,0);
    }

    incident_iterator edge_end() const{
      return IncidentIterator(graph_,index_,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if(graph_ == n.graph_ and index_ == n.index()){
      return true;
    } else{          // Quiet compiler warning
      return false;
    }
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
      if(graph_ == n.graph_ and index_ < n.index()){
      return true;
    } else{          // Quiet compiler warning
      return false;
    }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type index_;

    Node(const Graph* graph, size_type indexu)
       :graph_(const_cast<Graph*>(graph)), index_(indexu){
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    size_type n_nodes = node_list.size();
    return n_nodes;
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

    node_value_list.push_back(node_value);
    node_list.push_back(position);
    node_in = node_list.size()-1;
    adj_list.push_back({});
    // HW0: YOUR CODE HERE
   (void) position;      // Quiet compiler warning
    return Node(this,node_in);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if(n.graph_ == this ){
      return true;
    } else{      // Quiet compiler warning
    return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const{
    assert(i < size());
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    return Node(this, i);        // Invalid node
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
  class Edge :private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      Node temp1 = graph_ -> edge_list[edge_id][0];
      return temp1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      Node temp2 = graph_ -> edge_list[edge_id][1];
      return temp2;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (node1() == e.node1() and node2() == e.node2()){
      return true;
      } else if (node1() == e.node2() and node2() == e.node1()){
      return true;
      } else {           // Quiet compiler warning
      return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (node1() < e.node1() and node2() < e.node2()){
      return true;
      } else if (node1() < e.node2() and node2() < e.node1()){
      return true;
      } else {           // Quiet compiler warning
      return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type edge_id;
    
    Edge(const Graph* graph, size_type edgeid)
        :graph_(const_cast<Graph*>(graph)), edge_id(edgeid){
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    return Edge(this, i);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for(size_type j=0; j<num_edges(); ++j){
      if(edge_list[j][0] == a and edge_list[j][1] == b){
        return true;
      } else if (edge_list[j][1] == a and edge_list[j][0] == b){
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

    // Comments to HW0, somehow when I un-comment the lines, the edges disappear on the graph.
    // But when I comment these lines out, with no preconditions or whatever, the edges come back.
    // my guess is that I used the wrong operator "==", but I could not quite figure out a solution at the moment.
    // But I left it to you just to show the algorithm I proposed.

    // for(size_type j=0; j<num_edges(); ++j){
    //    if(edge_list[j][0].position() == a.position() and edge_list[j][1].position() == b.position()){
    //    edge_in = j;
    //  } else if (edge_list[j][1].position() == a.position() and edge_list[j][0].position() == b.position()){
    //    edge_in = j;
    //  } else {
        size_type a_idx = a.index();
        size_type b_idx = b.index();
        
        adj_list[a_idx].push_back(b_idx);
        adj_list[b_idx].push_back(a_idx);

        std::vector <Node> temp_t = {a,b};
        edge_list.push_back(temp_t);
    
        edge_in = edge_list.size()-1;
        (void) a, (void) b;   // Quiet compiler warning
        // Invalid Edge
     // }
     // }
      return Edge(this,edge_in);
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
   node_list = {};
   edge_list = {};
    // HW0: YOUR CODE HERE
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

    Node operator*() const{
      return graph_ -> node(nodep_);
    }

    NodeIterator& operator++(){
      if(nodep_ < graph_ -> size()){
      nodep_ = nodep_ + 1;
      }
      return *this;
    }

    bool operator==(const NodeIterator& ni) const{
      return nodep_ == ni.nodep_;
    }

   private:
    friend class Graph;
    
    Graph* graph_;
    size_type nodep_;
    
    NodeIterator(const Graph* graph, size_type nodep):graph_(const_cast<Graph*>(graph)),nodep_(nodep){
    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  node_iterator node_begin() const{
    return NodeIterator(this,0);  
  }

  node_iterator node_end() const{
    return NodeIterator(this,node_list.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
    
    Edge operator*() const{
    return graph_ -> edge(eid_);
    }

    IncidentIterator& operator++(){
    p_nidx_ = p_nidx_ + 1;
    return *this;
    }

    bool operator==(const IncidentIterator& ii) const{
    return(graph_ == ii.graph_ and eid_ == ii.eid_ and p_nidx_ == ii.p_nidx_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type eid_;
    size_type p_nidx_;

    IncidentIterator(const Graph* graph, size_type eid, size_type p_nidx): graph_(const_cast<Graph*>(graph)),eid_(eid),p_nidx_(p_nidx){
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
    Edge operator*() const{
      return Edge(graph_, eid_);
    }

    EdgeIterator& operator++(){
      eid_ = eid_ + 1;
      return *this;
    }

    bool operator==(const EdgeIterator& ei) const{
      return (graph_ == ei.graph_ and eid_ == ei.eid_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type eid_;

    EdgeIterator(const Graph* graph, size_type eid):
      graph_(const_cast<Graph*>(graph)), eid_(eid){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_list.size());
  }
 private:
   size_type node_in;
   size_type edge_in;
   std::vector<Point> node_list;
   std::vector<node_value_type> node_value_list;
   std::vector<std::vector<Node>> edge_list;
   std::vector<std::vector<size_type>> adj_list;
   
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
