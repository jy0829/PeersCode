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
  struct node_data;
  struct edge_data;
  std::vector<node_data> node_list;
  std::vector<edge_data> edge_list;
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
      return fetch_position(); 

    }
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    

    /** return reference to the value associated with the node */
    node_value_type& value(){
      return ((this->g)->node_list[this->node_id]).val; 
    }
    
    /** return const reference to the value associated with the node */
    const node_value_type& value() const {
      return ((this->g)->node_list[this->node_id]).val; 

    }

    /** return the number of edges connected to the calling node */
    size_type degree() const {
      return g->node_adj_list[node_id].size();
    }

    /** return the iterator pointing to the beginning of the collection of
     * edges that are adjacent to the calling node. 
     * @pre: cannot be dereferenced if the calling node is isolated */

    incident_iterator edge_begin() const {
      return IncidentIterator(g, node_id, g->node_adj_list[node_id].begin());
    }

    /** return the iterator pointing to one past the end of the collection
      * of edges adjacent to the calling node */

    incident_iterator edge_end() const {
      return IncidentIterator(g, node_id, g->node_adj_list[node_id].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->g == n.g && this-> node_id == n.index()){
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
      // HW0: YOUR CODE HERE
      // This way of comparison will maintain the order of nodes with in a
      // graph; for comparison of nodes between two graphs the ordering may
      // change
      if (this->g != n.g){
	return (this->g < n.g);
      }
      else{
	return (this-> node_id < n.index());
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* g;
    size_type node_id;
    
    Node(const Graph* graph, size_type uid)
        : g(const_cast<Graph*>(graph)), node_id(uid) {
    }

    size_type fetch_id() const {
        return this->node_id;
    }

    const Point& fetch_position() const {
        return ((this->g)->node_list[this->node_id]).p;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_list.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()){
    node_data new_node = {position, value};
    node_list.push_back(new_node);
    return Node(this, node_list.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.g == this && n.node_id < node_list.size()){
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
    // HW0: YOUR CODE HERE
    //every time this method is called a new Node object is returned, which should be ok (Node is basically pointer)
    if (node_list.size() > i){
        return Node(this, i); 
    }
    else{
        return Node();        // Invalid node
    }
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
      // HW0: YOUR CODE HERE
      return Node(g, _node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g, _node2_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->g == e.g && ((_node1_id == e._node1_id && _node2_id == e._node2_id) || (_node1_id == e._node2_id && _node2_id == e._node1_id)));
    }


    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // This way of comparison will maintain the order of edges with in a
      // graph; for comparison of edges of two graphs the ordering may
      // change
      if (this->g != e.g){
	return (this->g < e.g);
      }
      else{

      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* g;
    size_type _node1_id;
    size_type _node2_id;

    Edge(const Graph* graph, size_type node1_id, size_type node2_id)
      : g(const_cast<Graph*> (graph)), _node1_id(node1_id), _node2_id(node2_id) {
      } 
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
    return Edge(this, edge_list[i].node1_id, edge_list[i].node2_id); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (find_edge_index(a, b) != edge_list.size()){
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
    if (find_edge_index(a, b) != edge_list.size()){
      return edge(find_edge_index(a, b));
    }
    edge_data new_edge = {a.index(), b.index()};
    edge_list.push_back(new_edge);
    insert_into_adj_list(a, b, edge_list.size() - 1);
    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_list.clear();
    edge_list.clear();
    node_adj_list.clear();
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
    /** return a Node object corresponding to the node referred to by the Iterator.
     * same specification as the dereference operator for a standard forwardIterator*/
    Node operator*() const {
      return Node(this->_g, this->_node_id);
    }

    /** same specification as the ++ operator for forwardIterator */
    NodeIterator& operator++(){
      this->_node_id++;
      return *this;
    }

    /** return if this iterator and another iterator is equal. Two iterators
     * are equal if they refer to the same node. Must be in constant time. */

    bool operator==(const NodeIterator& n) const {
      return Node(this->_g, this->_node_id) == *n;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* _g;
    size_type _node_id;
    
    NodeIterator(const Graph* g, size_type ind)
      : _g(const_cast<Graph*>(g)), _node_id(ind) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** This returns an iterator pointing to the first element; should not be
   * dereferenced if the graph is emppty */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** This returns an iterator pointing one past the end of the node vector;
   * should not be dereferenced */
  node_iterator node_end() const {
    return NodeIterator(this, node_list.size());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** returns an Edge object corresponding the the edge the IncidentIterator
     * is pointing to. Does not have to work if the IncidentIterator if it
     * is pointing to one past the end. The node1() function of the returned Edge 
     * must return the incident node. */

    Edge operator*() const {
      size_type edge_id = *_incident_iter;
      size_type temp1 = (*_g).edge(edge_id).node1().index();
      size_type temp2 = (*_g).edge(edge_id).node2().index();
      if (temp1 == _node1_id){
	return Edge(_g, temp1, temp2);
      }
      else {
	return Edge(_g, temp2, temp1);
      }
    }

    /** returns an iterator that points to the next incident edge; if the
      * iterator is already one past the end then return the same iterator
      * (this follows from the vector::iterator) */

    IncidentIterator& operator++() {
      _incident_iter++;
      return *this;
    }

    /**  two IncidentIterators are equal if their all their members are equal */
    bool operator==(const IncidentIterator& iter) const {
      return (_g == iter._g && _node1_id == iter._node1_id && _incident_iter == iter._incident_iter);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* _g;
    size_type _node1_id;
    std::vector<size_type>::iterator _incident_iter;

    IncidentIterator(const Graph* graph, size_type node1, std::vector<size_type>::iterator iter)
      : _g(const_cast<Graph*>(graph)), _node1_id(node1), _incident_iter(iter) {
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
    /** the dereference operator should return an Edge object corresponding
     * to the edge the iterator is pointing to in the graph. If the iterator is
     * one past last element this doesn't have to work */ 

    Edge operator*() const {
      return Edge(_g, _g->edge_list[_edge_id].node1_id, _g->edge_list[_edge_id].node2_id);
    }

    /** this operator should return the iterator point to the next edge. */

    EdgeIterator& operator++() {
      _edge_id++;
      return *this;
    }

    /** this returns true if the two edges come from the same graph and has the
     * same index in edge_list */

    bool operator==(const EdgeIterator& iter) const {
      return (_g == iter._g && _edge_id == iter._edge_id);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* _g;
    size_type _edge_id;

    EdgeIterator(const Graph* graph, size_type id)
      : _g(const_cast<Graph*>(graph)), _edge_id(id) {
      }
  };

  // HW1 #5: YOUR CODE HERE
  /** return the EdgeIterator pointing to the start of edge_list */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** return the EdgeIterator for one past the last edge */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_list.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<std::vector<size_type>> node_adj_list;

  struct node_data{
    Point p;
    node_value_type val;
  };

  struct edge_data{
    size_type node1_id;
    size_type node2_id;
  };


  void insert_into_adj_list(const Node& a, const Node& b, size_type index){
    size_type max_ind = std::max(a.index(), b.index());
    if (node_adj_list.size() <= max_ind){
      while (max_ind >= node_adj_list.size()){
	std::vector<size_type> placeholder;
	node_adj_list.push_back(placeholder);
      }
    }
    node_adj_list[a.index()].push_back(index);
    node_adj_list[b.index()].push_back(index);
  } 

  /** find_edge_index is a function that takes in two nodes and finds the 
   * index for the corresponding edge. If such edge does not exist in the
   * graph, the index of one past last edge will be returned
   * @pre: @a n1 and @a n2 must be valid nodes in the graph
   * Note: Complexity is not O(1); it is O(n1.degree()) */

  size_type find_edge_index(const Node& n1, const Node& n2){
    if (node_adj_list.size() <= n1.index()) return edge_list.size();
    auto start = node_adj_list[n1.index()].begin(); 
    auto end = node_adj_list[n1.index()].end();
    for (auto it = start; it!= end; ++it){
      if (edge_list[*it].node1_id == n2.index() || edge_list[*it].node2_id == n2.index()){
	return *it;
      }
    }
    return edge_list.size();
  }
};

#endif // CME212_GRAPH_HPP
