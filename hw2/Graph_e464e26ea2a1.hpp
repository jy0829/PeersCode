#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <list>
#include <unordered_map>

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
  struct node_data;
  struct edge_data;

  using adj_list_outer_it = typename std::vector<std::vector<edge_data>>::iterator;
  using adj_list_outer_it_const = typename std::vector<std::vector<edge_data>>::const_iterator;

  using adj_list_inner_it = typename std::vector<edge_data>::iterator;
  using edge_value_it = typename std::list<E>::iterator;

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
  using edge_value_type = E;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    next_node_id = 0;
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
    Point& position() {
      return fetch_position();
    }
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch_position(); 

    }
    /** Return this node's index, a number in the range [0, graph_size). 
     *  complexity: on average constant but could be linear in worst case
     */

    size_type index() const {
      // HW0: YOUR CODE HERE
      return g->node_id_to_index[node_id];
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    

    /** return reference to the value associated with the node */
    node_value_type& value(){
      return (g->node_list[index()]).val; 
    }
    
    /** return const reference to the value associated with the node */
    const node_value_type& value() const {
      return (g->node_list[index()]).val; 

    }

    /** return the number of edges connected to the calling node */
    size_type degree() const {
      return g->node_adj_list[index()].size();
    }

    /** return the iterator pointing to the beginning of the collection of
     * edges that are adjacent to the calling node. 
     * @pre: cannot be dereferenced if the calling node is isolated */

    incident_iterator edge_begin() const {
      return IncidentIterator(g, g->node_adj_list.begin() + index(), g->node_adj_list[index()].begin());
    }

    /** return the iterator pointing to one past the end of the collection
      * of edges adjacent to the calling node */

    incident_iterator edge_end() const {
      return IncidentIterator(g, g->node_adj_list.begin() + index(), g->node_adj_list[index()].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->g == n.g && this-> node_id == n.node_id){
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
	return (this-> node_id < n.node_id);
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

    Point& fetch_position() const {
        return (g->node_list[index()]).p;
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
    size_type new_node_id = next_node_id;
    next_node_id ++;
    node_id_to_index.insert({new_node_id, node_list.size()});
    node_data new_node = {new_node_id, position, value};
    node_list.push_back(new_node);
    node_adj_list.push_back(std::vector<edge_data> ());
    return Node(this, new_node_id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1) on average, could be O(n) in worst case
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.g == this && node_id_to_index.find(n.node_id) != node_id_to_index.end()){
        return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1)
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //every time this method is called a new Node object is returned, which should be ok (Node is basically pointer)
    return Node(this, node_list[i].node_id);        // Invalid node
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

    /** returns the data associated with the edge
     *  complexity is O(node1.degree())
     */
    edge_value_type& value() {
      return g->get_edge_value(node1(), node2());
    }

    const edge_value_type& value() const {
      return g->get_edge_data_value(node1(), node2());
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
	return std::make_pair(std::min(_node1_id, _node2_id), std::max(_node1_id, _node2_id)) < std::make_pair(std::min(e._node1_id, e._node2_id), std::max(e._node1_id, e._node2_id));
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
    return edge_values.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less, current implementation is O(a.degree())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (auto it = a.edge_begin(); it!= a.edge_end(); ++it){
      if ((*it).node2() == b){
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
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less. Current implementation is O(a.degree()).
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if (has_edge(a, b)){
      return Edge(this, a.node_id, b.node_id);
    }

    edge_values.push_front(value);
    node_adj_list[a.index()].push_back({a.node_id, b.node_id, edge_values.begin()});
    node_adj_list[b.index()].push_back({b.node_id, a.node_id, edge_values.begin()});
    return Edge(this, a.node_id, b.node_id);
  }

  /** this function removes the specified node from node_list,
   *    node_id_to_index, and all edges connected to that node
   *  @return nothing meaningful
   *  @pre @a n must be a valid node of the graph
   *  @post all outstanding IncidentIterator, and
   *    EdgeIterator objects are invalidated; all Node and Edge
   *    objects that corresponds to the removed node and edges
   *    are invalidated; all other Edge and Node objects remain
   *    valid. For any remaining valid nodes n', if old n'.index() < 
   *    old n.index(), then new n'.index() == old n'.index(); 
   *    if old n'.index() > old n.index(), then new n'.index() == old n'.index - 1.
   *    The NodeIterators pointing to nodes before the removed node
   *    in @a node_list are still valid, other NodeIterators are also
   *    invalidated by this function.
   *
   *  run-time complexity is O(num_nodes() + n.degree()^2); since 
   *    we are assuming that the graph is sparce, n.degree()^2 << num_nodes
   *    and the complexity is roughly O(num_nodes()).
   */

  size_type remove_node(const Node& n) {
    size_type ind = n.index();
    while (!node_adj_list[ind].empty()){
      remove_edge(n, Node(this, node_adj_list[ind].back().node2_id));
    }
    node_list.erase(node_list.begin() + ind);
    node_adj_list.erase(node_adj_list.begin() + ind);
    
    // remove index of removed node and reassign index to remaining nodes
    node_id_to_index.erase(n.node_id);
    for (size_type it = ind; it != node_list.size(); ++it){
      node_id_to_index[node_list[it].node_id] = it;
    }
    return 0;
  }

  /** this function takes in a node_iterator and removes the corresponding node.
   *  @param[in] n_it node_iterator of the graph
   *  @return a node_iterator that points to the node one after the removed node
        in the old @a node_list (or one past the end if the removed node was the
   *    last node): *@return == old *++ @a n_it
   *  @pre n_it must be a valid, dereferencable node_iterator of the graph
   *  @post same as remove_node(const Node& n)
   *  complexity is the same as remove_node(const Node& n)
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  /** this remove_edge takes in two nodes and removes the edge 
   *  connecting them from the graph
   *  @return nothing meaningful
   *  @pre @a a and @a b must be valid nodes; has_edge(a, b) == true
   *  @post all outstanding EdgeIterator and IncidentIterator objects
   *    are invalidated after this function; all Edge objects corresponding
   *    to the removed edge are invalidated, but all other Edge objects
   *    remain valid. All invariants between the data structure of the
   *    graph should also hold.
   *  Run time complexity is O(a.degree() + b.degree())
   */

  size_type remove_edge(const Node& a, const Node& b) {
    for (auto it = node_adj_list[a.index()].begin(); it != node_adj_list[a.index()].end(); ++it){
      if ((*it).node2_id == b.node_id){
	edge_values.erase((*it).val_it);
	node_adj_list[a.index()].erase(it);
	break;
      } 
    }

    for (auto it = node_adj_list[b.index()].begin(); it != node_adj_list[b.index()].end(); ++it){
      if ((*it).node2_id == a.node_id){
	node_adj_list[b.index()].erase(it);
	break;
      }
    }
    return 0;
  }

  /** this remove_edge takes in an Edge object and removes
   *  the corresponding edge from the graph
   *  @return nothing meaningful
   *  @pre @a e must be a valid Edge object
   *  @post same as post-condition of remove_edge(const Node& a, const Node& b)
   *  complexity is the same as remove_edge(const Node& a, const Node& b)
   */

  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** this remove_edge takes in an EdgeIterator and removes
   *  the corresponding edge from the graph
   *  @return the first edge_iterator of the new graph after
   *    the edge is removed
   *  @pre @a e_it must be dereferencable
   *  @post same as the post-condition of remove_edge(const Edge& e)
   *  complexity is the same as remove_edge(const Edge e)
   */

  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return edge_begin();
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects and all iterators
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_list.clear();
    node_adj_list.clear();
    edge_values.clear();
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
      return Node(_g, _g->node_list[_index].node_id);
    }

    /** same specification as the ++ operator for forwardIterator */
    NodeIterator& operator++(){
      _index++;
      return *this;
    }

    /** return if this iterator and another iterator is equal. Two iterators
     * are equal if they refer to the same node. Must be in constant time. */

    bool operator==(const NodeIterator& n) const {
      return _g == n._g && _index == n._index;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* _g;
    size_type _index;
    
    NodeIterator(const Graph* g, size_type ind)
      : _g(const_cast<Graph*>(g)), _index(ind) {
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
      edge_data e = *_incident_iter;
      return Edge(_g, e.node1_id, e.node2_id);
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
      return (_g == iter._g && _list_iter == iter._list_iter && _incident_iter == iter._incident_iter);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    friend class EdgeIterator;
    Graph* _g;
    adj_list_outer_it _list_iter;
    adj_list_inner_it _incident_iter;

    IncidentIterator(const Graph* graph, adj_list_outer_it list_iter, adj_list_inner_it iter)
      : _g(const_cast<Graph*>(graph)), _list_iter(list_iter), _incident_iter(iter) {
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
      return Edge(_g, (*_out_it)[_ind].node1_id, (*_out_it)[_ind].node2_id);
    }

    /** this operator should return the iterator point to the next edge.
     *  should always be dereferencable if it is not pointing to one past
     *  the end (for which _out_it must be pointing to
     *  one past the end of node_adj_list 
     *  
     *  Will skip over edges for which node1_id > node2_id (node1 is the incident node of the edge)
     *  so each edge will be traversed only once
     */

    EdgeIterator& operator++() {
      if (_out_it == _g->node_adj_list.end()) return *this;
      ++ _ind;
      while ((*_out_it).size() <= _ind || (*_out_it)[_ind].node1_id > (*_out_it)[_ind].node2_id){
        if ((*_out_it).size() <= _ind){
          ++_out_it;
	  _ind = 0;
	  if (_out_it == _g->node_adj_list.end()) break;
        }
	else {
	  ++ _ind;
	}
      }
      return *this;
    }

    /** this returns true if the two edges come from the same graph and 
     *  point to the same edge_data */

    bool operator==(const EdgeIterator& iter) const {
      return (_g == iter._g && _out_it == iter._out_it && _ind == iter._ind);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* _g;
    adj_list_outer_it_const _out_it;
    size_type _ind;

    EdgeIterator(const Graph* graph, adj_list_outer_it_const out, size_type ind) 
      : _g(const_cast<Graph*>(graph)), _out_it(out), _ind(ind) {
      }
  };

  // HW1 #5: YOUR CODE HERE
  /** return the EdgeIterator pointing to the start of edge_list */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, node_adj_list.begin(), 0);
  }

  /** return the EdgeIterator for one past the last edge */
  edge_iterator edge_end() const {
    return EdgeIterator(this, node_adj_list.end(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  /** The data structures storing the data of the graph are
   *  @a node_list is a vector of nodes storing the nodes and
	data associated with every node
   *  @a node_id_to_index is a hash map for looking up node_id with the current
        index of the corresponding node in @a node_list
   *  @a node_adj_list is an adjacency list of the graph storing all
	the data for the edges; each edge's data is stored exactly
	twice in @a node_adj_list
   *  @a edge_values is a list of data for the edges. The data can be accessed
        via iterators.
  
   *  There are identification numbers associated with each node: its position in
        the @a node_list vector, which is node.index(), and node_id, which is a member
        of every Node object and is invariant for each node until destroyed. The mapping
        from index to node_id is stored in @a node_list, and the reverse mapping is
        stored in @a node_id_to_index
   *
   *  Invariants of this representation:
	@a node_adj_list.size() == @a node_list.size();

	For all 0 <= i < j < @a node_list.size(), if edge(node(i), node(j))
	  is in @a node_list[i], it is also in @a node_list[j]; if it is not
	  in @a node_list[i], it is not in @a node_list[j].

	For all 0 <= i < @a node_list.size(), for all 0 <= j  < @a node_adj_list[i].size(),
	  node_adj_list[i][j].node1_id == i

	2 * edge_values.size() == number of edge_data objects in node_adj_list.

	For each edge_data object e, its member @a e.val_it must be valid and dereferencable.

        node(i).index() == i for all i with 0 <= i < num_nodes().
 
        for any valid Node object n, node(n.index()) == n and if e is in node_adj_list[n.index()],
          e.node1_id == n.node_id

        for any valid node n, next_node_id > n.node_id
   */


  std::vector<std::vector<edge_data>> node_adj_list;
  std::vector<node_data> node_list;
  std::unordered_map<size_type, size_type> node_id_to_index;
  std::list<edge_value_type> edge_values;
  size_type next_node_id;
  
  struct node_data{
    size_type node_id;
    Point p;
    node_value_type val;
  };

  struct edge_data{
    size_type node1_id;
    size_type node2_id;
    edge_value_it val_it;
  };


  /** This function takes in two nodes and finds the 
   * data for the corresponding edge.
   * @pre: @a n1 and @a n2 must be valid nodes in the graph
   * Note: Complexity is not O(1); it is O(n1.degree()) */


  edge_value_type& get_edge_value(const Node& n1, const Node& n2){
    auto start = node_adj_list[n1.index()].begin(); 
    auto end = node_adj_list[n1.index()].end();
    for (auto it = start; it!= end; ++it){
      if ((*it).node2_id == n2.node_id){
	return *((*it).val_it);
      }
    }
    // returns a random value, this will never run by pre-condition
    return *((*start).val_it);
  }
};


#endif // CME212_GRAPH_HPP
