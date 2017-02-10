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

  /** Synonym for template parameter of graph, which allows nodes
  * to support a user specified value of type node_value_type
  */
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
  Graph() {
    numEdges = 0;
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph->nodes[this->index()].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return nid;
    }

    /** A function that allows the user to see and enter
     * a value for a node of type node_value_type.
     * @return node_value_type& represents value of node
     * @post   value is modifiable by user
     *
     * Complexity: O(1)
     */
    node_value_type& value(){
      return graph->nodes[this->index()].second;
    }

    /** A function that allows the user to see
     * a value for a node of type node_value_type.
     * @pre const node calls the function
     * @return const node_value_type& represents value of node
     *
     * Complexity: O(1)
     */
    const node_value_type& value() const{
      return graph->nodes[this->index()].second;
    }

    /** Returns the degree (number of adjacent nodes to n) of 'this' node
     * @return  a value of size_type indicating the degree of 'this' node
     *
     * Complexity: O(1)
     */
    size_type degree() const{
      return graph->adj_list[this->index()].size();
    }

    /** Returns an incident iterator which allows the user to begin
     * iterating over all of the adjacent nodes to a particular node in the graph
     * @return  an iterator of type incident_iterator
     *
     * The function uses the index of the current node and begins with
     * the first neighboring node.
     *
     * Complexity: O(1)
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph, this->index(), 0);
    }

    /** Returns an incident iterator which indicates to the user that
     *  all adjacent nodes to a particular node in the graph have been iterated over.
     * @return  an iterator of type incident_iterator
     *
     * The function uses the index of the current node and the degree of
     * the current node to indicate that iteration should stop.
     *
     * Complexity: O(1)
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph, this->index(), degree());
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (nid == n.index() && graph == n.graph);
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
      return ((graph < n.graph) || (nid < n.index() && graph == n.graph));
    }

   private:
    friend class Graph;     // Allow Graph to access Node's private member data and functions.

    /** Constructs a valid node within a specified graph
     * @param curr_graph  pointer to graph in which node is added to
     * @param i           id assigned to node
     *
     * @pre               @a graph != nullptr
     * @pre               @nid < number of nodes in @a graph
     */
    Node(const Graph* curr_graph, size_type i) {
       nid = i;
       graph = const_cast<Graph*>(curr_graph);
       assert(graph != nullptr);
       //assert(nid < (graph->num_nodes()+1));
    }

	  size_type nid;   // holds node id
	  Graph* graph;    // holds graph in which node lives
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    nodes.push_back(std::make_pair(position, val));	            // adds node to graph
    std::vector<size_type> empty_list;
    adj_list.push_back(empty_list);                             // adds node to adjacency list
    return Node(this, (size()-1));
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    size_type numNodes = this->size();
    return (!(numNodes < n.index()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == itihs
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph, nodeOne);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph, nodeTwo);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(node1() == e.node1() && node2() == e.node2() && graph == e.graph){
        return true;
      }
      if(node1() == e.node2() && node2() == e.node1() && graph == e.graph){
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
      return ((graph < e.graph) || (node1() < e.node1() && graph == e.graph));
    }

   private:
    friend class Graph;     // Allow Graph to access Edge's private member data and functions.

    /** Constructs a valid edge within a specified graph
     * @param curr_graph  pointer to graph in which edge is added to
     * @param n1          id for node 1 of edge
     * @param n2          if for node 2 off edge
     *
     * @pre               graph != nullptr
     * @pre               no self edges. in other words, @nodeOne != @nodeTwo
     */
    Edge(const Graph* curr_graph, size_type n1, size_type n2) {
       graph = const_cast<Graph*>(curr_graph);
       nodeOne = n1;
       nodeTwo = n2;
       assert(graph != nullptr);
       assert(nodeOne != nodeTwo);
    }

    Graph* graph;
    size_type nodeOne;
    size_type nodeTwo;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return numEdges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const __attribute__ ((deprecated)){
    return Edge(this, i, true);
  }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (size_t i = 0; i < adj_list[a.index()].size(); i++) {
      if(adj_list[a.index()][i] == b.index()){
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

    // If edge exists, return existing edge: O(num_edges())
    for (size_type i = 0; i < adj_list[a.index()].size(); i++) {
      if(adj_list[a.index()][i] == b.index()){
        return Edge(this, a.index(), b.index());
      }
    }

    // If edge does not exist, add to vector and return new edge
    adj_list[a.index()].push_back(b.index()); // (edge id, b)
    adj_list[b.index()].push_back(a.index()); // (edge id, a)
    numEdges++;

    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    numEdges = 0;
    nodes.clear();
    adj_list.clear();
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
    NodeIterator() {}

    /** Dereferences the NodeIterator
     * @return      returns object of type Node, referenced to by the NodeIterator
     *
     * Complexity: O(1)
     */
    Node operator*() const{
      return Node(currGraph, currIndex);
    }

    /** Increments the NodeIterator to point to the next node in the graph.
     * @return    returns an incremented NodeIterator
     * @post      NodeIterator must point to the immediate next node in @a currGraph
     *
     * Complexity: O(1)
     */
    node_iterator& operator++(){
      currIndex = currIndex + 1;
      return *this;
    }

    /** Provides conditions for checking whether two NodeIterators are 'equal'.
     *
     * @param otherIt    The other NodeIterator which is being compared to 'this' NodeIterator
     * @return bool      returns true if the two iterators being compared are 'equal', or
     *                   if(@a currIndex == @a otherIt.currIndex && currGraph == otherIt.currGraph)
     *
     * Complexity: O(1)
     */
    bool operator==(const node_iterator& otherIt) const{
      return (currIndex == otherIt.currIndex && currGraph == otherIt.currGraph);
    }

   private:
    friend class Graph;

    /** Constructs a NodeIterator
     * @param graphIn  pointer to graph in which NodeIterator iterates over
     * @param vIn      node id from which to start iterating from
     *
     */
    NodeIterator(const Graph* graphIn, int vIn){
        currGraph = const_cast<Graph*>(graphIn);
        currIndex = vIn;
    }

    int currIndex;
    Graph* currGraph;
  };

  /** Returns a node iterator which allows the user to begin iterating over the nodes in the graph
   * @return  an iterator of type node_iterator
   *
   * Complexity: O(1)
   */
  node_iterator node_begin() const{
   return NodeIterator(this, 0);
  }

  /** Returns a node iterator which indicates the end of the iteration process
   * @return  an iterator of type node_iterator
   *
   * Complexity: O(1)
   */
  node_iterator node_end() const{
   return NodeIterator(this, this->num_nodes());
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
    IncidentIterator() {}

    /** Dereferences the IncidentIterator
     * @return      returns object of type edge, referenced to by the IncidentIterator
     *              (where node1 = @a n1 and node2 = @a adj_list[n1][n2].second)
     *
     * Complexity: O(1)
     */
    Edge operator*() const{
      return Edge(currGraph, n1, currGraph->adj_list[n1][n2]);
    }

    /** Increments the IncidentIterator to point to the next neighbor of node @n1.
     * @return    returns an incremented IncidentIterator
     * @post      must return a neighbor to @n1
     *
     * Complexity: O(1)
     */
    IncidentIterator& operator++(){
      n2++;
      return *this;
    }

    /** Provides conditions for checking whether two IncidentIterators are 'equal'.
     *
     * @param otherIt    The other IncidentIterator which is being compared to 'this' IncidentIterator
     * @return bool      returns true if the two iterators being compared are 'equal', else false
     *                   (n2 == otherIt.n2 && n1 == otherIt.n1 && currGraph == otherIt.currGraph)
     *
     * @pre              n1 != n2, or no self edges
     *
     * Complexity: O(1)
     */
    bool operator==(const IncidentIterator& otherIt) const{
      return (n2 == otherIt.n2 && n1 == otherIt.n1 && currGraph == otherIt.currGraph);
    }

   private:
    friend class Graph;

    /** Constructs an IncidentIterator
     * @param graphIn  pointer to graph in which IncidentIterator iterates over
     * @param one      source node id (@a n1)
     * @param two      index referring to index of adjacent nodes to @n1
     *
     */
    IncidentIterator(const Graph* graphIn, int one, int two) {
      n1 = one;
      n2 = two;
      currGraph = const_cast<Graph*>(graphIn);
    }
    int n1;
    int n2;
    Graph* currGraph;
  };

  //
  // Edge Iterator
  //

  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    // Construct an invalid EdgeIterator.
    EdgeIterator(){}

    /** Dereferences the EdgeIterator
     * @return      returns object of type edge, referenced to by the EdgeIterator
     *              (where node1 = @a currNode and node2 = @a adj_list[currNode][currNeigh].second)
     *
     * Complexity: O(1)
     */
    Edge operator*() const{
      return Edge(currGraph, currNode, currGraph->adj_list[currNode][currNeigh]);
    }

    /** Increments the EdgeIterator to point to the next edge in @a currGraph.
     * @pre       an adjacency list is used to represent the edges and each edge
     *            is stored twice, such that both (a,b) and (b,a) can be found in
     *            the adjacency list.
     * @return    returns an incremented EdgeIterator to a unique edge that has not
     *            yet been visited. If n1 < n2, then the edge is returned.
     *
     * EdgeIterator iterates over the internal data structure, in this case,
     * an adjacency list (or a vector of vectors) until all unique edges have been
     * visited.
     *
     * Complexity: O(m)
     */
    EdgeIterator& operator++(){
      ++currNeigh;
      if(currNeigh >= currGraph->adj_list[currNode].size()){
        ++currNode;
        currNeigh = 0;
      }
      while(currNode < currGraph->adj_list.size() && (currNode > currGraph->adj_list[currNode][currNeigh])){
        if(currNeigh < currGraph->adj_list[currNode].size()-1){
          ++currNeigh;
        }
        else{
          ++currNode;
          currNeigh = 0;
        }
      }
      return *this;
    }

    /** Provides conditions for checking whether two EdgeIterators are 'equal'.
     *
     * @param otherIt    The other EdgeIterator which is being compared to 'this' EdgeIterator
     * @return bool      returns true if the two iterators being compared are 'equal', else false
     *                   ((currGraph == otherIt.currGraph) && (currNode == otherIt.currNode)
     *                        && (currNeigh == otherIt.currNeigh))
     *
     * Complexity: O(1)
     */
    bool operator==(const EdgeIterator& otherIt) const{
      return ((currGraph == otherIt.currGraph) && (currNode == otherIt.currNode) && (currNeigh == otherIt.currNeigh));
    }

   private:
    friend class Graph;
    Graph* currGraph;
    size_type currNode;
    size_type currNeigh;

    /** Constructs an EdgeIterator
     * @param g       pointer to graph in which EdgeIterator iterates over
     * @param node    source node id (@a currNode)
     * @param neigh   index referring to index of adjacent node to @a currNode
     *
     */
    EdgeIterator(const Graph* g, size_type node, size_type neigh){
      currGraph = const_cast<Graph*>(g);
      currNode = node;
      currNeigh = neigh;
    }

  };

  /** Returns an edge iterator which allows the user to begin iterating over the edges in the graph
   * @return  an iterator of type edge_iterator
   *
   * Creates an EdgeIterator which sets currNode = 0 and currNeigh = 0,
   * referring to the first edge in the adjacency list.
   *
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0, 0);
  }

  /** Returns an edge iterator which indicates the end of the iteration process
   * @return  an iterator of type edge_iterator
   * @post    the iterator refers to the first invalid point in the iteration process
   *
   * Creates an EdgeIterator in which currNode is set to (index of the last node + 1)
   * and currNeigh is set to 0.
   * Complexity: O(1)
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, adj_list.size(), 0);
  }

  /** Outputs the graph for testing purposes
   * @return  prints the adjacency list to std::cout
   *
   */
  void print_graph(){
    for(size_type i = 0; i < adj_list.size(); ++i){
      std::cout << i << ':';
      for(size_type j = 0; j < adj_list[i].size(); ++j){
        std::string x;
        std::cout << adj_list[i][j] << ' ';
      }
      std::cout << std::endl;
    }
  }

 private:
  size_type numEdges; // stores number of edges
  std::vector<std::pair<Point,node_value_type>> nodes; // stores nodes of graph
  std::vector<std::vector<size_type>> adj_list; // adjacency list that stores edges

};

#endif // CME212_GRAPH_HPP
