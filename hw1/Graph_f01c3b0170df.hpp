#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 * @author Greg DePaul
 */

#include <algorithm>
#include <vector>
#include <list>
#include <cassert>
#include <bitset>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

template <typename V>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 public:

  /**********************************************************************
  / PUBLIC TYPE DEFINITIONS
   **********************************************************************/

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;

  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Synonym for template V */
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

  /**********************************************************************
  / CONSTRUCTORS AND DESTRUCTOR
   **********************************************************************/

  /** Construct an empty graph. */
  Graph() {

  }

  /** Default destructor */
  ~Graph() = default;

  /**********************************************************************
  / NODE (INNER CLASS)
   **********************************************************************/

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */

  //NOTE: Size: 16 Bytes
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
    Node() : myIndex(size_type(-1)) {
      myParent = NULL; 
    }

    /** Construct an valid Node. */
    Node(size_type newIndex, Graph<node_value_type>* parent) {   
      myParent = parent;
      myIndex = newIndex;	
    }

    /** Return this node's position. */
    const Point& position() const {
      return *(myParent->myPositions.at(myIndex));
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return myIndex;
    }

    /** Return this node's degree, a number in the range [0, graph_size). */
    size_type degree() const {
      return (myParent->myConnections.at(myIndex))->size();
    }

    /** Return an iterator that points to the first incident edge of this node. */
    IncidentIterator edge_begin() const {
      return IncidentIterator(myParent, myIndex);
    }
    
    /** Return an iterator that points to the end incident edge of this node. */
    IncidentIterator edge_end() const {
      return IncidentIterator(myParent, myIndex, myParent->myConnections.at(myIndex)->end());
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    node_value_type& value() {
      return *(myParent->myNodeValues.at(myIndex));
    }

    /** Return this node's value of type V */
    const node_value_type& value() const {
      return *(myParent->myNodeValues.at(myIndex));
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return n.index() == index() && n.position() == position() && n.myParent == myParent;
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
      return n.myParent == myParent && n.index() < index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph<node_value_type>* myParent; 
    size_type myIndex;  
 
  };

  /**********************************************************************
  / GRAPH CLASS METHODS ON NODES
   **********************************************************************/

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return myPositions.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new nodes's stored value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& newValue = node_value_type()) {
    Node newNode = add_node_helper(position);
    *(myNodeValues.at(newNode.index())) = newValue;
    return newNode;
  }

  /** Add node helper function. Maintains the functionality of 
   * the previous version's add_node function.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node_helper(const Point& position) {
    myPositions.push_back(new Point(position));
    myConnections.push_back(new std::list<size_type>());
    myNodeValues.push_back(new node_value_type());
    return node(num_nodes()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.index() < num_nodes() && n == node(n.index());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(i, const_cast<Graph<node_value_type>*>(this));
  }

  /**********************************************************************
  / EDGE (INNER CLASS)
   **********************************************************************/

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */

  //NOTE: Size: 24 Bytes
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      myNode1 = (size_type)-1;
      myNode2 = (size_type)-1;
      myParent = NULL;
    }

    /** Construct an valid Edge. */
    Edge(size_type nodeAIndex, size_type nodeBIndex, Graph<node_value_type>* parent) {
	myNode1 = nodeAIndex;
	myNode2 = nodeBIndex; 
	myParent = parent; 
    }

    /** Return a node of this Edge */
    Node node1() const {
      return myParent->node(myNode1);     
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return myParent->node(myNode2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (myParent == e.myParent) && ((node1() == e.node1() && node2() == e.node2()) || (node1() == e.node2() && node2() == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return myParent == e.myParent && node1().index() < e.node1().index() && node2().index() < e.node2().index();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    size_type myNode1;
    size_type myNode2;
    Graph<node_value_type>* myParent;

  };

  /**********************************************************************
  / GRAPH CLASS METHOD ON EDGES
   **********************************************************************/

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  // NOTE: Complexity: O(1)
  size_type num_edges() const {
    return myEdges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  // NOTE: Complexity: O(1)
  Edge edge(size_type i) const {
    return Edge(myEdges.at(i)[0], myEdges.at(i)[1], const_cast<Graph<node_value_type>*>(this));     
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  // NOTE: Complexity: O(num_nodes() + num_edges())
  bool has_edge(const Node& a, const Node& b) const {
    std::list<size_type>::iterator currNode; 
    for(currNode = myConnections.at(a.index())->begin(); currNode != myConnections.at(a.index())->end(); currNode++) {
      if(a.myParent->node(*currNode) == b) {
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
 
  // NOTE: O(1) amortized time [ O(log(num_nodes)) depending on the std::vector resizing factor]
  Edge add_edge(const Node& a, const Node& b) {
    Edge newEdge = Edge(a.index(), b.index(), this); 
    if(!(a==b) && !has_edge(a, b)) {
      myConnections.at(a.index())->push_back(b.index());
      myConnections.at(b.index())->push_back(a.index());
      size_type* newArray = (size_type *)malloc(2*sizeof(size_type));
      newArray[0] = a.index();
      newArray[1] = b.index();
      myEdges.push_back(newArray);
      return newEdge;
    } 
    else {
      return newEdge;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    
    // Clear myConnections
    for(std::vector<std::list<size_type>*>::iterator it = myConnections.begin(); it != myConnections.end(); it++) {
      delete(*it);
    }

    // Clear myPositions
    for(std::vector<Point*>::iterator it = myPositions.begin(); it != myPositions.end(); it++) {
      delete(*it);
    }

    // Clear myNodeValues
    for(typename std::vector<node_value_type*>::iterator it = myNodeValues.begin(); it != myNodeValues.end(); it++) {
      delete(*it);
    }

    // Clear myEdges
    for(std::vector<size_type*>::iterator it = myEdges.begin(); it != myEdges.end(); it++) {
      delete(*it);
    }
  }

  /** Print out the statuses of the graph's nodes */
  void printNodes() {
    std::cout << "Index" << ":\t" << "Position" << "\t" << "Value" << "\t" << "Degree" << std::endl;
    for(auto ni = node_begin(); ni != node_end(); ++ni) {
        auto currNode = *ni;
	std::cout << currNode.index() << ":\t" << currNode.position() << "\t\t" << currNode.value() << "\t" << currNode.degree() << std::endl;
    }
  }

  /**********************************************************************
  / NODE ITERATOR (INNER CLASS)
   **********************************************************************/

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
      myIndex = size_type(-1);
      myParent = NULL;
    }

    /** Construct a valid NodeIterator. */
    NodeIterator(Graph<node_value_type>* newParent, size_type newIndex = (size_type)0) {
      myIndex = newIndex;
      myParent = newParent;
    }
    
    // Return the node this iterator points to
    Node operator*() const {
      return myParent->node(myIndex);
    }
    
    // Increment the iterator to the next node in graph
    NodeIterator& operator++() {
      myIndex++;
      return *this;
    }
    
    // Check equality of two iterators (specifically, do they point to the same node)
    bool operator==(const NodeIterator& other) const {
      return myParent == other.myParent && myIndex == other.myIndex;
    }

   private:
    friend class Graph;
    size_type myIndex;
    Graph<node_value_type>* myParent;
  };

  /**********************************************************************
  / GRAPH CLASS METHODS ON NODE ITERATOR
   **********************************************************************/

  // Return an iterator that points to the begining of a container of nodes
  NodeIterator node_begin() const {
    return NodeIterator(const_cast<Graph<node_value_type>*>(this));
  }
  
  // Return an iterator that points to the end of a container of nodes
  NodeIterator node_end() const {
    return NodeIterator(const_cast<Graph<node_value_type>*>(this), num_nodes());
  }

  /**********************************************************************
  / INCIDENT ITERATOR (INNER CLASS)
   **********************************************************************/

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
      myNodeIndex = (size_type)(-1);
      myParent = NULL;
      currEdge = std::list<size_type>::iterator();
    }

    /** Construct a valid IncidentIterator that points to the edge
       specified by the user. */
    IncidentIterator(Graph<node_value_type>* newParent, size_type newIndex, std::list<size_type>::iterator newEdge) {
      myNodeIndex = newIndex;
      myParent = newParent;
      currEdge = newEdge;
    }

    /** Construct a valid IncidentIterator that points to the beginning
       of a list of edges incident to a node. Can't abuse default operators
       because the default is dependent of the other parameters. */
    IncidentIterator(Graph<node_value_type>* newParent, size_type newIndex) {
      myNodeIndex = newIndex;
      myParent = newParent;
      currEdge = newParent->myConnections.at(newIndex)->begin();
    }

    // Return the edge this iterator points to
    Edge operator*() const {
      return Edge(myNodeIndex, (*currEdge), myParent);
    }
    
    // Increment this iterator to the next edge in the node's adjacency list
    IncidentIterator& operator++() {
      currEdge++;
      return *this;
    }
    
    // Check to see if this iterator points to the same edge as the other iterator
    bool operator==(const IncidentIterator& other) const {
      return myParent == other.myParent && currEdge == other.currEdge;
    }

   private:
    friend class Graph;
    std::list<size_type>::iterator currEdge;
    Graph<node_value_type>* myParent;
    size_type myNodeIndex;
  };

  /**********************************************************************
  / EDGE ITERATOR (INNER CLASS)
   **********************************************************************/

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
	myParent = NULL;
	myEdgeIndex = (size_type)(-1);
    }

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(size_type startIndex, Graph<node_value_type>* parent) {
	myParent = parent;
	myEdgeIndex = startIndex;
    }

    // Return the edge this iterator points to
    Edge operator*() const {
	return myParent->edge(myEdgeIndex);
    }
    
    // Increment this iterator to the next edge in the global edge container
    EdgeIterator& operator++() {
	myEdgeIndex++;
	return *this;
    }
    
    // Check to see if this iterator points to the same edge as the other iterator
    bool operator==(const EdgeIterator& other) const {
	return myEdgeIndex == other.myEdgeIndex;
    }

   private:
    friend class Graph;
    size_type myEdgeIndex;
    Graph<node_value_type>* myParent; 
  };

  /**********************************************************************
  / GRAPH CLASS METHODS ON EDGE ITERATOR
   **********************************************************************/

  // Return an iterator that points to the first edge in the graph's edge list
  edge_iterator edge_begin() const {
	return EdgeIterator(0, const_cast<Graph<node_value_type>*>(this));
  }

  // Return an iterator that points to the last edge in the graph's edge list
  edge_iterator edge_end() const {
	return EdgeIterator(num_edges(), const_cast<Graph<node_value_type>*>(this));
  }

 private:

  /**********************************************************************
  / DATASTRUCTURES
   **********************************************************************/
  std::vector<size_type*> myEdges;
  std::vector<std::list<size_type>*> myConnections;
  std::vector<Point*> myPositions;
  std::vector<node_value_type*> myNodeValues;

};

#endif // CME212_GRAPH_HPP
