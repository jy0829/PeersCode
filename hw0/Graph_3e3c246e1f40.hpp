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

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


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

  /** Predeclaration of Edge type. */
  class Edge;

  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /**********************************************************************
  / CONSTRUCTORS AND DESTRUCTOR
   **********************************************************************/

  /** Construct an empty graph. */
  Graph() {
    mySize = 0;
    Graph::myPositions = std::vector<Point*>();
    edgeSize = 0;
    Graph::myConnections = std::vector<std::list<int>*>(); 
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
  class Node {
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

    /** Return this node's position. */
    const Point& position() const {
      return *(myParent->myPositions.at(myIndex));
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return myIndex;
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
      return n.position() != position() && n.index() < index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    const Graph* myParent; 
    size_type myIndex; 

    Node(size_type newIndex, const Graph* parent) {
      myParent = parent;
      myIndex = newIndex;	
    }
 
  };

  /**********************************************************************
  / GRAPH CLASS METHODS ON NODES
   **********************************************************************/

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return mySize;
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
  Node add_node(const Point& position) {
    Graph::myPositions.push_back(new Point(position));
    Graph::myConnections.push_back(new std::list<int>());
    return Node(mySize++, this);
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
    return Node(i, this);
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      myNode1 = (size_type)-1;
      myNode2 = (size_type)-1;
      myParent = NULL;
      myIndex = (size_type)-1; 
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
      return myParent == e.myParent && node1() == e.node1() && node2() == e.node2();
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return node1().index() < e.node1().index() && node2().index() < e.node2().index();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    size_type myNode1;
    size_type myNode2;
    Graph* myParent;
    size_type myIndex;
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
    return edgeSize;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  // NOTE: Complexity: O(1)
  Edge edge(size_type i) const {
    return myEdges.at(i);     
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  // NOTE: Complexity: O(num_nodes() + num_edges())
  bool has_edge(const Node& a, const Node& b) const {
    std::list<int>::iterator currNode; 
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
    Edge newEdge = Edge();
    newEdge.myNode1 = a.index();
    newEdge.myNode2 = b.index();
    newEdge.myParent = this;
    newEdge.myIndex = edgeSize; 
    if(!(a==b) && !has_edge(a, b)) {
      myConnections.at(a.index())->push_back(b.index());
      myConnections.at(b.index())->push_back(a.index());
      edgeSize++; 
      myEdges.push_back(newEdge);
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
    for(std::vector<std::list<int>*>::iterator it = myConnections.begin(); it != myConnections.end(); it++) {
      delete(*it);
    }
    myEdges.clear(); 
    for(std::vector<Point*>::iterator it = myPositions.begin(); it != myPositions.end(); it++) {
      delete(*it);
    }
  }

 private:

  /**********************************************************************
  / DATASTRUCTURES
   **********************************************************************/
  size_type mySize;
  size_type edgeSize;
  std::vector<Edge> myEdges;
  std::vector<std::list<int>*> myConnections;
  std::vector<Point*> myPositions;

};

#endif // CME212_GRAPH_HPP
