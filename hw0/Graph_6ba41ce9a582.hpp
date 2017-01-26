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
class Graph {

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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;


 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Declaring private attributes points, edge class and vector

 std::vector<Point> Points;


 struct edge_items
 {
     size_type NodeId1;
     size_type NodeId2;
 
     edge_items(size_type id1,size_type id2)

     {
          NodeId1=id1;
          NodeId2=id2;
     }
 };

 std::vector<edge_items> Edges;

 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //


  /** Construct an empty graph. */

  Graph() : Points(),Edges() {}
 

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

      
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
        return GraphPointer->Points[NodeId];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(GraphPointer->Points.size()>NodeId);
      return NodeId;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return GraphPointer == n.GraphPointer && index() == n.index();
      }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const 
    {

        if (GraphPointer == n.GraphPointer) 
            return NodeId < n.NodeId;
        else 
            return GraphPointer < n.GraphPointer;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Constructor that initalizes Node object

    Node(const Graph* currentgraph, size_type index)
    {
        NodeId = index;
        GraphPointer = currentgraph;
    }

    // Declaring private variables

    size_type NodeId;
    const Graph* GraphPointer;
    

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const 
  {
      return Points.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const 
  {
      return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) 
  {

      // Pushing back position to points vector

      Points.push_back(position); // Add x,y,z to point vector
      Node NodeObject(this, size()-1);
      return NodeObject;

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */

  bool has_node(const Node& n) const 
  {
      size_type numofnodes = size();          
      
      if (numofnodes < n.index())
      {
        return false;
      }
    
      else
      {
          return true;
      }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */

  Node node(size_type i) const 
  {
      assert(Points.size()>i); // Asserting numnodes > i
      return Node(this,i);       
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
  class Edge {
   public:

    /** Construct an invalid Edge. */
    Edge() {
    }

    Edge(const Graph* currentgraph, size_type id1, size_type id2)
    {
        NodeId1=id1;
        NodeId2=id2;
        GraphPointer=currentgraph;
    }

    /** Return a node of this Edge */
    Node node1() const 
    {
      return Node(GraphPointer,NodeId1);      
    }

    /** Return the other node of this Edge */
    Node node2() const 
    {
      return Node(GraphPointer,NodeId2);   
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const 
    {
   
    // Checking for equality
 
    return (GraphPointer==e.GraphPointer &&
            std::min(NodeId1,NodeId2)==std::min(e.NodeId1,e.NodeId2) &&
            std::max(NodeId1,NodeId2)==std::max(e.NodeId1,e.NodeId2));
     
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const 
    {

      return (GraphPointer==e.GraphPointer && 
              std::min(NodeId1,NodeId2)<=std::min(e.NodeId1,e.NodeId2) &&
              std::max(NodeId1,NodeId2)<=std::max(e.NodeId1,e.NodeId2));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    const Graph* GraphPointer;
    size_type NodeId1;
    size_type NodeId2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  size_type num_edges() const 
  {
    return Edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const 
  {

    assert(Edges.size()>i); //Asseting that i < number of edges
    Edge EdgeObject(this,Edges[i].NodeId1,Edges[i].NodeId2);
    return EdgeObject;    

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const 
  {

    assert(a.GraphPointer== this and b.GraphPointer == this); // asserting nodes in graph
    assert(a.NodeId < size() && b.NodeId < size()); // asserting nodes are valid
    

    for (unsigned int ind = 0; ind < Edges.size(); ind++) // checking for edge
    {
        if (std::min(a.NodeId,b.NodeId)==std::min(Edges[ind].NodeId1,Edges[ind].NodeId2) &&
           std::max(a.NodeId,b.NodeId)==std::max(Edges[ind].NodeId1,Edges[ind].NodeId2))
           {
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
  Edge add_edge(const Node& a, const Node& b) 
  {

    assert(this==a.GraphPointer && this == b.GraphPointer && // asserting preconditions
    a.NodeId < Points.size() && b.NodeId < Points.size() && a.NodeId!=b.NodeId);
    
    if (has_edge(a,b))
    {
        Edge EdgeObject(this,a.NodeId,b.NodeId);
        return EdgeObject;
    }  
   
    
        
    edge_items edgeData(a.NodeId,b.NodeId);
    Edges.push_back(edgeData);
    Edge EdgeObject(this,a.NodeId,b.NodeId);
    return EdgeObject;
    

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Points.clear();
    Edges.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP



