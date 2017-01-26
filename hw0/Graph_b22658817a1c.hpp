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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
  //It has private vector v_set, length v_num, vector e_set, length e_num. It represents the mathematical graph.

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

     //Node has pointer to graph @a graph, and graph index @a idx, which is read by the index() function.
    Node() {
      // HW0: YOUR CODE HERE
      /**Creates an invalid node, which may have no connection to any graph, 
      or have no valid point assigned to it.
      */ 

      
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      /**Looks into its parent graph vertex set @a v_set, finds the point corresponding to the node's index, returns point.
      *@return Point @a p.
      *@post Point @a p == graph.node(@a this->idx)
      */
      return this->graph->get_pt(this->index());
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      /*Returns the index of the node in the graph.
      *@return size_type @a idx
      *@post idx belongs to interval [0,Node->graph.graph_size).
      */
      return idx;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      /**Checks if the node equals node @a n by matching their graphs and indices.
      *@param n: Node we compare with.
      *@return bool @a b.
      *@post b=1 if (@a this->graph == @a n.graph and @a this->index() == @a n.index()), b=0 otherwise.
      */
      if ((this->graph == n.graph) && (this->index() == n.index())) {
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
      /**Checks if the node is less than node @a n by comparing their indices.
      *@pre @a this->graph == @a n.graph.
      *@param n: Node we compare with.
      *@return bool @a b.
      *@post b=1 if (@a this->index() < @a n.index()), b=0 otherwise.
      */
      if (this->index() < n.index()) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    const Graph* graph;
    size_type idx;
    Node(const Graph* g,size_type i) : graph(g), idx(i) {
    /** Creates a valid node with pointer @a graph pointing to @a g and @a this->index() equal to @a i.
    *@param g: pointer to parent graph.
    *@param i: index of the node corresponding to the point in graph
    *@return None. This is a constructor.
    *@post @a this->graph == @a g
    *@post @a this->index() == @a i
    */
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    /*Returns the number of points in the graph vertex set.
    *@return @a l.
    *@post @a l == @a this->length.
    */
    return this->v_num;
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
    // HW0: YOUR CODE HERE
    this->v_set.push_back(position);
    Node node(this,this->num_nodes());
    this->v_num+=1;
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      /*Look into @a graph.e_set and collect first node corresponding to @a idx.*/
      return graph->node(graph->e_set[idx].first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph->node(graph->e_set[idx].second);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      /*Checks if two edges have the same graph and then if the node indices match.
      *@param e, the edge with which we compare @a this.
      *@return bool @a b, which signals equality.
      *@post @a b == 1 if @a this->graph == @a e.graph and {@a node1(),@a node2()} == {@a e.node1(),@a e.node2()}, 
      *@a b = 0 otherwise.
      */
      if (graph == e.graph) {
        if (graph->e_set[idx].first+graph->e_set[idx].second == e.graph->e_set[idx].first+e.graph->e_set[idx].second) {
          if (graph->e_set[idx].first*graph->e_set[idx].second == e.graph->e_set[idx].first*e.graph->e_set[idx].second) {
            return true;
          }
        }
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      /*Given two edges of the same graph, returns the value of a mathematically imposed ordering based on its indices.
      *@param e, the edge with which we compare @a this.
      *@pre @a this->graph == @a e->graph.
      *@return bool @a b, which signals comparison.
      *@post < satisfies trichotomy as defined for any mathematical order relation. 
      */
      if (graph->e_set[idx].first+graph->e_set[idx].second < e.graph->e_set[idx].first+e.graph->e_set[idx].second) {
        return true;
      }
      if (graph->e_set[idx].first+graph->e_set[idx].second == e.graph->e_set[idx].first+e.graph->e_set[idx].second) {
        if (graph->e_set[idx].first*graph->e_set[idx].second < e.graph->e_set[idx].first*e.graph->e_set[idx].second) {
          return true;
        }
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const Graph* graph;
    size_type idx;
    Edge(const Graph* g, size_type i) : graph(g), idx(i) {}; 
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return e_num;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type i1 = a.index();
    size_type i2 = b.index();
    size_type j1,j2;
    for (auto itr = e_set.begin();itr != e_set.end();itr++) {
      j1 = itr->first;
      j2 = itr->second;
      if (i1*i2 == j1*j2 && i1+i2 == j1+j2) {
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
    size_type i1 = a.index();
    size_type i2 = b.index();
    size_type j1,j2;
    size_type idx = 0;
    for (auto itr = e_set.begin();itr != e_set.end();itr++) {
      j1 = itr->first;
      j2 = itr->second;
      if (i1*i2 == j1*j2 && i1+i2 == j1+j2) {
        return Edge(this,idx);
      }
      idx+=1;
    }
    std::pair<size_type,size_type> e_pair;
    e_pair = std::make_pair(i1,i2);
    e_set.push_back(e_pair);
    e_num+=1;
    return Edge(this,idx);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    v_set.clear();
    e_set.clear();
    v_num = 0;
    e_num = 0;
  }
  
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  /*@a v_set is the vertex set of the graph.*/
  std::vector<Point> v_set;
  /*@a v_num is the number of vertices in the graph.*/
  size_type v_num = 0;


  /*@a v_set is the edge set of the graph. Edges stored as two node indices.*/
  std::vector<std::pair<size_type,size_type>> e_set;
  /*@a e_num is the number of vertices in the graph.*/
  size_type e_num = 0;

  const Point& get_pt(size_type i) const {
  /**Returns the address to the point at index @a i. 
  This is a layer of abstraction to access the Graph data structures.
  *@param i is the index of the point.
  *@return @a p, a reference to the point in the graph vertex set.
  *@post @a p == @a this->v_set[i].
  */
  return this->v_set[i];
  }
  
};

#endif // CME212_GRAPH_HPP
