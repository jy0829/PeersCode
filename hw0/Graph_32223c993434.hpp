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
  struct internal_node;
  struct internal_edge;
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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
        : nodes_(), edges_(), size_node_(0), next_node_uid_(0),size_edge_(0), next_edge_uid_(0)  {
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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_ -> nodes_[node_uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_uid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n; // Quiet compiler warning
      return (n.graph_ == graph_ and n.index() == node_uid_);           
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
      if (n.graph_ == graph_ and n.index() == node_uid_)
        return false;
      else 
        return node_uid_ < n.index();

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_; //Pointer back to the Graph container
    size_type node_uid_;// This element's unique identification number 

    Node(const Graph* graph, size_type node_uid)
        : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid){
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
    return size_node_;
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
    (void) position;      // Quiet compiler warning
    internal_node new_node; // Defining a new node
    new_node.node_uid = next_node_uid_; //Assigning the right attributes
    new_node.position = position; 
    nodes_.push_back(new_node);// Adding this node to our vector
    ++size_node_;//Incrementing the size and the next node's index
    ++next_node_uid_;
    return Node(this, next_node_uid_ -1);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    if(size_node_ == 0) 
      return false;// Graph empty 
    else 
      return n.graph_ == this; //Checking if the graphs attached to the two nodes are the same
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
    return Node(this,i); //Return the node from the graph that has index i (see constructor) 
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
      return graph_ -> edges_[edge_uid_].node_1;  
    // HW0: YOUR CODE HERE
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_ -> edges_[edge_uid_].node_2;      
// HW0: YOUR CODE HERE
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if (! (e.graph_ == graph_)) //Check if the two edges belong to the same graph
        return false;
     else
      return (e.graph_ == graph_ and e.edge_uid_ == edge_uid_);

    }
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return edge_uid_ < e.edge_uid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type edge_uid_;

    Edge(const Graph* graph, size_type edge_uid)
      : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid){
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
    return size_edge_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
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
    (void) a; (void) b;   // Quiet compiler warning
    if (edge_helper_.count(a) > 0){
      if(edge_helper_.at(a).count(b) > 0)
        return true;
    }
    if (edge_helper_.count(b) > 0){
      if(edge_helper_.at(b).count(a) > 0)
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
    if(has_edge(a,b)){
      return Edge(this, edge_helper_[a][b]);
    }
    else{ 
      internal_edge new_edge;
      new_edge.node_1 = a;
      new_edge.node_2 = b; 
      new_edge.edge_uid = next_edge_uid_;
      edges_.push_back(new_edge);
      edge_helper_[a][b] = next_edge_uid_;
      edge_helper_[b][a] = next_edge_uid_;
      ++size_edge_;
      ++next_edge_uid_;
      return Edge(this, next_edge_uid_ -1);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    edge_helper_.clear();
    size_edge_ = 0; 
    size_node_ = 0;
    next_node_uid_ = 0;
    next_edge_uid_ = 0;
    // HW0: YOUR CODE HERE
  }

 private:
  struct internal_node {
    Point position; //position (coordinates) of the node
    size_type node_uid; //id of the node: the unique identification for an element
 };

  struct internal_edge {
    Node node_1; // One of the two node that the edge 
    Node node_2;
    size_type edge_uid; //The unique identification for an element
};  

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
//If node a and node b are connected by edge uid, then edge_helper_[a][b] = uid
// and edge_helper_[b][a] = uid 
  std::map<Node, std::map<Node, size_type>> edge_helper_;
  size_type size_node_;
  size_type size_edge_;
  size_type next_node_uid_;
  size_type next_edge_uid_;



  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
