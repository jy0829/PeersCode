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
  struct internal_node;
  struct internal_edge;

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
//Do I need to specify the initial parameters?
  Graph()
  : internal_nodes(), internal_edges(), n_nodes(0), n_edges(0), next_node_uid(0), next_edge_uid(0) {
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
//Should do this later
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch().index;
    }

    size_type uid() const {
      // HW0: YOUR CODE HERE
      return fetch().uid;
    }

    // Return the pointer to the graph
    const Graph* graph_pointer() const {
      // HW0: YOUR CODE HERE
      return graph_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((n.index()==index()) && (n.graph_pointer()==graph_pointer())){
        return true;
      } else {
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
      if (index()<n.index()){
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // A pointer to the graph
    Graph* graph_;
    // A unique identification number
    size_type uid_;

    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    // Fetches the internal node corresponding to the node
    internal_node& fetch() const {
      for (size_type i = 0; i < graph_->num_nodes(); ++i)
        if (graph_->internal_nodes[i].uid == uid_)
          return graph_->internal_nodes[i];
      assert(false);
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    internal_node new_internal_node = internal_node(position, next_node_uid, n_nodes);
    internal_nodes.push_back(new_internal_node);
    ++next_node_uid;
    ++n_nodes;
    return Node(this, next_node_uid-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_pointer() == this){
      for (unsigned int i = 0; i < n_nodes; i++){
        if (internal_nodes[i].uid == n.uid()){
          return true;
        }
      }
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
    return Node(this, internal_nodes[i].uid);        // Invalid node
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

    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch().index;
    }

    size_type uid() const {
      // HW0: YOUR CODE HERE
      return fetch().uid;
    }

    // Return the pointer to the graph
    Graph* graph_pointer() {
      // HW0: YOUR CODE HERE
      return graph_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, (fetch().pointer1)->uid);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, fetch().pointer2->uid);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==( Edge& e) const {
      if (e.index()==index() && e.graph_pointer()==graph_){
        return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (index()<e.index()){
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    size_type uid_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    // Fetches the internal edge corresponding to the edge
    internal_edge& fetch() const {
      for (size_type i = 0; i < graph_->num_edges(); ++i)
        if (graph_->internal_edges[i].uid == uid_)
          return (graph_->internal_edges)[i];
      assert(false);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return n_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, internal_edges[i].uid);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (unsigned int i = 0; i < n_edges; ++i){
      const internal_edge* edge = &internal_edges[i];
      size_type uid1 = a.uid();
      size_type uid2 = b.uid();
      if (edge->pointer1->uid == uid1 && edge->pointer2->uid == uid2){
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
    internal_edges.push_back(internal_edge(n_edges, next_edge_uid, &(a.fetch()), &(b.fetch())));
    ++next_edge_uid;
    ++n_edges;
    return Edge(this, next_edge_uid-1);        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    internal_nodes.clear();
    internal_edges.clear();
    n_nodes = 0;
    n_edges = 0;
  }

 private:
   struct internal_node {
     internal_node( Point point_, size_type index_, size_type uid_){
       point = point_;
       index = index_;
       uid = uid_;
     }
    Point point;         // The position of the node
     size_type index;     // The index of the node
     size_type uid;      // The unique identifcation for the node
   };

   struct internal_edge {
     internal_edge(size_type index_, size_type uid_, internal_node* pointer1_, internal_node* pointer2_){
       index = index_;
       uid = uid_;
       pointer1 = pointer1_;
       pointer2 = pointer2_;
     }
     size_type uid;     //The unique identification for the edge
     size_type index;   //The index of the edge
     internal_node* pointer1;      // Pointer to node1
     internal_node* pointer2;  // Pointer to node2
   };


  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Define a vector of nodes and edges, as well as the number of nodes, edges
  // and the next unique identifier
  std::vector<internal_node> internal_nodes;
  std::vector<internal_edge> internal_edges;
  size_type n_nodes;
  size_type n_edges;
  size_type next_node_uid;
  size_type next_edge_uid;


};

#endif // CME212_GRAPH_HPP
