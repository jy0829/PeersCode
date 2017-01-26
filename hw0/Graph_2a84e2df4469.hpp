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
    Node() {
      // HW0: YOUR CODE HERE

    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE

      return graph_->pts[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE

      return uid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //compare 2 graphs - same uid and node on same graph

      return ((n.graph_ ==  graph_) && (n.uid_ == uid_));
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

      //COMPARE IF NOT THE SAME GRAPH - ADD IT
      //There are 2 cases - if on different graphs and if on same graph
      if (n.graph_ != graph_) return (n.graph_ < graph_);
      if (n.graph_ == graph_) return (n.uid_ < uid_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // A node is defined by its pointer to the graph, and its id
    Graph* graph_;
      // This element's unique identification number
      size_type uid_;
      /** Private Constructor */
      Node(const Graph* graph, size_type uid)
              : graph_(const_cast<Graph*>(graph)), uid_(uid) {
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return pts.size();
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

    //when we add a node, we need to add to the - vector of nodes, but also
    // to the edges via the adjoined List
    pts.push_back(position);
    adjList.push_back(std::vector<size_type>());
    return Node(this, (size()-1));
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    return (n.graph_ == this && (n.uid_ >=0 && n.uid_<=size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

    assert(i < size());
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(uid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //written code though no indication
      //conditions: belong to same graph and see if the nodes id's match (check both directions)

      return (((uid1_== e.uid1_) && (uid2_ == e.uid2_) && (graph_ == e.graph_)) ||
              ((uid1_== e.uid2_) && (uid2_ == e.uid1_) && (graph_ == e.graph_)));

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //no need to write code here - by the hw's indications; wrote it anyways

      //IDEA - look, in this order, at both graphs and nodes
      //If different graphs, look directly at the graphs id's
      // If same graphs, check first by minumum, and if same minimum, by maximum

      //If different graph, simiply look at their graphs ids
      if (graph_ != e.graph_) {return graph_ < e.graph_;}

      // if on same graph, check the min and max between 2 nodes
      if (graph_ == e.graph_)
         {
            if (std::min(e.uid1_, e.uid2_) != std::min(uid1_,uid2_))
                {return std::min(e.uid1_, e.uid2_) < std::min(uid1_,uid2_);}
            else
                {return std::max(e.uid1_, e.uid2_) < std::max(uid1_,uid2_);}
         }

    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // An Edge is defined by its pointer to a gprah, and its 2 nodes(id's of nodes)
    Graph* graph_;
    // This element's unique identification number
    size_type uid1_;
    size_type uid2_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type uid1, size_type uid2)
           : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
      }

      //
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE

    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    auto id1 = edges[i].first;
    auto id2 = edges[i].second;
    return Edge(this, id1, id2);

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //chevk via the adjointed list, where the adjointed list is like a vector
    // of vector; to remember later in terms of x[i][j], in terms of the search
    for (int i=0; i<adjList[a.uid_].size(); i++) {
        if (adjList[a.uid_][i]==b.uid_) return true;
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
    //Case 1: if the add we want to add already exists
    if (has_edge(a, b)) {return Edge(this, a.uid_, b.uid_);}
    // Case 2: if the edge does not exist
    else {
          // Add to the adjointed list
          adjList[a.uid_].push_back(b.uid_);
          adjList[b.uid_].push_back(a.uid_);
          // Add also to the pairs, in a sorted way
          if (a.uid_ < b.uid_) edges.push_back(std::make_pair(a.uid_, b.uid_));
          else edges.push_back(std::make_pair(b.uid_, a.uid_));

          return Edge(this, a.uid_, b.uid_);
          }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    edges.clear();
    pts.clear();
    adjList.clear();
  }

 private:

    std::vector<Point> pts; //vector to store the nodes(points)
    std::vector<std::vector<size_type>> adjList; //Store the nodes in adjointed list format
    //aka for each node, we all know all nodes to which it is connected
    std::vector<std::pair<size_type, size_type>> edges; //store simple pairs of nodes, aka
    //"simple edges" - (i,j)
};

#endif // CME212_GRAPH_HPP
