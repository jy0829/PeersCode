#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP


/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>

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
    size_ = 0;
    next_id_ = 0;
    size_edges_ = 0;
    next_e_id_ = 0;
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[id_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.index() != id_) {
        return false;
      }
      if (n.graph_ != graph_) {
        return false;
      }
      return true;
    }

    /** Test whether this node and @a n are inequal.
     *
     * inequal nodes either don't have the same graph or
     * don't have the same index.
     */
    bool operator!=(const Node& n) const{
      if (n.index() != id_) {
        return true;
      }
      if (n.graph_ != graph_) {
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
      if (graph_ < n.graph_) {
        return true;
      }
      if (n.graph_ < graph_) {
        return false;
      }
      if (id_ < n.index()) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    const graph_type* graph_;
    size_type id_;
    /** Construct a valid node.
     */
    Node(const graph_type* graph, size_type id) {
      graph_ = graph;
      id_ = id; 
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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
    Node n(this, next_id_);
    ++next_id_;
    ++size_;
    nodes_.push_back(position);
    std::vector<std::pair<size_type, size_type>> v;
    edge_adj_.push_back(v);
    return n;     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this) {
      if (size_ > n.id_) {
        return true;
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, end_points_.first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, end_points_.second);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ != e.graph_) {
        return false;
      }
      if ((end_points_.first == e.end_points_.first) and 
        (end_points_.second == e.end_points_.second)) {
        return true;
      }
      if ((end_points_.first == e.end_points_.second) and 
        (end_points_.second == e.end_points_.first)) {
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
      if (graph_ < e.graph_) {
        return true;
      }
      if (e.graph_ < graph_) {
        return false;
      }
      if (id_ < e.id_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const graph_type* graph_;
    size_type id_;
    std::pair<size_type, size_type> end_points_;
    /** Create valid edge
     */
    Edge(const graph_type* graph, size_type id, Node a, Node b) {
      id_ = id;
      graph_ = graph;
      end_points_.first = a.id_;
      end_points_.second = b.id_;
    }
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return size_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    Node n1(this, edges_[i].first);
    Node n2(this, edges_[i].second);
    Edge e(this, i, n1, n2);
    return e;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned i = 0; i < edge_adj_[a.id_].size(); i++) {
      if (edge_adj_[a.id_][i].first == b.id_) {
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
    if (has_edge(a, b)) {
      size_type index = find_edge_index(a, b); 
      return Edge(this, index, a, b);
    }
    else {
      Edge e(this, next_e_id_, a, b);
      std::pair<size_type, size_type> p1;
      std::pair<size_type, size_type> p2;
      std::pair<size_type, size_type> p3;
      p1 = std::make_pair(b.id_, next_e_id_);
      p2 = std::make_pair(a.id_, next_e_id_);
      p3 = std::make_pair(a.id_, b.id_);
      edge_adj_[a.id_].push_back(p1);
      edge_adj_[b.id_].push_back(p2);
      edges_.push_back(p3);
      ++next_e_id_;
      ++size_edges_; 
      return e;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    size_ = 0;
    next_id_ = 0;
    next_e_id_ = 0;
    size_edges_ = 0;
    nodes_.clear();
    edges_.clear(); 
    edge_adj_.clear();
  }

 private:
  size_type size_;
  size_type next_id_;
  size_type next_e_id_;
  size_type size_edges_;
  std::vector<Point> nodes_;
  std::vector<std::vector<std::pair<size_type, size_type>>> edge_adj_;
  std::vector<std::pair<size_type, size_type>> edges_;
  /** Find the index of an edge in a graph.
  * @pre @a a and @a b are distinct valid nodes of this graph for which 
  * has_edge(@a a, @a b) == true 
  * @return @a i, where edges_[i] is a pair containing @a a.id_ and @a b.id_ 
  */
   size_type find_edge_index(const Node& a, const Node& b) {
     for (unsigned i = 0; i < edge_adj_[a.id_].size(); i++) {
       if (edge_adj_[a.id_][i].first == b.id_) {
         return edge_adj_[a.id_][i].second;
       }
     }
     return size_edges_ + 1;
   }
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
