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

  //node internal struct
  struct node_info;

  //edge interal struct
  struct edge_info;

  //node array for graph
  std::vector<node_info> nodes_vec;

  //edge array for graph
  std::vector<edge_info> edges_vec;

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
  Graph() : nodes_vec(), edges_vec() {
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
      //no additiontional code needed to construct an invalid node
    }

    /** Return this node's position. */
    const Point& position() const { 
     //get node struct from graph and extract position   
      return n_graph_->nodes_vec[node_id_].P;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //return node_id (which is a member of the Node class)
      return node_id_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // check if they have the same graph and index
      if (n.node_id_ == node_id_ && n.n_graph_ == n_graph_) {
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
      //if they are in the same graph, order based on node id
      if ( n.n_graph_ == n_graph_) {
        if (node_id_ < n.node_id_) {
          return true;
        }
      }
      /*if different graphs, global order will depend on graph location 
      in memory, i.e. if mem location of graph of node n is larger than the
      mem location of this graph, node n is larger in global order, regardless 
      of node id */
      else if (n_graph_ <  n.n_graph_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    //pointer back to the Graph (8 bits)
    const Graph* n_graph_;

    //node id (8 bits)
    size_type node_id_;

  }; //end Node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    //get length of vector holding unique node info
    size_type total_nodes = nodes_vec.size();
    return total_nodes;
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
    //create and populate struct for node being added
    node_info current_node;
    current_node.P = position; 
    //add this node struct to the graph
    nodes_vec.push_back(current_node);
    //create node to be returned and populate
    Node N;
    N.n_graph_ = this;
    //get node id
    N.node_id_ = nodes_vec.size()-1;
    
    return N;     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //current number of nodes in graph
    size_type no_nodes = nodes_vec.size();
    //make sure graph is the same, evaluate if index is larger than vector size
    if (no_nodes > n.node_id_ && n.n_graph_ == this) {
      return true;
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
    //node to be returned
    Node node_at_i;
    //assign node id
    node_at_i.node_id_ = i;
    //assign graph
    node_at_i.n_graph_ = this;

    return node_at_i;   
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
      //no additional code needed to construct an invalid edge
    }

    /** Return a node of this Edge */
    Node node1() const {
      //construct node 1
      Node n1;
      n1.node_id_ = n_1_id_;
      n1.n_graph_ = e_graph_;
      return n1; 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //construct node 2
      Node n2;
      n2.node_id_ = n_2_id_;
      n2.n_graph_ = e_graph_;
      return n2;   
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //check if (a,b)=(a,b) or (a,b)=(b,a) using min and max node indices
      size_type min_node = std::min(e.n_1_id_, e.n_2_id_);
      size_type max_node = std::max(e.n_1_id_, e.n_2_id_);
      if (std::min(n_1_id_, n_2_id_) == min_node && std::max(n_1_id_, n_2_id_) == max_node) {
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
      //for edges in the same graph, will be determined by smallest node index. 
      //If they have the same min node index, will then be decided by smallest
      //second node.
      if (e_graph_ == e.e_graph_) {
        if (std::min(n_1_id_, n_2_id_) < std::min(e.n_1_id_, e.n_2_id_)) {
          return true;
        }
        else if (std::min(n_1_id_, n_2_id_) == std::min(e.n_1_id_, e.n_2_id_)) {
          if (std::max(n_1_id_, n_2_id_) < std::max(e.n_1_id_, e.n_2_id_)) {
            return true;
          }
        }
      }
      /* if they aren't in the same graph, use smallest graph memory location
      (similar to node global order */
      else if (e_graph_ < e.e_graph_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    //pointer to graph (8 bits)
    const Graph* e_graph_;

    //node 1 id (8 bits)
    size_type n_1_id_;

    //node 2 id (8 bits)
    size_type n_2_id_;

  }; //end Edge class

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
   //get length of vector holding unique node info
    size_type total_edges = edges_vec.size();
    return total_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // get edge struct from graph (for edge i)
    edge_info edge_at_i = edges_vec[i];
    //create edge from the struct
    Edge E;
    //add nodes
    E.n_1_id_ = edge_at_i.N_1_id;
    E.n_2_id_ = edge_at_i.N_2_id;
    //add graph
    E.e_graph_ = this;
    return E;    
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //scan all existing edges
    for (size_type i=0; i < edges_vec.size(); i++) {
      //get node ids from struct in edges vector
      edge_info vec_edge = edges_vec[i];
      //minimum node id of this pair
      size_type n_min = std::min(vec_edge.N_1_id, vec_edge.N_2_id);
      //maximum node id of this pair
      size_type n_max = std::max(vec_edge.N_1_id, vec_edge.N_2_id);
      //check if this is equal to the nodes in the input edge
      if (std::min(a.node_id_, b.node_id_) == n_min && std::max(a.node_id_, b.node_id_) == n_max) {
        //case where the two input nodes connect an edge in the graph
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
    //check if edge currently exists
    if (has_edge(a, b)) {
      //create edge from these nodes
      Edge e_return;
      e_return.n_1_id_ = a.node_id_;
      e_return.n_2_id_ = b.node_id_;
      e_return.e_graph_ = this;
      return e_return;
    }
    //if not, build edges struct and add it to the graph
    edge_info edge_to_add;
    edge_to_add.N_1_id = a.node_id_;
    edge_to_add.N_2_id = b.node_id_;
    edges_vec.push_back(edge_to_add);
    //create edge from these nodes
    Edge e_new;
    e_new.n_1_id_ = a.node_id_;
    e_new.n_2_id_ = b.node_id_;
    e_new.e_graph_ = this;

    return e_new;        
  }

 
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // clear the vectors
    nodes_vec.clear();
    edges_vec.clear();
  }

 private:

  //interal type for node
  struct node_info {

    Point P; //Cartesian location of the node
  };

  //interal type for edge
  struct edge_info {
    
    size_type N_1_id; //id of node 1 of the edge
    size_type N_2_id; //id of node 2 of the edge
   
  };

};

#endif // CME212_GRAPH_HPP
