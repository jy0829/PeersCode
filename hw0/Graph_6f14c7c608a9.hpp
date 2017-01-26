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
    
    class internal_node; // vector 
    class internal_edge; // 
    

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
      return fetch().point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        if (n.graph_ == graph_ && n.uid_ == uid_) {
            return true;
        }
      //(void) n;          // Quiet compiler warning
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
        if (uid_ < n.uid_)
            return true;
      (void) n;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
      graph_type* graph_;
      size_type uid_;
      
      Node(const graph_type* graph, size_type uid) 
          : graph_(const_cast<Graph*>(graph)), uid_(uid) {
      }
      // Helper method to return the appropriate node
      internal_node& fetch() const {
        return graph_->node_vec[uid_];
      }
      
      
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
    // HW0: YOUR CODE HERE
      internal_node node_new;
      node_new.point_ = position;
      node_new.uid_ = num_nodes();//+1;
      node_vec.push_back(node_new);
      size_++;    
    (void) position;      // Quiet compiler warning
    return Node(this,size_-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    internal_node temp_node;
    //temp_node.point_ = node_vec[n.uid_].point_;
    temp_node.uid_ = n.index();
      
    auto is_bool = std::find(node_vec.begin(), node_vec.end(), temp_node); 
    if (is_bool != node_vec.end() ) {
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
    // Line 91 from proxy example
    //(void) i;             // Quiet compiler warning
    return Node(this,i);        // Invalid node
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
      size_type node1_uid = graph_->edge_vec[uid_].node1;
      Node node1 = graph_->node(node1_uid);  
      return node1;
      //return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type node2_uid = graph_->edge_vec[uid_].node2;
      Node node2 = graph_->node(node2_uid);  
      return node2;
      //return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ == e.graph_ && uid_ == e.uid_) {
          return true;
      }
      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (uid_ < e.uid_) {
            return true;
      }
      (void) e;           // Quiet compiler warning
      return false;
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
    
    // Constructor (SimpleElement...)
    Edge(const Graph* graph , size_type uid) 
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
      
      
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    return Edge(this,i);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
      internal_edge temp_edge;
      temp_edge.node1 = a.index();
      temp_edge.node2 = b.index();
  
      auto is_true = std::find(edge_vec.begin(), edge_vec.end(), temp_edge); 
          
      if (is_true != edge_vec.end() ) {
          return true;
      }
      return false;
      (void) a; (void) b;   // Quiet compiler warning
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
    internal_edge edge_new;
    edge_new.node1 = a.index();
    edge_new.node2 = b.index();  
      
    auto is_bool = std::find(edge_vec.begin(), edge_vec.end(), edge_new); 
    if (is_bool == edge_vec.end() ) {
        edge_new.uid_ = edge_vec.size();
        edge_vec.push_back(edge_new);
        return Edge(this, edge_new.uid_);
    }  
      
    // if (!(a.index() > size_ || b.index() > size_) ) {
    //    edge_new.uid_ = edge_vec.size();
    //    edge_vec.push_back(edge_new);
    //    return Edge(this, edge_new.uid_);
    // } 

    (void) a, (void) b;   // Quiet compiler warning
    internal_edge edge_return = edge_vec[size_-1];  
    return Edge(this, edge_return.uid_);        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    size_ = 0;
    node_vec.clear();
    edge_vec.clear();

  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
    
  class internal_edge {
     public:
     size_type node1;
     size_type node2;
     size_type uid_;
     
      bool operator==(const internal_edge& e) const {
          if ( (node1 == e.node1 && node2 == e.node2) || (node2 == e.node1 && node1 == e.node2) ) {
              return true;
          }
          return false;
      } 
  };
    
  class internal_node {
      public:
      Point point_;
      size_type uid_;
      bool operator==(const internal_node& n) const {
        if ( uid_ == n.uid_ ) {
            return true;
        }
        return false;
      } 
  };
    
  std::vector<internal_node> node_vec;
  std::vector<internal_edge> edge_vec;
  size_type size_;
              
};

#endif // CME212_GRAPH_HPP
