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

// Citations/Collaboration:
// Worked with Ian Shaw, Jonathon Roth, and Amery Martin directly
// Attended office hours with Amy Shoemaker (TA)
// Referenced (and followed the structure closely) of proxy_example.cpp
// Also referenced public github reposititories for previous students in 
// Harvard CS207: Lan, Tian; Piasevoli, Filip; Tran, Dustin; Zacarias, Lisa; 
// **Code was not copied direcly, only referenced for ideas when I was stuck

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
    
class Graph {
 private:
  // Predeclaring the internal structures: internal_node and internal_edge
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
  Graph() {
    // Construct a Graph with zero size
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
        return fetch().point;
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
      // Check for same graph and uid
        if (n.graph_ == graph_ && n.uid_ == uid_)
            return true;
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
      // Assign order by uid
        if (uid_ < n.uid_)
            return true;
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // This space declares private data members and methods for Node
    // that are not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // **Using proxy_example.cpp as guide**:
    // Pointer back to the graph_type container
      graph_type* graph_;
    // This element's unique identification number
      size_type uid_;
    // Private Constructor 
      Node(const graph_type* graph, size_type uid)
          : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
      }
    // Get internal node that is stored in graph
    /** Helper method to return the appropriate node.
     * This checks if the uid is less than graph size and returns the
     * proper uid
     */
      internal_node& fetch() const {
          // We do not wish to loop over all nodes
          // We have commented out the section from proxy_example.cpp
          // we eliminate the for loop here
          //for (size_type i = 0; i < set_->size(); ++i)
              //if (set_->nodes_[i].uid == uid_)
                  //return set_->nodes_[i];
          //assert(false);
          if (uid_ < graph_->size())
              return graph_->nodes[uid_];
          assert(false);
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  // Return Graph's size
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
    // std vector of nodes
      internal_node new_node;
    // Each element of the node vector has two features: point & uid
      new_node.point = position;
      new_node.uid = size_;
    // Use push_back function to append the new node to vector nodes
      nodes.push_back(new_node);
    // Increment our size 
      ++size_;
    // Returns a node that points to the new node
      return Node(this, size_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // n.uid is unsigned, so we need not check non-negativity here
      if ((n.uid_ < nodes.size()) && (n.index() < size())) 
          return true;
      return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // Check that i is appropriate
    //  assert (i < size());
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
      // No code needed at this time
    }

    /** Return a node of this Edge */
    Node node1() const {
      // Return appropriate node of the edge corresponding to uid
        return graph_->node(graph_->edges[uid_].node1); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Return appropriate node2 
        return graph_->node(graph_->edges[uid_].node2); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Follow structure for equality operator in Node
        if (graph_ == e.graph_ && uid_ == e.uid_)
            return true;
        return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //  Follow structure for ordering as in class Node
      //  Order by uid
        if (uid_ < e.uid_)
            return true;
        return false;
    }

   private:
      friend class Graph;
    // Pointer back to the Graph container
      graph_type* graph_;
    // This edge's unique identification number
      size_type uid_;
    // Private Constructor
      Edge(const graph_type* graph, size_type uid)
          : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // Return the total number of edges
      return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Check if i is appropriate
    //  assert(i < size()); // (from proxy_example)
    // Return new Edge
      return Edge(this, i);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Define temporary internal_edge vector to see if a and b form an edge
    // and call this vector is_edge
      internal_edge is_edge;
    // Assign index of a and b to node1 and node2 respectively in the is_edge vector
      is_edge.node1 = a.index();
      is_edge.node2 = b.index();
    // Search the vector edges to see if is_edge exists
      if (std::find(edges.begin(), edges.end(), is_edge) == edges.end())
          return false;
      return true;
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
    // Returns an Edge that points to a new edge
    // Initialize new_edge of type internal_edge and assign node1 and node2 
    // the values of a and b respectively
      internal_edge new_edge;
      new_edge.node1 = a.index();
      new_edge.node2 = b.index();
    // Initialize our vector iterator of type internal_edge
      std::vector<internal_edge>::iterator it;
    // Check if edge is in edges using find function
      it = std::find(edges.begin(), edges.end(), new_edge);
    // Add edge if it does not exist
      if (it == edges.end()){
        // If iterator is not in edges, add it to the edge list
          new_edge.uid = edges.size();
        // Append the new edge to vector edges
          edges.push_back(new_edge);
        // Return the updated edge class and appropriate uid
          return Edge(this, new_edge.uid);
      }
      else{
        // If the new edge is not a new edge, update to the pointer of this edge
          return Edge(this, (*it).uid);
      }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      nodes.clear();
      edges.clear();
      size_ = 0;
  }

 private:
  // Internal type for set nodes
    struct internal_node {
        Point point;
        size_type uid;
    };
  // Internal type for set edges
    struct internal_edge {
        public:
            size_type node1;
            size_type node2;
            size_type uid;
            bool operator == (const internal_edge& n) const {
                if ((node1 == n.node1 && node2 == n.node2) || 
                    (node1 == n.node2 && node2 == n.node1))
                    return true;
                return false;
            }
    };
        

    std::vector<internal_node> nodes;
    std::vector<internal_edge> edges;
    size_type size_;

  // Disable copy and assignment of a Graph (from proxy_example.cpp)
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;
};

#endif // CME212_GRAPH_HPP
