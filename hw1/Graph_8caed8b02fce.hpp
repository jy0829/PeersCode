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
template <typename V>
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
  using graph_type = Graph<V>;
 
  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Type of Node value which the graph can accept */
  using node_value_type = V; 

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

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
  class Node : private totally_ordered<Node> {
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
    Node(){ 
      // HW0: YOUR CODE HERE
    }
 
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g_->points[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value, of type node_value_type. */
    node_value_type& value() {
      return g_->values[index_];
    };
    
    /** Return this node's value, of type node_value_type. Similar to above, 
    *   but used whenever this node is a constant object.
    */
    const node_value_type& value() const {
      return g_->values[index_];
    };

    /** Return this node's degree */
    size_type degree() const {
      return g_->sp_adj_matrix[index_].size();
    };
    
    /** Create an iterator object to loop through the edges in this node's
    *   neighborhood. Returns an iterator pointing to the first edge in the
    *   neighborhood.
    */
    IncidentIterator edge_begin() const {
      return IncidentIterator(g_, index_, 0);
    };
    
    /** Create an iterator object to loop through the edges in this node's 
    *   neighborhood. Returns an iterator pointing to the last edge in this 
    *   node's neighborhood.
    */
    IncidentIterator edge_end() const {
      return IncidentIterator(g_, index_, degree());
    };

    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return g_ == n.g_ && index_ == n.index_;
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
      if(g_ != n.g_) {
        return g_ < n.g_;
      } else {
        return index_ < n.index_;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(const Graph* g, size_type id): g_(const_cast<Graph*>(g)), index_(id) {}
    //Store a graph pointer and the index of the node in the graph
    Graph* g_;
    size_type index_; 
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return points.size();
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
  Node add_node(const Point& position, 
                  const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    points.push_back(position);
    values.push_back(val);
    std::vector<size_type> neighborhood_vec;
    sp_adj_matrix.push_back(neighborhood_vec);
    return Node(this, points.size()-1);     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.g_ == this && n.index_ <= points.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i);        // Invalid node
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
      return Node(g_, index1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g_, index2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool same_nodes = (index1_ == e.index1_ && index2_ == e.index2_) ||
      (index1_ == e.index2_ && index2_ == e.index1_);
      return g_ == e.g_ && same_nodes;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(g_ != e.g_) {
        return g_ < e.g_;
      } else {
        if(std::min(index1_, index2_) != std::min(e.index1_, e.index2_)) {
          return std::min(index1_, index2_) < std::min(e.index1_, e.index2_);
        } else {
          return std::max(index1_, index2_) < std::max(e.index1_, e.index2_);
        }
      }  
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* g, Node n1, Node n2): g_(g), index1_(n1.index()), 
       index2_(n2.index()) {}
    const Graph* g_;
    size_type index1_;
    size_type index2_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return std::distance(edge_begin(), edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for(uint k = 0; k < sp_adj_matrix[a.index_].size(); ++k) {
      if(sp_adj_matrix[a.index_][k] == b.index_) return true;
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
    if(not has_edge(a, b)) {
      sp_adj_matrix[a.index_].push_back(b.index_);
      sp_adj_matrix[b.index_].push_back(a.index_);
    }
    return Edge(this, a, b);    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    points.clear();
    sp_adj_matrix.clear();
    values.clear();
  }
  

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    
    /** Deference node iterator. Returns the node object at the current
    *   position.
    */
    Node operator*() const {
      return Node(gptr_, index_);
    }

    /** Increment current position. */
    NodeIterator& operator++() {
      ++index_;
      return *this;
    }
  
    /** Test iterators for equality. Returns true if iterators belong to 
    *   the same graph and have the same current position.
    */
    bool operator==(const NodeIterator& iter) const {
      return gptr_ == iter.gptr_ && index_ == iter.index_;
    }
  
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    /** Construct iterator object. Maintains a graph pointer and a pointer
    *   to the current position.
    */
    NodeIterator(const Graph* g, size_type node_id): gptr_(g), index_(node_id)  {}
    const Graph* gptr_;
    size_type index_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /** Create a node iterator. Returns an iterator pointing to a node object
  *   with index 0.
  */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Create a node iterator. Returns a node iterator pointing to a node 
  *   node object with index equal to num_nodes().
  */
  NodeIterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    
    /** Dereference iterator. Returns edge object pointed to by current 
    *   position in adjacency list.
    */
    Edge operator*() const {
      size_type idx2 = gptr_->sp_adj_matrix[nodeidx_][iteridx_];
      return Edge(gptr_, Node(gptr_, nodeidx_), Node(gptr_, idx2));
    };
    
    /** Increment this iterator's current position and return this iterator. */
    IncidentIterator& operator++() {
      ++iteridx_;
      return *this;
    }
  
    /** Compare two iterator objects. Returns true if the iterators belong to 
    *   the same graph and have the same current position.
    */
    bool operator==(const IncidentIterator& iter) const {
      return gptr_ == iter.gptr_ && nodeidx_ == iter.nodeidx_ 
                    && iteridx_ == iter.iteridx_; 
    } 

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    /** Construct incidence iterator objects. Maintains a pointer to the parent
    *   graph and the current position as a row and column index for the 
    *   edge adjacency list. 
    */
    IncidentIterator(const Graph* g, size_type node_id, size_type starter_id):
         gptr_(g), nodeidx_(node_id) {iteridx_ = starter_id;}
    const Graph* gptr_;
    size_type nodeidx_;
    size_type iteridx_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    /** Dereference edge iterator. Returns edge object pointed to by this 
    *   iterator's current position.
    */
    Edge operator*() const {
      size_type indx2 = gptr_->sp_adj_matrix[nodeidx_][iteridx_]; 
      return Edge(gptr_, Node(gptr_, nodeidx_), Node(gptr_, indx2)); 
    } 
    
    /** Increment this iterator's current position and return this iterator. */
    EdgeIterator& operator++() {
      ++iteridx_;
      fix();
      return *this;
    }
  
    /** Compare two iterators. Returns true if iterators belong to the same
    *   graph and have the same current position.
    */
    bool operator==(const EdgeIterator& iter) const {
      return gptr_ == iter.gptr_ && nodeidx_ == iter.nodeidx_ 
                      && iteridx_ == iter.iteridx_;
    }
     
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    /** Construct edge iterator objects. Maintain a graph pointer and a 
    *   current position as a row and column index for the adjacency list.
    */
    EdgeIterator(const Graph* g, size_type node_id, size_type iter_id) : 
                   gptr_(g), nodeidx_(node_id), iteridx_(iter_id) {fix();}
    const Graph* gptr_;
    size_type nodeidx_;
    size_type iteridx_;
    void fix() {
      while(true) {
        if(nodeidx_ == gptr_->sp_adj_matrix.size()) {
          break;
        } else if(iteridx_ >= gptr_->sp_adj_matrix[nodeidx_].size()) {
          iteridx_ = 0;
          ++nodeidx_;
        } else if(nodeidx_ < gptr_->sp_adj_matrix[nodeidx_][iteridx_]) {
          ++iteridx_;
        } else {
          break;
        };
      };  
    };
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  /** Returns an edge iterator object pointing to the edge whose node
  *   endpoints have the smallest indices in lexicographical order. 
  */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  };

  /** Returns an edge iterator object pointing to the edge whose node 
  *   endpoints have the largest indices in lexicographical order.
  */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, num_nodes(), 0);
  };

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  //The Graph class has two private data members. A vector of Point objects
  //stores the locations of the Nodes in the graph using the proxy design 
  //pattern. An adjacency list (vector of vectors) stores the indices of
  //the nodes which are connected in the graph.
   std::vector<Point> points;
   std::vector<std::vector<size_type>> sp_adj_matrix;
   std::vector<node_value_type> values;
};

#endif // CME212_GRAPH_HPP
