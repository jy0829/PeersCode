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
template <typename V, typename E>
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
  using graph_type = Graph<V, E>;
 
  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Type of Node value which the graph can accept */
  typedef V node_value_type;

  /** Predeclaration of Edge type. */
  class Edge;

  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  
  /** Type of Edge value which the graph can accept */
  typedef E edge_value_type;

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
    }
    
    /** Return this node's current position */
    Point& position() {
      //assert(valid());
      return g_->node_list[uid_].point_;
    }
     
    /** Return this node's position. */
    const Point& position() const {
      return g_->node_list[uid_].point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return g_->node_list[uid_].idx_;
    }

    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value, of type node_value_type. */
    node_value_type& value() {
      return g_->node_list[uid_].val_;
    };
    
    /** Return this node's value, of type node_value_type. Similar to above, 
    *   but used whenever this node is a constant object.
    */
    const node_value_type& value() const {
      return g_->node_list[uid_].val_;
    };

    /** Return this node's degree */
    size_type degree() const {
      return g_->adj_matrix[g_->node_list[uid_].idx_].size();
    };
    
    /** Create an iterator object to loop through the edges in this node's
    *   neighborhood. Returns an iterator pointing to the first edge in the
    *   neighborhood.
    */
    IncidentIterator edge_begin() const {
      return IncidentIterator(g_, uid_, 0);
    };
    
    /** Create an iterator object to loop through the edges in this node's 
    *   neighborhood. Returns an iterator pointing to the last edge in this 
    *   node's neighborhood.
    */
    IncidentIterator edge_end() const {
      return IncidentIterator(g_, uid_, degree());
    };

    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return g_ == n.g_ && uid_ == n.uid_;
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
        return uid_ < n.uid_;
      }
    }

    /** Check whether this node has not been invalidated by the remove method */
    bool valid() const {
      return uid_ >= 0 && uid_ < g_->node_list.size() 
             && g_->node_list[uid_].idx_ < g_->i2u.size()
             && g_->i2u[g_->node_list[uid_].idx_] == uid_;
    }
    
    Graph* g_;
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(const Graph* g, size_type id): g_(const_cast<Graph*>(g)), uid_(id) {}
    //Store a graph pointer and the index of the node in the graph
    size_type uid_; 
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u.size();
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
    node_list.push_back(NodeInfo(position, val, num_nodes()));
    //Initialize i2u vector s.t. i2u[i] = i
    i2u.push_back(node_list.size()-1);
    adj_matrix.push_back(std::vector<EdgeInfo>());
    return Node(this, i2u[num_nodes()-1]);     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.valid(); 
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    //assert(node_list[i2u[i]].idx_ == i);
    return Node(this, i2u[i]);
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(g_, uid1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g_, uid2_);      // Invalid Node
    }

    /** Return the value of this edge */
    edge_value_type& value() {
      size_type outeridx = g_->node_list[std::max(uid1_, uid2_)].idx_;
      size_type inner = std::min(uid1_, uid2_);
      for(uint i = 0; i < g_->adj_matrix[outeridx].size(); ++i) 
        if(g_->adj_matrix[outeridx][i].uid2_ == inner) 
          return g_->adj_matrix[outeridx][i].val_;
    }
 
    /** Return the value of this edge */
    const edge_value_type& value() const {
      size_type outeridx = g_->node_list[std::max(uid1_, uid2_)].idx_;
      size_type inner = std::min(uid1_, uid2_);
      for(uint i = 0; i < g_->adj_matrix[outeridx].size(); ++i) 
        if(g_->adj_matrix[outeridx][i].uid2_ == inner) 
          return g_->adj_matrix[outeridx][i].val_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool same_nodes = (uid1_ == e.uid1_ && uid2_ == e.uid2_) ||
      (uid1_ == e.uid2_ && uid2_ == e.uid1_);
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
        if(std::min(uid1_, uid2_) != std::min(e.uid1_, e.uid2_)) {
          return std::min(uid1_, uid2_) < std::min(e.uid1_, e.uid2_);
        } else {
          return std::max(uid1_, uid2_) < std::max(e.uid1_, e.uid2_);
        }
      }  
    }
   
    /** Return the length of the edge defined as the Euclidean distance 
    *   between its incident nodes.
    */
    double length() const;
 
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* g, Node n1, Node n2): g_(const_cast<Graph*>(g)), 
         uid1_(g_->i2u[n1.index()]), uid2_(g_->i2u[n2.index()]) {}
    Graph* g_;
    size_type uid1_;
    size_type uid2_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
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
    for(uint k = 0; k < adj_matrix[a.index()].size(); ++k) {
      if(adj_matrix[a.index()][k].uid2_ == b.uid_) return true;
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
  Edge add_edge(const Node& a, const Node& b, 
                  const edge_value_type& val = edge_value_type()) {
    if(not has_edge(a, b)) {
      adj_matrix[a.index()].push_back(EdgeInfo(b.uid_, val));
      adj_matrix[b.index()].push_back(EdgeInfo(a.uid_, val));
    }
    return Edge(this, a, b);    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_list.clear();
    i2u.clear();
    adj_matrix.clear();
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

    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    
    /** Deference node iterator. Returns the node object at the current
    *   position.
    */
    Node operator*() const {
      return Node(gptr_, uid_);
    }

    /** Increment current position. */
    NodeIterator& operator++() {
      ++uid_;
      return *this;
    }
  
    /** Test iterators for equality. Returns true if iterators belong to 
    *   the same graph and have the same current position.
    */
    bool operator==(const NodeIterator& iter) const {
      return gptr_ == iter.gptr_ && uid_ == iter.uid_;
    }
  
   private:
    friend class Graph;
    /** Construct iterator object. Maintains a graph pointer and a pointer
    *   to the current position.
    */
    NodeIterator(const Graph* g, size_type node_id): gptr_(g), uid_(node_id)  {}
    const Graph* gptr_;
    size_type uid_;
  };

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
      size_type idx2 = gptr_->adj_matrix[gptr_->node_list[uid_].idx_]
                                        [iteridx_].uid2_;
      return Edge(gptr_, Node(gptr_, uid_), Node(gptr_, idx2));
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
      return gptr_ == iter.gptr_ && uid_ == iter.uid_ 
                    && iteridx_ == iter.iteridx_; 
    } 

   private:
    friend class Graph;
    /** Construct incidence iterator objects. Maintains a pointer to the parent
    *   graph and the current position as a row and column index for the 
    *   edge adjacency list. 
    */
    IncidentIterator(const Graph* g, size_type node_id, size_type starter_id):
         gptr_(g), uid_(node_id), iteridx_(starter_id) {}
    const Graph* gptr_;
    size_type uid_;
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

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    /** Dereference edge iterator. Returns edge object pointed to by this 
    *   iterator's current position.
    */
    Edge operator*() const {
      size_type id2 = gptr_->adj_matrix[nodeidx_]
                                       [iteridx_].uid2_; 
      return Edge(gptr_, Node(gptr_, gptr_->i2u[nodeidx_]), Node(gptr_, id2)); 
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
    /** Construct edge iterator objects. Maintain a graph pointer and a 
    *   current position as a row and column index for the adjacency list.
    */
    EdgeIterator(const Graph* g, size_type node_idx, size_type iter_id) : 
                   gptr_(g), nodeidx_(node_idx), 
                   iteridx_(iter_id) {fix();}
    const Graph* gptr_;
    size_type nodeidx_;
    size_type iteridx_;
    /** Move the iterator position to valid indices in the adjacency list. */
    void fix() {
      while(true) {
        if(nodeidx_ == gptr_->num_nodes()) {
          break;
        } else if(iteridx_ >= gptr_->adj_matrix[nodeidx_].size()) {
          iteridx_ = 0;
          ++nodeidx_;
        //Increment position only if indices lie in lower triangle of list
        } else if(gptr_->i2u[nodeidx_] > gptr_->adj_matrix[nodeidx_][iteridx_].uid2_) {
          ++iteridx_;
        } else {
          break;
        };
      };  
    };
  };

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
  

  /** Remove the given node n.
   *  @param[in, out] n
   *  @post n.valid() != true 
   *  @post new num_nodes() == old num_nodes() - 1  
   *  @post adj_matrix[n.index()] == adj_matrix[old num_nodes()-1]
   *  @post i2u[n.index()] == i2u[old num_nodes()-1]
   *  @result == 1 if @in n.valid() == true
   *  @result == 0, otherwise
   */
  size_type remove_node(const Node& n) {
    if(has_node(n)) {
      for(auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
        size_type neighidx = (*iit).node2().index();
        for(uint i = 0; i < adj_matrix[neighidx].size(); ++i) {
          if(adj_matrix[neighidx][i].uid2_ == n.uid_) {
            adj_matrix[neighidx][i] = adj_matrix[neighidx][adj_matrix[neighidx].size()-1];
            adj_matrix[neighidx].pop_back();
          }
        }
      }    
      adj_matrix[n.index()] = adj_matrix[num_nodes()-1];
      adj_matrix.pop_back();
      i2u[n.index()] = i2u[num_nodes()-1];
      node_list[i2u[num_nodes()-1]].idx_ = n.index();
      i2u.pop_back();
      for(uint i = 0; i <  num_nodes(); i++) 
        assert(node_list[i2u[i]].idx_ == i); 
      assert(!n.valid());
      return 1;
    }
    return 0;
  }

   /** Remove the given node n.
   *  @param[in, out] n
   *  @pre n.valid() == true
   *  @post n.valid() != true 
   *  @post new num_nodes() == old num_nodes() - 1  
   *  @post adj_matrix[n.index()] == adj_matrix[old num_nodes()-1]
   *  @post i2u[n.index()] == i2u[old num_nodes()-1]
   *  @result == node_begin() if @in n.valid() == true
   *  @result == node_end() otherwise
   */
  NodeIterator remove_node(NodeIterator n_it) {
    Node n = (*n_it);
    if(has_node(n)) {
      for(auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
        size_type neighidx = (*iit).node2().index();
        for(uint i = 0; i < adj_matrix[neighidx].size(); ++i) {
          if(adj_matrix[neighidx][i].uid2_ == n.uid_) {
            adj_matrix[neighidx][i] = adj_matrix[neighidx][adj_matrix[neighidx].size()-1];
            adj_matrix[neighidx].pop_back();
          }
        }
      }    
      adj_matrix[n.index()] = adj_matrix[num_nodes()-1];
      adj_matrix.pop_back();
      i2u[n.index()] = i2u[num_nodes()-1];
      node_list[i2u[num_nodes()-1]].idx_ = n.index();
      i2u.pop_back();
      for(uint i = 0; i <  num_nodes(); i++) 
        assert(node_list[i2u[i]].idx_ == i); 
      assert(!n.valid());
      return node_begin();
    }
    return node_end();
  }
 
  /** Remove a given edge from the graph. 
   *  @param[in] n1, n2
   *  @post new num_edges() == old num_edges() - 1
   */ 
  size_type remove_edge(const Node& n1, const Node& n2) {
    if(has_edge(n1, n2)) {
      for(uint i = 0; i < adj_matrix[n1.index()].size(); ++i) {
        if(adj_matrix[n1.index()][i].uid2_ == n2.uid_) {
          adj_matrix[n1.index()][i] = adj_matrix[n1.index()][n1.degree()-1];
          adj_matrix[n1.index()].pop_back();
        }
       }
       for(uint i = 0; i < adj_matrix[n2.index()].size(); ++i) {
         if(adj_matrix[n2.index()][i].uid2_ == n1.uid_) { 
          adj_matrix[n2.index()][i] = adj_matrix[n2.index()][n2.degree()-1];
          adj_matrix[n2.index()].pop_back();
         }
       }
       return 1;
    }
    return 0;
  }
 /** See above */ 
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }
  
  /** See above */
  edge_iterator remove_edge(edge_iterator e_it) {
    Node n1 = (*e_it).node1();
    Node n2 = (*e_it).node2();
    if(has_edge(n1, n2)) {
     for(uint i = 0; i < adj_matrix[n1.index()].size(); ++i) {
        if(i == n2.uid_) {
          adj_matrix[n1.index()][i] = adj_matrix[n1.index()][n1.degree()-1];
          adj_matrix[n1.index()].pop_back();
        }
       }
       for(uint i = 0; i < adj_matrix[n2.index()].size(); ++i) {
         if(i == n1.uid_) { 
          adj_matrix[n2.index()][i] = adj_matrix[n2.index()][n2.degree()-1];
          adj_matrix[n2.index()].pop_back();
         }
       }
       return edge_begin();
    }
    return edge_end();
  }

 
 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  //The Graph class has three private data members. A vector of NodeInfo 
  //objects stores the data associated to each node in the graph using the 
  //proxy design pattern. An adjacency list (vector of vectors) stores the 
  //indices of the nodes which are adjacent in the graph. A vector of
  //EdgeInfo objects stores the data associated to each edge. Edges are 
  //identified uniquely via the EdgeIterator, where each edge is always 
  //chosen to lie on the lower triangle of the adjacency list.
  
  /** This object stores the data associated to each node */
  struct NodeInfo {
    NodeInfo(Point point, node_value_type val, size_type idx) : 
      point_(point), val_(val), idx_(idx) {}
    Point point_;
    node_value_type val_;
    size_type idx_;
  };
 
  struct EdgeInfo {
    size_type uid2_;
    edge_value_type val_;
    EdgeInfo(size_type uid2, edge_value_type val) : uid2_(uid2), val_(val) {}
  };
   
   std::vector<NodeInfo> node_list;
   std::vector<size_type> i2u;
   std::vector<std::vector<EdgeInfo>> adj_matrix;
};

#endif // CME212_GRAPH_HPP

