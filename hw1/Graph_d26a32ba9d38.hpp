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

template<typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)


  /**It has data variables:
   * vector<Point> v_set: stores the Point variables corresponding to nodes.
   * vector<V> v_val: stores values assigned to nodes.
   * vector<Edge> e_set: stores the Edge objects.
   * struct Arr adj: An adjacency list for the graph.
   * struct Arr e_ind: same as adjacency list but stores edge indices. 
   * The data variables represent the mathematical graph.
   */

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

  using node_value_type = V;
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
  class Node: private totally_ordered<Node> {
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
      *@post Point @a p == @a graph.v_set[index()]
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

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /**Returns reference to the value of a given Node object.
    *@returns a reference to value of Node.
    */
    node_value_type& value() {
      return const_cast<Graph*>(graph)->v_val[index()];
    }
    /**Returns reference to the value of a given Node object.
    *@returns a reference to value of Node, but the value cannot be changed.
    */
    const node_value_type& value() const {
      return graph->v_val[index()];
    }
    /**Returns the degree of the node.
    *@returns a size_type value.
    *@post Let @a S be the set of all edges @a e such that 
           @a this == @a e.node1() or @a this == @a e.node2(),
           Then @a result = number of elements in @a S.      
    */
    size_type degree() const {
      return this->graph->adj(index());
    }    
    /**Returns the an incident iterator at the start of the neighbors of the node.
    * @returns an IncidentIterator
    * @post @a result->node1() == @a this.      
    */
    IncidentIterator edge_begin() const {
      //We start here at 1 because the first value is always 0.
      return IncidentIterator(this->graph,*this,0);
    }

    /**Returns the an incident iterator at the end of the neighbors of the node.
    * @returns an incident_iterator.
    * @post @a result->node1() == @a this.     
    */
    IncidentIterator edge_end() const {
      return IncidentIterator(this->graph,*this,this->degree());
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
      /**Checks if the node is less than node @a n by comparing their graphs and indices.
      *@param n: Node we compare with.
      *@return bool @a b, where
             b = true iff (@a (long int)graph < @a (long int)n.graph) or
                   (@a graph == @a n.graph) and @a this->index() < @a n.index()).
      *@post a global trichotomy of nodes is established.
      */
      if ((long int)graph < (long int)(n.graph)) return true;
      if (graph == n.graph && index() < n.index()) {
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
    return v_set.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    v_set.push_back(position);
    Node node(this,num_nodes()-1);
    v_val.push_back(val);
    adj.augment();
    e_ind.augment();
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      /*Return first node.*/
      return graph->node(i1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph->node(i2);
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
        if (i1+i2 == e.i1+e.i2) {
          if (i1*i2 == e.i1*e.i2) {
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
      *@return bool @a b, which signals comparison.
      *@post < satisfies trichotomy as defined for any mathematical order relation. 
      */
      if ((long int)graph < (long int)e.graph) return true;
      if (graph == e.graph) {
        if (i1+i2 < e.i1+e.i2) {
          return true;
        }
        if (i1+i2 == e.i1+e.i2) {
          if (i1*i2 < e.i1*e.i2) {
            return true;
          }
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
    size_type idx,i1,i2;
    
    Edge(const Graph* g, size_type i, size_type j1, size_type j2) : graph(g), idx(i), i1(j1), i2(j2) {}; 
  
    void switch_nodes() {
      size_type i;
      i = i2;
      i2 = i1;
      i1 = i;
    }
  };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return e_set.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return e_set[i];
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
    return adj.find(i1,i2)>0;
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
    size_type i,k, e_num = num_edges();
    k = adj.find(i1,i2);
    if (k > 0) {
      i=e_ind(i1,k);
      Edge e = edge(i);
      if (e.node1()!= a) e.switch_nodes();
      return e;
    }
    Edge e(this,e_num,i1,i2);
    update_edge_data(e);
    return e;
  }
  /*This is for preserving the RI of graph that
     * e_set is vector of edges. 
     * adj is the adjacency list,
     * e_ind is adjacency list with edge indices.
  */
  void update_edge_data(const Edge& e) {
    e_set.push_back(e);
    adj.add(e.i1,e.i2);
    e_ind.add(e.i1,e.idx);
    if (e.i2!=e.i1) {
      adj.add(e.i2,e.i1);
      e_ind.add(e.i2,e.idx);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    v_set.clear();
    e_set.clear();
    adj.clear();
    e_ind.clear();
  }
//<<<<<<< HEAD
  
//=======

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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
    
    /**Returns the Node being pointed by the node iterator.
    */
    value_type operator*() const {
      return Node(graph,idx);
    }
    /**Returns the Node being pointed by the node iterator.
    *@returns a NodeIterator.
    *@post @a this->index()+1 = @a result->index();
    */
    NodeIterator& operator++() {
      idx++;
      return *this;
    }
    /**
    *Checks if two iterators are equal
    *@param nIter: input iterator
    *@returns True or False.
    *
    *@post @a result == True if: 
         @a this->index() == @a n->index() and
         @a this->graph == @a nIter->graph.
    */
    bool operator==(const NodeIterator& nIter) const {
      return this->graph == nIter.graph && this->idx == nIter.idx;
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* graph;
    size_type idx;
    NodeIterator(const Graph* g,size_type i): graph(g), idx(i) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
  /** Create an iterator from the starting node.
  * @returns a node_iterator
  * @post @a result->index() == 0
  */
  node_iterator node_begin() const {
    return node_iterator(this,0);
  }
  /** Create an iterator at the point beyond the last node.
  * It does not point to a real node.
  * @returns a node_iterator
  */
  node_iterator node_end() const {
    return node_iterator(this,num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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
    
    /** Return edge pointed by iterator.
    * @returns an Edge.
    * @post @result.node1() == @a this->node()
    */
    value_type operator*() const {
      /*
      * e_idx is index of edge, obtained from our e_ind matrix,
      * by looking up our node n.index(), and number iterator over nbrs idx.
      * Create edge e from e_idx, switch nodes if necessary, and return e.
      */
      size_type e_idx;
      e_idx = graph->e_ind(n.index(),idx);
      Edge e;
      e = graph->edge(e_idx);
      if (e.node1() != n) e.switch_nodes();
      return e;
    }
    /** Increments the iterator.
    * @returns an IncidentIterator reference to itself.
    */
    IncidentIterator& operator++() {
      idx++;
      return *this;
    }
    /** Checks equality by checking if the graphs, node and indices are equal.
    * @param The iterator to compare with.
    * @returns bool.
    */
    bool operator==(const IncidentIterator& iter) const {
      if (graph == iter.graph && n == iter.n && idx == iter.idx) return true;
      return false;
    }

    /** Returns reference to node which we are using.
    */
    Node& node() const {
      return n;
    }
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graph;
    Node n;
    size_type idx;
    IncidentIterator(const Graph* g, Node node, size_type i): graph(g), n(node), idx(i) {}
    
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private totally_ordered<EdgeIterator> {
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

    /*Returns edge pointed by EdgeIterator.*/
    value_type operator*() const {
      return graph->edge(idx);
    }
    /**Returns EdgeIterator pointing at next edge.
    * @post @a result->idx == @a this->idx+1
    */
    EdgeIterator& operator++() {
      idx++;
      return *this;
    }

    /**Check EdgeIterator equality
    * @param eIter, the iterator to compare with.
    * @result bool which equals true if
        @a this->graph == @a eIter.graph and @a this->idx == @a eIter.idx 
    */
    bool operator==(const EdgeIterator& eIter) const {
      return graph == eIter.graph && idx == eIter.idx;
    }
    
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph;
    size_type idx;
    EdgeIterator(const Graph* g, size_type i): graph(g), idx(i) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  
  /* Returns an iterator at the start of the edge set.
  */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /* Returns an iterator at the end of the edge set.
  */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

//>>>>>>> CME212/master
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  /*@a v_set is the vertex set of the graph.*/
  std::vector<Point> v_set;

  /*@a v_val is the vertex value set of the graph.*/
  std::vector<V> v_val;
  
  struct Arr {
    std::vector<std::vector<size_type>> A;
    size_type operator()(size_type i, size_type j) const {
      return A[i][j];

    }
    size_type operator()(size_type i) const {
      /* Return degree, i.e. number of neighbors. */
      return A[i].size();
    }
    size_type find(size_type i, size_type j) const {
      size_type k = 0;
      for (auto iter = A[i].begin(); iter!= A[i].end(); iter++) {
        if (*iter == j) return k;
        k++;
      }
      return 0;
    }
    void clear() {
      A.clear();
    }
    void add(size_type i, size_type j) {
      A[i].push_back(j);
    }
    void augment() {
      std::vector<size_type> r(0,0);
      A.push_back(r);
    }
  };
  /*@a v_set is the edge set of the graph. Edges stored as two node indices.*/
  std::vector<Edge> e_set;
  Arr adj;
  Arr e_ind;
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
