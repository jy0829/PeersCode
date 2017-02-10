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
// Defining template
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Define an internal vector to store the nodes
  

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;
  using node_value_type = V;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : internal_nodes(), internal_edges(), adjacency() {
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
    Node() {
      // HW0: YOUR CODE HERE
    } // Construction of invalid node

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->internal_nodes[uid_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    // node_value_type
    node_value_type& value() {
      return graph_->internal_nodes[uid_].second;
    }

    const node_value_type& value() const {
      return graph_->internal_nodes[uid_].second;
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->graph_ == n.graph_)
        if (this->uid_ == n.uid_)
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
      // HW0: YOUR CODE HERE
      // Comparing nodes in same graph 
      if (this->graph_ == n.graph_)
        if (this->uid_ < n.uid_)
          return true;

      // Comparing nodes in different graphs
      return (this->graph_ < n.graph_);
    }

    /** Return the degree of this node
    */
    size_type degree() const {
      return graph_->adjacency[uid_].size();
    }

     /** Return an incident iterator to the first edge incident on a node
    */
    IncidentIterator edge_begin() const {
      return IncidentIterator(graph_,uid_,0);
    }

    /** Return an incident iterator that points past the last edge incident on 
    a node
    */
    IncidentIterator edge_end() const {
      return IncidentIterator(graph_,uid_,degree());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_; // Pointer to Graph class
    size_type uid_; // ID number    
    Node(const Graph* graph, size_type uid)
       : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    } // Valid Constructor

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return internal_nodes.size();
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
  // Default constructor of node_value_type is used to assign a default value
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    // HW0: YOUR CODE HERE
    internal_nodes.push_back(std::make_pair(position,node_value_type())); // Adding nodes to the internal_nodes vector
    adjacency.push_back(std::vector<std::pair<size_type, size_type>> ()); // Empty adjacency vector 

    return Node(this,internal_nodes.size()-1);   // Return node to the user
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_ == this)
      if (n.uid_ < n.graph_-> num_nodes())
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
    // HW0: YOUR CODE HERE
    return Node(this,i);  // Node with index i
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_,this->node1_);      // Node 1
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_,this->node2_);     // Node 2
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {  
      if (std::min(node1(),node2()) == std::min(e.node1(),e.node2()))
        if (std::max(node1(),node2()) == std::max(e.node1(),e.node2()))
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Comparing edges in same graph
      if (this->graph_ == e.graph_)
        if (this->uid_ < e.uid_)
          return true;

      // Comparing edges between different graphs
      if (this->graph_ < e.graph_)
        return true;

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
    size_type uid_,node1_,node2_; // edge id, node ids
    Edge(const Graph* graph, size_type uid, size_type node1, size_type node2)
       : graph_(const_cast<Graph*>(graph)), uid_(uid),node1_(node1),node2_(node2) {
    } // Valid Constructor

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return internal_edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this,i,internal_edges[i].first,internal_edges[i].second);       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    auto size_ad = adjacency[a.uid_].size();
    for (size_type i = 0; i<size_ad; i++){
      if (adjacency[a.uid_][i].first == b.uid_)
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
    // Return existing edge
    for (size_type i=0; i<adjacency[a.uid_].size();++i)
      if (adjacency[a.uid_][i].first == b.uid_)
        return Edge(this,adjacency[a.uid_][i].second,a.uid_,b.uid_);       

    // add new edge
    internal_edges.push_back(std::make_pair(a.uid_,b.uid_));
    adjacency[a.uid_].push_back(std::make_pair(b.uid_,internal_edges.size()-1));
    adjacency[b.uid_].push_back(std::make_pair(a.uid_,internal_edges.size()-1));
    return Edge(this,internal_edges.size()-1,a.uid_,b.uid_);          
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    internal_nodes.clear();
    internal_edges.clear();
    adjacency.clear();
    return;
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

    /** Deference node iterator 
    *@pre NodeIterator points to a valid Node object of this graph
    *@post returns Node object, such that
    *Node.graph_ = NodeIterator.graph_;
    *Node.uid_ = NodeIterator.uid_;
    */
    Node operator*() const {
      // Dereferencing operator
      return Node(graph_,uid_);
    } 

    /** Increment operator
    *@pre NodeIterator points to a valid Node object of this graph
    *@post returns NodeIterator, such that
    *NodeIterator_return.graph_ = NodeIterator.graph_;
    *NodeIterator_return.uid_ = NodeIterator.uid_ +1;
    */
    NodeIterator& operator ++(){
      // Increment operator
      uid_++;
      return *this;
    }

    /** Test whether this NodeIterator is equal to @a node_iterator in a global order.
    */    
    bool operator ==(const NodeIterator& node_iterator) const {
      return ((graph_ == node_iterator.graph_)
        && (uid_ == node_iterator.uid_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_; // Pointer to Graph class
    size_type uid_; // ID number    
    NodeIterator(const Graph* graph, size_type uid)
       : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    } // Valid Constructor

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns a NodeIterator to the start of the container containing nodes
  */
  NodeIterator node_begin () const {
    return NodeIterator(this,0);
  }

  /** Returns a NodeIterator to one past the end of the container
  * containing nodes
  */
  NodeIterator node_end () const {
    return NodeIterator(this,num_nodes());
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

    /** Deference IncidentIterator
    *@pre IncidentIterator points to a valid Edge object of this graph
    *@post returns Edge object
    */
    Edge operator *() const {
      return Edge(graph_,graph_->adjacency[node_id_][uid_].second, 
        node_id_, graph_->adjacency[node_id_][uid_].first);
    }

      /** Increment operator
    *@pre IncidentIterator points to a valid Edge object of this graph
    *@post returns IncidentIterator
    */
    IncidentIterator& operator ++() {
      // Increment
      ++uid_;
      return *this;
    }

    /** Test whether this IncidentIterator is equal to @a iit in a global order.
    */
    bool operator ==( const IncidentIterator & iit ) const {
      return ((graph_ == iit.graph_) && (node_id_ == iit.node_id_)
        && (uid_ == iit.uid_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_id_, uid_; // edge id, node ids
    IncidentIterator(const Graph* graph, size_type node_id, size_type uid)
       : graph_(const_cast<Graph*>(graph)), node_id_(node_id), uid_(uid) {
    } // Valid Constructor
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
    EdgeIterator()  {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Deference node iterator 
    *@pre EdgeIterator points to a valid Edge object of this graph
    *@post returns edge object
    */
    Edge operator*() const {
      // Dereferencing operator
      return Edge(graph_,uid_,graph_->internal_edges[uid_].first,graph_->internal_edges[uid_].second);
    } 

       /** Increment operator
    *@pre EdgeIterator points to a valid Edge object of this graph
    *@post returns EdgeIterator
    */
    EdgeIterator& operator ++(){
      // Increment operator
      uid_++;
      return *this;
    }

    /** Test whether this EdgeIterator is equal to @a edge_iterator in a global order.
    */
    bool operator ==(const EdgeIterator& edge_iterator) const {
      return ((graph_ == edge_iterator.graph_)
        && (uid_ == edge_iterator.uid_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_; // Pointer to Graph class
    size_type uid_; // ID number    
    EdgeIterator(const Graph* graph, size_type uid)
       : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    } // Valid Constructor

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

   /** Returns a EdgeIterator to the start of the container containing edges
  */

  EdgeIterator edge_begin () const {
    return EdgeIterator(this,0);
  }

  /** Returns a EdgeIterator to one past the end of the container
  * containing edges
  */
  EdgeIterator edge_end () const {
    return EdgeIterator(this,num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<std::pair<Point,node_value_type>> internal_nodes;
  std::vector<std::pair<size_type,size_type>> internal_edges;

  // For efficient comparison in has_edge() method, we store the adjacency matrix as well
  std::vector<std::vector<std::pair<size_type, size_type>>> adjacency;

};

#endif // CME212_GRAPH_HPP
