#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <utility>
#include <vector>

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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** Interpretation of nodes contained in graph. */
  using node_value_type = V;

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

  struct vertex {
    Point p;
    node_value_type value;
  };

  /** Vector of Points that represent each node in graph */
  std::vector<vertex> nodes;
  /** Vector containing indices of nodes that form each edge */
  std::vector<std::pair<size_type, size_type>> edges;
  /** Map containing strings of adjacent nodes */
  std::map<size_type, std::vector<size_type>> connections;

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
    Node() { 
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes[node_index_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return @a value of node at @a index() */
    node_value_type& value() {
      return graph_->nodes[index()].value; 
    }
    /** Return @a value of node at @a index() in const graph */
    const node_value_type& value() const {
      return graph_->nodes[index()].value; 
    }
    /** Return number of nodes connected to current node (the degree) */
    size_type degree() const {
      return graph_->connections[index()].size(); 
    }
    /** Return IncidentIterator that starts at the first adjacent node to current node */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index(), 0);
    }
    /** Return IncidentIterator that signifies one position past last valid edge */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index(), degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ( (this->node_index_ == n.node_index_) && (n.graph_ == this->graph_) ) {
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
      if ((this->graph_ == n.graph_) && (this->node_index_ < n.node_index_)) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  
    /** Pointer back to the Graph container */
    graph_type *graph_;
    /** This node's index in Graph container */
    size_type node_index_;

    Node(const graph_type* g, size_type nn) 
      : graph_(const_cast<graph_type*>(g)), node_index_(nn) {

    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
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
    vertex vert;
    vert.p = position;
    vert.value = val;
    //nodes.push_back(position);
    nodes.push_back(vert);
    return Node(this, nodes.size() - 1);    
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /** Using the node invariant property node(i) => node.node_index_ = i */
    size_type idx = n.node_index_;
    /** Check if graph has at least this many indices and check node equality */
    if ( (num_nodes() > idx) && (node(idx) == n) ) {
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() { 

    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->graph_, uid_a);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(this->graph_, uid_b);   
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_) {
        return false;
      }
      Node a = node1();
      Node b = node2();

      Node e1 = e.node1();
      Node e2 = e.node2();
      if ( (a == e1 && b == e2) || (a == e2 && b == e1) ) {
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
      if (  (this->graph_ == e.graph_) && (this->uid_a < e.uid_a) && (this->uid_b < e.uid_b) ) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type *graph_;
    size_type uid_a;
    size_type uid_b;

    Edge(const graph_type* g, size_type na, size_type nb) 
      : graph_(const_cast<graph_type*>(g)), 
        uid_a(na),
        uid_b(nb) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, edges[i].first, edges[i].second);    
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // iterate through edges and check if there exists an edge
    // that connects @a a and @a b.
    /** Use Connections map for complexity: O(num_adjacent_nodes) */
    if (connections.find(a.node_index_) == connections.end() 
      || connections.find(b.node_index_) == connections.end() ) {
      return false;
    }

    std::vector<size_type> a_adj = connections.at(a.node_index_);
    std::vector<size_type> b_adj = connections.at(b.node_index_);

    /** Search for a in b's list and vica-versa */
    if(std::find(a_adj.begin(), a_adj.end(), b.node_index_) != a_adj.end()
      && std::find(b_adj.begin(), b_adj.end(), a.node_index_) != b_adj.end() ) {
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
    if ( has_edge(a,b) ) {
      return Edge(this, a.node_index_, b.node_index_);
    }
    std::pair<size_type, size_type> ab = std::make_pair(a.node_index_, b.node_index_);
    edges.push_back(ab);
    connections[a.node_index_].push_back(b.node_index_);
    connections[b.node_index_].push_back(a.node_index_);
    
    return edge(edges.size() - 1);  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges.clear();
    nodes.clear();
    connections.clear();
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
    NodeIterator(const graph_type *g, size_type x) : graph_(const_cast<graph_type*>(g)), 
                                                     index_(x) {

    }

    /** Return Node that belongs to graph @a graph_ and has index @a index_.
     *
     * Must be a valid node.
    */
    Node operator*() const {
      return Node(graph_, index_);
    }

    /** Increment node index and return pointer to iterator.
     *
     * Nodes must be totally ordered.
    */
    NodeIterator& operator++() {
      this->index_++;
      return *this;
    }

    /** Return true if two node iterators are equal, false otherwise.
     *
     * Nodes belonging to same graph must be totally ordered.
    */
    bool operator==(const NodeIterator& iter) const {
      if(this->graph_ == iter.graph_ && index_ == iter.index_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    /** Pointer back to the Graph container */
    graph_type* graph_;
    /** This node's index in Graph container */
    size_type index_;
  };

  /** Return NodeIterator pointing to the start of node container.
    *
    * Node container must be nonempty.
  */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return NodeIterator that points to one position past last valid Node
    *
    * This node iterator must never be dereferenced.
  */
  node_iterator node_end() const {
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
    IncidentIterator(const graph_type* g, size_type i, size_type j)
      : graph_(const_cast<graph_type*>(g)), 
        central_node(i),
        adj_node(j) {

    }

    /** Return valid Edge from graph.
     *
     *  Returned edge connects node located at index @a central_node in
     *  adjacency matrix to the node with index stored at @a adj_node in 
     *  adjacency matrix.
     */
    Edge operator*() const {
      return Edge(graph_, central_node, graph_->connections[central_node][adj_node]);
    }

    /** Increment the index of the adjacent node to node at @a central_node
     *  and return IncidentIterator.
     *
     *  adj_node must always be valid when incremented. So it cannot equal
     *  the edge_end() that returns an IncidentIterator.
     *
     */
    IncidentIterator& operator++() {
      adj_node++;
      return *this;
    }

    /** IncidentIterator equality comparator function. 
     *  
     *  IncidentIterators are only equal if their graphs are equal,
     *  their @a central_nodes are equal, and their @a adj_nodes are equal.
     *
     *  Returns boolean value.
     */
    bool operator==(const IncidentIterator& iter) const {
      if(this->graph_ == iter.graph_ 
        && central_node == iter.central_node
        && adj_node == iter.adj_node) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    /** Pointer back to the Graph container */
    graph_type *graph_;
    /**  central and adjacent nodes */
    size_type central_node;
    size_type adj_node;
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
    EdgeIterator(const graph_type *g, size_type x) : graph_(const_cast<graph_type*>(g)), 
                                                                              index_(x) {

    }

    /** Dereference EdgeIterator and return valid Edge. 
     *  
     *  Returns edge at index @a index_. There is no guarantee
     *  on the ordering of the nodes for a given edge.
     *
     *  EX: if Edge(a,b) is returned then node1() and node2()
     *  will each return one or the other. 
     */
    Edge operator*() const {
      return graph_->edge(index_);
    }

    /** Increment and return EdgeIterator to point at the next valid edge. 
     *  
     *  The next valid edge must be valid to be dereferenced. 
     */
    EdgeIterator& operator++() {
      index_++;
      return *this;
    }

    /** Equality comparator for EdgeIterators. 
     *  
     *  EdgeIteratos are equal if their @a graphs_ are the same and 
     *  their @a index_s are equal. 
     */
    bool operator==(const EdgeIterator& ej) const {
      if (graph_ == ej.graph_ && index_ == ej.index_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    graph_type *graph_;
    size_type index_;
  };

 /** Returns EdgeIterator that points to the start of a container of edges. 
  *  
  *  Edge at this index must be valid.
  */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
    
 /** Returns EdgeIterator that points to one index past the last valid edge.
  *  
  *  Returned edge must not be dereferenced. This edge is not valid
  */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->num_edges());
  }

 private:


};

#endif // CME212_GRAPH_HPP
