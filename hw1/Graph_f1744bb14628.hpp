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

 public:

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

  // type of the node value
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): nodevec(), nodevalue(), edgevec(), adjcence(){
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

    /** Return this node's position. */
    const Point& position() const {
      return graph->nodevec[ind];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return ind;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    //return the value of the node
    node_value_type& value() {
      return graph->nodevalue[ind];
    }

    //return the value of the node, which could be modified
    const node_value_type& value() const {
      return graph->nodevalue[ind];
    }
    
    //return the degree of the node
    size_type degree() const {
      return graph->adjcence[ind].size();
    }

    //return an iterator pointing to the first incident edge 
    incident_iterator edge_begin() const {
      return IncidentIterator(graph, ind, 0);
    }

    //return an iterator pointing to the last incident edge 
    incident_iterator edge_end() const {
      return IncidentIterator(graph, ind, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph == n.graph and ind == n.ind);
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
      if(graph == n.graph){
          return ind < n.ind;
      }
      return graph < n.graph; 
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    graph_type* graph;

    // This node's index
    size_type ind;

    Node(const graph_type* graph_, size_type index)
        : graph(const_cast<graph_type*>(graph_)), ind(index) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type node_size() const {
    return nodevec.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return nodevec.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    nodevec.push_back(position);
    nodevalue.push_back(value);
    return Node(this, nodevec.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph == this;
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
 
  void setValue(size_type index, node_value_type value){
    nodevalue[index] = value;
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
      return Node(graph, graph->edgevec[ind][0]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph, graph->edgevec[ind][1]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph == e.graph and ind == e.ind); 
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(graph == e.graph){
          return ind < e.ind;
      }
      return graph < e.graph;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    graph_type* graph;

    // This edge's index
    size_type ind;

    Edge(const graph_type* graph_, size_type index)
        : graph(const_cast<graph_type*>(graph_)), ind(index) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edgevec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(adjcence[a.index()].size() < adjcence[b.index()].size()){
        for(size_type i = 0; i < adjcence[a.index()].size(); i++){
            if(adjcence[a.index()][i][0] == b.index()){
                return true;
            }
        }
        return false;
    }
    for(size_type i = 0; i < adjcence[b.index()].size(); i++){
        if(adjcence[b.index()][i][0] == a.index()){
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
    if(adjcence[a.index()].size() < adjcence[b.index()].size()){
        for(size_type i = 0; i < adjcence[a.index()].size(); i++){
            if(adjcence[a.index()][i][0] == b.index()){
                return Edge(this, adjcence[a.index()][i][1]);
            }
        }
    } 
    for(size_type i = 0; i < adjcence[b.index()].size(); i++){
        if(adjcence[b.index()][i][0] == a.index()){
            return Edge(this, adjcence[b.index()][i][1]);
        }
    }
    edgevec.push_back({a.index(), b.index()});
    if(a.index() >= adjcence.size()){
        for(size_type i = adjcence.size(); i <= a.index(); i++){
            adjcence.push_back({});
        }
    }
    adjcence[a.index()].push_back({b.index(), size_type(edgevec.size() - 1)});
    if(b.index() >= adjcence.size()){
        for(size_type i = adjcence.size(); i <= b.index(); i++){
            adjcence.push_back({});
        }
    }
    adjcence[b.index()].push_back({a.index(), size_type(edgevec.size() - 1)});
    return Edge(this, edgevec.size() - 1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodevec.clear();
    edgevec.clear();
    adjcence.clear();
  }

  //
  // Node Iterator
  //


  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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

    //Return the node pointed by the iterator
    Node operator*() const {
      return Node(graph, ind);      
    }
    
    //Move the index to the next point
    NodeIterator& operator++() {
      ind++;
      return *this;
    }

    //Test whether two iterators are equal
    bool operator==(const NodeIterator& n) const {
      return (graph == n.graph and ind == n.ind);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // Pointer back to the Graph container
    graph_type* graph;

    // This node's index
    size_type ind;

    NodeIterator(const graph_type* graph_, size_type index)
        : graph(const_cast<graph_type*>(graph_)), ind(index) {
    }

    
    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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

    // return the incident edge 
    Edge operator*() const {
      return Edge(graph, graph->adjcence[node_ind][edge_ind][1]);
    }

    // move the index to the next position
    IncidentIterator& operator++() {
      edge_ind++;
      return *this;
    }

    // test whether the two incidnet iterator are the same
    bool operator==(const IncidentIterator& e) const {
      return (graph == e.graph and node_ind == e.node_ind and edge_ind == e.edge_ind);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    // Pointer back to the Graph container
    graph_type* graph;

    // This node's index
    size_type node_ind;

    // Incident edge's index
    size_type edge_ind;

    IncidentIterator(const graph_type* graph_, size_type node_index, size_type edge_index)
        : graph(const_cast<graph_type*>(graph_)), node_ind(node_index), edge_ind(edge_index) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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

    // return the edge pointed by the iterator
    Edge operator*() const {
      return Edge(graph, ind);
    }

    // move the index to nect position
    EdgeIterator& operator++() {
      ind++;
      return *this;
    }

    // test whether two iterators are the same
    bool operator==(const EdgeIterator& e) const {
      return (graph == e.graph and ind == e.ind);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
    // Pointer back to the Graph container
    graph_type* graph;

    // This edge's index
    size_type ind;

    EdgeIterator(const graph_type* graph_, size_type index)
        : graph(const_cast<graph_type*>(graph_)), ind(index) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  //
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  //
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:
  
  std::vector<Point> nodevec;
  std::vector<node_value_type> nodevalue;
  std::vector<std::vector<size_type>> edgevec;
  std::vector<std::vector<std::vector<size_type>>> adjcence;

};

#endif // CME212_GRAPH_HPP
