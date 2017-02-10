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
template <typename V, typename E = double>
class Graph {
private:
  // Added private variables to end of Graph
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

  /** Type of values within each node */
  using node_value_type = V;
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    edge_count = 0;
    node_count = 0;
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
   * @invariant graph must be a valid object t
   * @invariant the node id cannot be larger than the graph size
   */
  class Node : private totally_ordered<Node> {
  public:
    // Construct an invalid node.
    Node() {}

    /** Construct a valid node with given graph and id */
    Node(const Graph* g, size_type id) {
      graph = const_cast<Graph*>(g);
      uid = id;
      assert(graph != nullptr);
      assert(uid < g->size());
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph->nodes[uid].first;
    }
    
    /** Return this node's position as a reference. */
    Point& position() {
      return graph->nodes[uid].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid;
    }

    /** Return the nodes value */
    node_value_type& value(){
      return graph->nodes[uid].second;
    }

    /** Return node value as a const */
    const node_value_type& value() const{
      return graph->nodes[uid].second;
    }

    /** Return degree of node */
    size_type degree() const{
      return graph->adj_list[uid].size();
    }

    /** Return iterator that points to first edge of a node */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph, uid, 0); 
    }

    /** Return iterator that points past last edge of a node */
    incident_iterator edge_end() const{
      return IncidentIterator(graph, uid, degree());
    }

    /** Test whether this node and @a n are equal.
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return n.graph == graph && n.uid == uid;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      return (uid < n.uid && graph ==  n.graph) || (graph < n.graph);
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer to a main graph class holding pointers
    Graph* graph;

    /** Positive integer representing unique id of node */
    size_type uid; 
  };

  /** Return the number of nodes in the graph.
   * Complexity: O(1).
   */
  size_type size() const {
    return node_count;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, with a specified value and
    return the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = 
      node_value_type()) {
    nodes.push_back(std::pair<Point, node_value_type>(position,value));

    // Add node to adjacency with initial empty vector
    adj_list.push_back(std::vector<std::pair<size_type, 
          edge_value_type>>());
    node_count++;  
    return Node(this, node_count-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
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
   * @invariant A node cannot have an edge to itself
   * @invariant graph != nullptr
   */
  class Edge : private totally_ordered<Edge>{
  public:

    /** Construct an invalid Edge. */
    Edge() {}

    /** Construct valid edge with a given graph and two nodes
     */
    Edge(const Graph* g, size_type one, size_type two) {
      graph = const_cast<Graph*>(g);
      n1 = one;
      n2 = two;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph, n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return Node(graph, n2);
    }

    edge_value_type& value(){
      if(n1 > n2){
        size_type temp = n2;
        n2 = n1;
        n1 = temp;
      }
      for(unsigned i = 0; i < graph->adj_list[n1].size(); ++i){
        if(graph->adj_list[n1][i].first == n2){
          return graph->adj_list[n1][i].second;
        }
      }
      return graph->adj_list[0][0].second;
    }
 
    const edge_value_type& value() const{
      for(unsigned i = 0; i < graph->adj_list[n1].size(); ++i){
        if(graph->adj_list[n1][i].first == n2){
          return graph->adj_list[n1][i].second;
        }
      }
      return edge_value_type();
    }
    /** Test whether this edge and @a e are equal.
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return e.graph == graph && ((e.n1 == e.n1 && e.n2 == n2) || 
        (e.n1 == n2 && e.n2 == n1)); 
    }

    /** Test whether this edge is less than @a e in a global order.
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (n1 < e.n1 && graph == e.graph) || (graph < e.graph);
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer to graph class containing edge data
    Graph* graph;

    // Unique id of edge object to gather data
    size_type n1;
    size_type n2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    EdgeIterator ei = edge_begin();
    while(i > 0){
      ++ei;
      --i;
    }
    return *ei;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for(auto  p : adj_list[a.uid]){
      if(p.first == b.uid){ 
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
   *       Else,                        new num_edges() == old num_edges()+1
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, 
      const edge_value_type& value = edge_value_type()) {
    // Check if graph already contains edge
    for(auto p : adj_list[a.index()]){
      if(p.first == b.index()){
        return Edge(this, a.index(), b.index());
      }
    }

    // Add edge to both nodes of adj list
    adj_list[a.index()].push_back(std::pair<size_type, edge_value_type>(b.index(), value));
    adj_list[b.index()].push_back(std::pair<size_type, edge_value_type>(a.index(), value));
    edge_count++;
    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    adj_list.clear();
    node_count = 0;
    edge_count = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. 
   * @invariant uid <= g->size()
   */
  class NodeIterator : private totally_ordered<NodeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Construct valid NodeIterator with graph pointer and id */
    NodeIterator(const Graph* g,  size_type id){
      uid = id;
      graph = const_cast<Graph*>(g);
      assert(uid <= graph->size());
    }

    /** Returns the Node object that the iterator points to */
    Node operator*() const{
      return Node(graph, uid);
    }

    /** Increments NodeIterator by one and returns itself*/
    NodeIterator& operator++(){
      uid++;
      return *this; 
    }

    /** Determines whether two node iterators are equal */
    bool operator==(const NodeIterator& iter) const{
      return graph == iter.graph && uid == iter.uid;
    }

  private:
    friend class Graph;
    Graph* graph;

    // Id for which node it refers to
    size_type uid; 
  };

  /** Returns the iterator that points to the first node in the graph */
  node_iterator node_begin() const{
    return node_iterator(this, 0);
  }

  /** Return iterator that points to one after the last node in the graph */
  node_iterator node_end() const{
    return node_iterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   * @invariant node_id and the edge_id are both valid
   */
  class IncidentIterator : private totally_ordered<IncidentIterator> {

  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** Construct a valid IncidentIterator. */
    IncidentIterator(const Graph* g, size_type node, size_type edge){
      graph = const_cast<Graph*>(g);
      node_id = node;
      edge_id = edge;
      assert(node_id <= graph->size());
      assert(edge_id <= graph->adj_list[node_id].size());
    }

    /** Return the edge that the IncidentIterator poitns to */
    Edge operator*() const{
      size_type return_index = graph->adj_list[node_id][edge_id].first;
      return Edge(graph, node_id, return_index);
    }

    /** Increment the iterator by one and return a reference to itself */
    IncidentIterator& operator++(){
      edge_id++;
      return *this;
    }

    /** Test whether two incident iterators are equal */
    bool operator==(const IncidentIterator& iter) const{
      return graph == iter.graph && node_id == iter.node_id &&
        edge_id == iter.edge_id;
    }

  private:
    friend class Graph;
    Graph* graph;

    // UID of node in question
    size_type node_id;
    
    // Index of edge in adj list
    size_type edge_id;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator.
   * @invariant node index and edge index are valid
   */
  class EdgeIterator : private totally_ordered<EdgeIterator> {

  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}
    
    /** Construct an EdgeIterator using provided node and edge index as well
      * as graph pointer
      */
    EdgeIterator(const Graph* g, size_type node, size_type edge){
      graph = const_cast<Graph*>(g);
      node_id = node;
      edge_id = edge;
      assert(node_id <= graph->size());
      assert(edge_id <= graph->adj_list[node_id].size());
    }

    // Reference back to edge object
    Edge operator*() const{
      return Edge(graph, node_id, graph->adj_list[node_id][edge_id].first);
    }

    /** Increments EdgeIterator by one. If there are no more edge for a 
      * given node, it will increment to next node. Only counts each edge
      * once.
      */
    EdgeIterator& operator++(){
      // Increment by one initially
      if(graph->adj_list[node_id].size() > edge_id+1){
        ++edge_id;
      } else{
        ++node_id;
        edge_id = 0;
      }
      
      /** Skip over times when node is less than other node of edge
          to prevent duplicates */
      while(node_id != graph->num_nodes() && 
          graph->adj_list[node_id][edge_id].first < node_id){ 
        if( graph->adj_list[node_id].size() > edge_id+1){
          ++edge_id;
        } else{
          ++node_id;
          edge_id = 0;
        }
      }
      return *this;
    }

    /** Test for equality between two EdgeIterators */
    bool operator==(const EdgeIterator& iter) const{
      return graph == iter.graph && node_id == iter.node_id 
        && edge_id == iter.edge_id;
    }

  private:
    friend class Graph;
    Graph* graph;

    // Node index for one side of edge
    size_type node_id;

    // Edge index of adjancency list
    size_type edge_id;
  };

  /** Return edge iterator that points to first edge */
  edge_iterator edge_begin() const{
    return edge_iterator(this, 0, 0);
  }

  /** Return edge iterator that points to one after the last edge */
  edge_iterator edge_end() const{
    return edge_iterator(this, num_nodes(), 0);
  }

private:
  // Carries position and value data of nodes
  std::vector<std::pair<Point, node_value_type>> nodes;

  // Adjacency list as an alternative representation of edges
  std::vector<std::vector<std::pair<size_type, edge_value_type>>> adj_list;

  // Graph node count
  size_type node_count;

  // Graph edge count 
  size_type edge_count;
};

#endif // CME212_GRAPH_HPP
