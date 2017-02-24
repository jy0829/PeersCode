#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <set>
#include <utility>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <cmath>

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

  struct NodeData;
  
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

  /** Type of node value */
  using node_value_type = V;
  using edge_value_type = E;
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
  Graph():next_id_(0){
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
    Node(): id_(0), graph_(nullptr) {
      // HW0: YOUR CODE HERE
    }

    int fetch() const{
      size_t low = 0;
      size_t high = graph_ -> nodes.size() - 1;
      size_t mid = (id_ <= high? id_ : (low+high)/2);
      while (low <= high){
        if ( graph_ -> nodes[mid].id_ < id_){
          low = mid + 1;
          mid = (low + high)/2;
        } else if (graph_ -> nodes[mid].id_ > id_){
          high = mid - 1;
          mid = (low + high)/2;
        } else return mid;
      }
      return -1;
    }
          
      
    /** Return this node's position. */
    Point& position() const {
      // HW0: YOUR CODE HERE
      #ifdef ENABLE_DEBUG
      if (graph_ == nullptr || id_ >= graph_ -> nodes.size()) {
        std::cerr << "Invalid Node" << std::endl;
      }
      
      #endif
      return graph_ -> nodes[fetch()].position_;
      
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch();
    }
    
    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR  
      return (id_ == n.id_ && graph_ == n.graph_);
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
      
      return (id_ < n.id_);
    }

    node_value_type& value(){
      return graph_ -> nodes[fetch()].value_;
    }

    const node_value_type& value() const{
      return graph_ -> nodes[fetch()].value_;
    }

    size_type degree() const{
      return graph_ -> incident_map[id_].size();
    }

    incident_iterator edge_begin() const{
      return IncidentIterator(graph_ -> incident_map[id_].begin());
    }
    
    incident_iterator edge_end() const{
      return IncidentIterator(graph_ -> incident_map[id_].end());
    }

    
        
    

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class Graph::NodeIterator;
    friend class Graph::Edge;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type id_;
    Graph* graph_;


    // private valid constructor
    Node (const Graph* graph, size_type id):
      id_(id), graph_(const_cast<Graph*> (graph)){
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& nv = node_value_type()) {
    // HW0: YOUR CODE HERE
    Node newNode(this, next_id_);
    nodes.push_back(NodeData(next_id_, position, nv));
    ++ next_id_;
    return newNode;   
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.fetch() == -1) return false;
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, nodes[i].id_);
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
    Edge():graph_(nullptr),first(0),second(0) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      Node n1(graph_, first);
      //      return graph_ -> node(n1.fetch());
      return n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      Node n2(graph_, second);
      //      return graph_ -> node(n2.fetch());
      return n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ != e.graph_) return false; 
      size_type this1=first, this2=second, e1=e.first, e2=e.second;
      if (this1 > this2) std::swap(this1, this2);
      if (e1    > e2   ) std::swap(e1   , e2   );
      return (this1 == e1 && this2 == e2) ;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type this1=first, this2=second, e1=e.first, e2=e.second;
      if (this1 > this2) std::swap(this1, this2);
      if (e1    > e2   ) std::swap(e1   , e2   );
      
      if (this1 < e1) return true;
      if (this1 > e1) return false;
      if (this2 < e2) return true;
      return false;
    }

    /*
     * returns the current length
     */
    double length() const{
      auto n1 = node1();
      auto n2 = node2();
      return norm(n1.position() - n2.position());
    }


    /* return custom value */
    edge_value_type& value(){
      return graph_ -> edge_values_[*this];
    }

    const edge_value_type& value() const{
      return graph_ -> edge_values_[*this];
    }
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    std::size_t first;
    std::size_t second;
    
    // valid constructor
    Edge (const Graph* graph, size_type first_arg, size_type second_arg):
      graph_(const_cast<Graph*> (graph)),first(first_arg), second(second_arg){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
      auto it = edges.begin();
      std::advance(it, i);
      return *it;
      
    }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if ( edges.find(Edge(this, a.id_ , b.id_))== edges.end()){
      return false;
    } else {
      return true;
    }
      
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
    Edge newedge(this, a.id_, b.id_);
    Edge newedge2(this, b.id_, a.id_);
    if ( edges.find(newedge) == edges.end()){
      edges.insert (newedge);
      edge_values_[newedge] = edge_value_type();
      if (incident_map.find(a.id_) == incident_map.end()){
        std::vector<Edge> newvec;
        newvec.push_back(newedge);
        incident_map[a.id_]=newvec;
      } else {
        incident_map[a.id_].push_back(newedge);
      }
      if (incident_map.find(b.id_) == incident_map.end()){
        std::vector<Edge> newvec;
        newvec.push_back(newedge2);
        incident_map[b.id_] = newvec;
      } else {
        incident_map[b.id_].push_back(newedge2);
      }
    }
    return newedge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    // HW0: YOUR CODE HERE
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<Graph::NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator():it_(nullptr) {
    }

    NodeIterator(const NodeIterator& other):it_(other.it_),graph_(other.graph_){} 

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    value_type operator*() const{
      return Node(graph_, it_ -> id_);
    }
    
    // NodeIterator& operator++()
    NodeIterator& operator++(){
      (it_)++;
      return *this;
    }
    
    // bool operator==(const NodeIterator&) const
    bool operator== (const node_iterator& other) const{
      return (it_ == other.it_);
    }
    
   private:
    typename std::vector<NodeData>::iterator it_;
    graph_type *graph_;
    NodeIterator (typename std::vector<NodeData>::iterator it, const graph_type* graph):it_(it), graph_(const_cast<graph_type*>(graph)){
    }
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  node_iterator node_begin() const {
    return NodeIterator(const_cast<Graph*>(this)->nodes.begin(), this);
  }
  // node_iterator node_end() const

  node_iterator node_end() const{
    return NodeIterator(const_cast<Graph*>(this)->nodes.end(), this);
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
  
    IncidentIterator(){
    }

    IncidentIterator(typename std::vector<Edge>::iterator i): it(i){
    }
    
    value_type operator*() const{
      return *it;
    }

    IncidentIterator& operator++(){
      ++it;
      return *this;
    }

    bool operator== (const IncidentIterator & iit) const{
      return (it == iit.it);
    }
    
    
    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    private:
    typename std::vector<Edge>::iterator it;
 
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator :private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator():it_(nullptr) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const{
      return *it_;
    }
    
      
    // EdgeIterator& operator++()
    EdgeIterator& operator++(){
      ++it_;
      return *this;
    }
      
    // bool operator==(const EdgeIterator&) const
    bool operator==(const EdgeIterator& e) const{
      return (it_ == e.it_);
    }
    

   private:
    friend class Graph;
    typename std::set<Edge>::iterator it_;
    // private constructor
    EdgeIterator(typename std::set<Edge>::iterator it): it_(it){
    }
    
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  edge_iterator edge_begin() const{
     return edge_iterator(edges.begin());
  }
  // edge_iterator edge_end() const
  edge_iterator edge_end() const{
    return edge_iterator(edges.end());
  }

  size_type remove_edge(const Edge& e){
    // remove edge @a e from edges, incident_map, edge_values
    if (!has_edge(e.node1(), e.node2())) return 0;
    edges.erase(e);
    for (auto n : std::vector<std::size_t>{e.first, e.second}){
      auto & vec = incident_map[n];
      auto it = std::find(vec.begin(), vec.end(), e);
      std::iter_swap(it, vec.rbegin());
      vec.pop_back();
    }
    edge_values_.erase(e);
    return 1;
  }

  size_type remove_edge(const Node& n1, const Node& n2){
    // remove edge @a e from edges, incident_map, edge_values
    if (!has_edge(n1, n2)) return 0;
    Edge e(this, n1.id_, n2.id_);
    edges.erase(e);
    for (auto n : std::vector<std::size_t>{e.first, e.second}){
      auto & vec = incident_map[n];
      auto it = std::find(vec.begin(), vec.end(), e);
      std::iter_swap(it, vec.rbegin());
      vec.pop_back();
    }
    edge_values_.erase(e);
    return 1;
  }

  edge_iterator remove_edge(edge_iterator e_it){
    // remove edge @a e from edges, incident_map, edge_values
    Edge e = *e_it;
    e_it = edges.erase(e);
    for (auto n : std::vector<std::size_t>{e.first, e.second}){
      auto & vec = incident_map[n];
      auto it = std::find(vec.begin(), vec.end(), e);
      std::iter_swap(it, vec.rbegin());
      vec.pop_back();
    }
    edge_values_.erase(e);
    return e_it;
  }

  size_type remove_node(const Node& n){
    // remove node and connected edges
    if (!has_node(n)) return 0;   
    for (auto ptr = incident_map[n.id_].rbegin(); ptr != incident_map[n.id_].rend(); ++ptr){
      remove_edge(*ptr);
    }
    incident_map.erase(n.id_);
    nodes.erase(nodes.begin() + n.fetch());
    return 1;
  }

  node_iterator remove_node(node_iterator n_it){
    // remove node and connected edges
    Node n = *n_it;
    if (!has_node(n)) return n_it;   
    for (auto ptr = incident_map[n.id_].rbegin(); ptr != incident_map[n.id_].rend(); ++ptr){
      remove_edge(*ptr);
    }
    incident_map.erase(n.id_);
    auto it = nodes.erase(nodes.begin() + n.fetch());
    return NodeIterator(it, this);
  }                                           
    
  

  
 private:
  
  struct EdgeHasher{
    std::size_t operator()(Edge e) const{
      return std::hash<std::size_t>()(e.first) ^ std::hash<std::size_t>()(e.second);
    }
  };

  struct NodeData{
    size_t id_;
    Point position_;
    node_value_type value_;
    NodeData(size_t id, Point position, node_value_type value):id_(id), position_(position), value_(value){}
  };
  std::vector<NodeData> nodes;
  std::set<Edge> edges;
  std::unordered_map<unsigned long, std::vector<Edge> > incident_map;
  std::map<Edge, edge_value_type> edge_values_;
  size_t next_id_;
  //    std::unordered_map<Edge, edge_value_type, EdgeHasher> edge_values_;
  
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    

};

#endif // CME212_GRAPH_HPP
