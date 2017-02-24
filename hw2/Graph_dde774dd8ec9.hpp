 #ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>


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
  /** Allow user-specified node value. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Allow user-specified node value. */
  using edge_value_type = E;

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
    // graph_edge_index_ = 0;
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
      graph_ = NULL;
      node_index_ = 0;
    }

    /** Return this node's position. */
    Point& position() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> nodes[node_index_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_ -> nodes[node_index_].idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return this node's value*/
    node_value_type& value(){
      assert(graph_);
      return graph_ -> nodes[node_index_].value_;
    }

    /** Return this const node's const value*/
    const node_value_type& value() const{
      assert(graph_);
      return graph_ -> nodes[node_index_].value_;     
    }

    /** Return this node's degree*/
    size_type degree() const {
      assert(graph_);
      return graph_ -> nodes[node_index_].valid_adjacents_.size();
    }

    /** Return this node's begin incident iterator*/
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, node_index_, graph_ -> nodes[node_index_].valid_adjacents_.begin());
    }

    /** Return this node's end incident iterator*/
    incident_iterator edge_end() const {
      // assert(node_index_ < graph_->num_nodes());
      return incident_iterator(graph_, node_index_, graph_ -> nodes[node_index_].valid_adjacents_.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // return true when the nodes have the same graph and the same index;
      if (graph_ && n.graph_) {
        return graph_ == n.graph_ && node_index_ == n.node_index_;
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

      return std::tie(graph_, node_index_) < std::tie(n.graph_, n.node_index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type node_index_;
    Node(const graph_type* graph, size_type node_index) {
      graph_ = const_cast<graph_type*>(graph);
      node_index_ = node_index;
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return valid_nodes.size();
    // return 0;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {

    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[] node_value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE    
    size_type node_index = nodes.size();
    size_type idx = valid_nodes.size();
    valid_nodes.push_back(node_index);
    nodes.push_back(node_info(position, node_value, idx));
    // std::cout << "node(" << nodes.size() - 1 << ") value is: " << nodes[nodes.size() - 1].value_ << std::endl;
    
    // adjacents.push_back(std::vector<size_type>());
    // edge_value.push_back(std::map<size_type, edge_value_type>());
    return node(idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (this == n.graph_ && nodes[n.node_index_].idx_ != -1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

    return Node(this, valid_nodes[i]);        // Invalid node
  }
  /* Remove a node from graph
  *@pre @a n is a node
  *@post if node is in the graph, remove the node from the graph, remove all edges that contains node from the graph
  *@return whether if the remove operation succeed
  *Complexity: O(d * log(num_edges())) = O(d*log(num_nodes()))
  */
  bool remove_node(const Node& n) {
    if (!has_node(n)) return 0;
    while(n.degree()){
      remove_edge(*(n.edge_begin()));
    }
    valid_nodes[n.index()] = valid_nodes[valid_nodes.size() - 1];
    nodes[valid_nodes[n.index()]].idx_ = n.index();

    valid_nodes.pop_back();   
    nodes[n.node_index_].idx_ = -1;
    return 1;
  }
  /* Remove a node from graph
  *@pre @a n_it is a node iterator
  *@post if node that @a n_it points to is in the graph, remove the node from the graph, remove all edges that contains node from the graph
  *@return the iterator
  *Complexity: O(d * log(num_edges())) = O(d*log(num_nodes()))
  */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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
    /** Return the node that the iterator points to. */
    Node operator*() const {

      return value_type(graph_, graph_ -> valid_nodes[idx_]);
    }

    /** Point to next node and return the node. */
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /** Return whether two nodes are same or not: graph and index should be the same. */
    bool operator==(const NodeIterator& nit) const {
      if (graph_ && nit.graph_) {
        return graph_ == nit.graph_ && idx_ == nit.idx_;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type idx_;
    NodeIterator (const Graph* graph, size_type idx) : graph_(const_cast<Graph*>(graph)), idx_(idx){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return the begin iterator of nodes*/
  node_iterator node_begin() const {

    return NodeIterator(this, 0);
  }
  /** Return the end iterator of nodes*/
  node_iterator node_end() const {

    return NodeIterator(this, valid_nodes.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
    /** Return the edge that the incidentIterator points to. */
    Edge operator*() const {

      return edge_type(graph_, node_index_, it_ -> first);
    }

    /** Return the next incidentIterator. */
    IncidentIterator& operator++() {
      ++it_;
      return *this;
    }
    /** Return if two interators are point to the same thing or not*/
    bool operator==(const IncidentIterator& iit) const{
      if (graph_ && iit.graph_) {
        return graph_ == iit.graph_ && node_index_ == iit.node_index_ && it_ == iit.it_;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_index_;
    typename std::unordered_map<size_type, edge_value_type>::iterator it_;
    IncidentIterator(const Graph* graph, size_type node_index, typename std::unordered_map<size_type, edge_value_type>::iterator it){
       graph_ = const_cast<Graph*>(graph);
       node_index_ = node_index;
       it_ = it;
    }
  };



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
      graph_ = NULL;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      assert(graph_);
      return graph_ -> node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // return true when graph_, node1_, and node2_ are the same;
      if (graph_ && e.graph_) {
        return graph_ == e.graph_ && node1_ == e.node1_ && node2_ == e.node2_;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // same as node
      return std::tie(graph_, node1_, node2_) < std::tie(e.graph_, e.node1_, e.node2_);
    }

    /** Return the length of edge
     *
     */
    double length() const {

      return norm(node1().position() - node2().position());
    }

    edge_type otherSide() const {

      return edge_type(graph_, node2_, node1_);
    }

    edge_value_type& value() {

      return graph_ -> nodes[node1_].valid_adjacents_[node2_];
    }

    const edge_value_type value() const {

      return graph_ -> nodes[node1_].valid_adjacents_[node2_];
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type node1_;
    size_type node2_;


    Edge(const graph_type* gragh, size_type node1, size_type node2) {
      graph_ = const_cast<graph_type*>(gragh);
      node1_ = std::min(node1, node2);
      node2_ = std::max(node1, node2);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // if (edge_index_nodes.size() != edge_nodes_index.size())
      // std::cout << "WRONG!!!!!! " << edge_nodes_index.size() " VS " << edge_index_nodes.size() << endl;
    return edge_nodes_index.size();
  }

  // size_type num_edges_node() const {
  //   return edge_nodes_index.size();
  // }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  // Edge edge(size_type i) const __attribute__((deprecated)) {
  //   // HW0: YOUR CODE HERE
  //   return *std::next(edge_begin(), i);
  // }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (b.node_index_ < a.node_index_) return has_edge(b, a);
    if (num_edges() == 0 || !has_node(a) || !has_node(b)) return false;
    return edge_nodes_index.count(std::pair<size_type, size_type>(a.node_index_, b.node_index_));
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE

    assert(has_node(a) && has_node(b));
    if (b.node_index_ < a.node_index_) return add_edge(b, a, value);
    if (has_edge(a, b)) return edge_type(this, a.node_index_, b.node_index_);   
    nodes[a.node_index_].valid_adjacents_[b.node_index_] = value;
    nodes[b.node_index_].valid_adjacents_[a.node_index_] = value;
    edge_nodes_index.insert(std::pair<std::pair<size_type, size_type>, size_type> (std::pair<size_type, size_type>(a.node_index_, b.node_index_), edge_nodes_index.size()));
    edge_index_nodes.push_back(std::pair<size_type, size_type>(a.node_index_, b.node_index_));
    return edge_type(this, a.node_index_, b.node_index_);
  }

  Edge edge(size_type i) {
    return edge_type(this, edge_index_nodes[i].first, edge_index_nodes[i].second);
  }

  /* Remove an edge from graph
  *@pre @a n1 and @a n2 are nodes in the graph
  *@post if @a n1 and @a n2 has edge in the graph, then remove the edge
  *@return whether if the remove operation succeed
  *Complexity: O(log(num_edges()))
  */

  bool remove_edge(const Node& n1, const Node& n2) {
    
    if (n2.node_index_ < n1.node_index_) return remove_edge(n2, n1); 
    if (!has_edge(n1, n2)) return 0;
    // std::cout << "before edges number:" << num_edges() << std::endl;
    nodes[n1.node_index_].valid_adjacents_.erase(n2.node_index_);
    nodes[n2.node_index_].valid_adjacents_.erase(n1.node_index_);
    size_type idx = edge_nodes_index[std::pair<size_type, size_type>(n1.node_index_, n2.node_index_)];
    if (idx != edge_index_nodes.size() - 1) {
      edge_index_nodes[idx] = edge_index_nodes[edge_index_nodes.size() - 1];
      edge_nodes_index[edge_index_nodes[idx]] = idx;
    }
    edge_index_nodes.pop_back();
    edge_nodes_index.erase(std::pair<size_type, size_type>(n1.node_index_, n2.node_index_));
    // std::cout << "end edges number:" << num_edges() << std::endl;
    return 1;

  }

  /* Remove an edge from graph
  *@pre @a e is an edge
  *@post if edge is in the graph, remove the edge from the graph
  *@return whether if the remove operation succeed
  *Complexity: O(log(num_edges()))
  */
  bool remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /* Remove an edge from graph
  *@pre @a e_it is an edge iterator
  *@post if edge that @a e_it points to is in the graph, remove the edge from the graph
  *@return the iterator
  *Complexity: O(log(num_edges()))
  */

  edge_iterator remove_edge(edge_iterator e_it) {
    return remove_edge(*e_it);
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    for (auto it = nodes.begin(); it != nodes.end(); ++it) {
      it -> valid_adjacents_.clear();
    }
    nodes.clear();
    valid_nodes.clear();
    edge_nodes_index.clear();
    edge_index_nodes.clear();
    // adjacents.clear();
  }





  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {

      // graph_ = NULL;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the edge that the iterator points to*/
    Edge operator*() const {

      return value_type(graph_, graph_ -> edge_index_nodes[idx_].first, graph_ -> edge_index_nodes[idx_].second);
    }

    /** Increase iterator and return the new iterator*/
    EdgeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /** Return if two edge iterator point to the same thing or not*/
    bool operator==(const EdgeIterator& eit) const {
      if (graph_ && eit.graph_) {
        return graph_ == eit.graph_ && idx_ == eit.idx_;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type idx_; 

    EdgeIterator(const graph_type* graph, size_type idx) {
      graph_ = const_cast<graph_type*> (graph);
      idx_ = idx;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return the begin iterator of edges*/
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);

  }
  /** Return the end iterator of edges*/
  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_index_nodes.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct node_info {
    Point position_;
    node_value_type value_;
    int idx_;
    std::unordered_map<size_type, edge_value_type> valid_adjacents_;
    node_info(const Point& position, const node_value_type& value, int idx) {
      position_ = position;
      value_ = value;
      idx_ = idx;
    }
  };

  // struct edge_info {
  //   edge_value_type value_;
  //   int idx_;
  //   edge_info (const edge_value_type& value, int idx) : value_(value), idx_(idx) {}
  // };

  std::vector<node_info> nodes;
  std::vector<size_type> valid_nodes;
  std::vector<std::pair<size_type, size_type>> edge_index_nodes;
  std::map<std::pair<size_type, size_type>, size_type> edge_nodes_index;
  // std::vector<std::vector<size_type>> adjacents;

  // std::vector<std::map<size_type, edge_value_type>> edge_value;
  // size_type graph_edge_index_;

};

#endif // CME212_GRAPH_HPP
