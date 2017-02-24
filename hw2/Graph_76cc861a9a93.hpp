#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iterator>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

//using namespace std;

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

  // using size_type = unsigned;


  // std::vector<size_type[2]> edges;

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
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */

  // HW0: YOUR CODE HERE
  Graph() : nodes_(), num_edges_(0), adj_() {}

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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].position_;
      //return Point();
    }

    /** Return this node's position. */
    Point& position() {
      return graph_->nodes_[uid_].position_;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
      //return size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    node_value_type& value(){
      return graph_->nodes_[uid_].value_;
    }
    const node_value_type& value() const{
      return graph_->nodes_[uid_].value_;
    }

    size_type degree() const {
      return graph_->adj_[uid_].size();
    }

    incident_iterator edge_begin() const{
      return incident_iterator(graph_, uid_, 0);
    }
    incident_iterator edge_end() const{
      return incident_iterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // if (graph_ == n.graph_ && uid_ == n.uid_)
      //   return true;

      // //(void) n;          // Quiet compiler warning
      // return false;
      return (graph_ == n.graph_ and uid_ == n.uid_);
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
      std::ptrdiff_t pd = graph_ - n.graph_;
      return (pd < 0) or (pd ==0 and uid_ < n.uid_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // friend class Edge;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type uid_;

    //private constructor
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid){
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();
    //return 0;
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

  // Node add_node(const Point& position) {
  //   // HW0: YOUR CODE HERE
  //   nodes_.push_back(position);
  //   return Node(this, num_nodes()-1);
  // }

  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    nodes_.push_back(nodeinfo(position, value));
    adj_.push_back(std::vector<edgeinfo>());
    return Node(this, num_nodes()-1);
  } 

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // if (n.graph_ == this)
    //   return true;
    // //(void) n;            // Quiet compiler warning
    // return false;
    if (size() == 0) return false;
    return (n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i);
    //(void) i;             // Quiet compiler warning
    //return Node();        // Invalid node
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
      assert(uid1_>=0 and uid1_<graph_->num_nodes());
      return Node(graph_, uid1_);
      // return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      assert(uid2_>=0 and uid2_<graph_->num_nodes());
      return Node(graph_, uid2_);
      //return Node();      // Invalid Node
    }

    /** Return the length of the edge
     */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    edge_value_type& value_helper() {
      size_type ida = graph_->adj_[uid1_].size() < graph_->adj_[uid2_].size()? uid1_ : uid2_;
      size_type idb = std::max(uid1_, uid2_) - ida + std::min(uid1_, uid2_);
      size_type i;
      for (i = 0; i < graph_->adj_[ida].size(); ++i) {
        if (graph_->adj_[ida][i].pid_ == idb) {
          break;
        }      
      }
      assert(i < graph_->adj_[ida].size());
      return graph_->adj_[ida][i].value_;
    }

    edge_value_type& value() {
      return value_helper();
    }

    const edge_value_type& value() const {
      return value_helper();    
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_ and uid1_ == e.uid1_ and uid2_ == e.uid2_) or
       (graph_ == e.graph_ and uid1_ == e.uid2_ and uid2_ == e.uid1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      std::ptrdiff_t pd = graph_ - e.graph_;
      return (pd < 0) or (pd == 0 and uid1_ < e.uid1_) or 
        (pd == 0 and uid1_ == e.uid1_ and uid2_ < e.uid2_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type uid1_;
    size_type uid2_;

    //private constructor
    Edge(const Graph* graph, size_type uid1, size_type uid2)
        : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
    }   
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
    //return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * here I don't use std::next() due to compiler issues
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //return Edge(this, edges_[i].first, edges_[i].second);
    size_type id = 0;
    auto it = edge_begin();    
    while (id++ < i) {
      ++it;
    }
    return *it;
    //(void) i;             // Quiet compiler warning
    //return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  //here we search based on the adjacent matrix representation, which is much 
  //faster than vector representation
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // for (size_type i = 0; i < num_edges(); ++i) {
    //   if ((a.uid_ == edges_[i].first and b.uid_ == edges_[i].second) or
    //       (a.uid_ == edges_[i].second and b.uid_ == edges_[i].first)) {
    //     return true;
    //   }
    // }
    size_type ida = adj_[a.uid_].size() < adj_[b.uid_].size()? a.uid_ : b.uid_;
    size_type idb = std::max(a.uid_, b.uid_) - ida + std::min(a.uid_, b.uid_);
    for (size_type i = 0; i < adj_[ida].size(); ++i) {
      if (adj_[ida][i].pid_ == idb) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if(!has_edge(a,b)) {
      adj_[a.uid_].push_back(edgeinfo(b.uid_, value));
      adj_[b.uid_].push_back(edgeinfo(a.uid_, value));
      ++num_edges_;
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    adj_.clear();
    num_edges_ = 0;
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

    /** Return the node that the iterator points to*/
    value_type operator*() const{
      return graph_->node(idx_);
      // return Node(graph_,idx_);
    }

    /** Increase iterator to next object and return the new iterator*/
    node_iterator& operator++(){
      ++idx_;
      return *this;
    }

    /** Return if two edge iterator point to the same object or not*/
    bool operator==(const node_iterator& rhs) const{
      return graph_ == rhs.graph_ and idx_ == rhs.idx_;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type idx_;

    //private constructor
    NodeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx){
    }    
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return the begin iterator of edges*/
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  /** Return the end iterator of edges*/
  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
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
    // IncidentIterator operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return the edge that the iterator points to*/
    value_type operator*() const{
      return Edge(graph_, uid_, graph_->adj_[uid_][idx_].pid_);
    }

    /** Increase iterator to next object and return the new iterator*/
    incident_iterator& operator++(){
      ++idx_;
      return *this;
    }

    /** Return if two edge iterator point to the same object or not*/
    bool operator==(const incident_iterator& rhs) const{
      return graph_ == rhs.graph_ and uid_ == rhs.uid_ and idx_ == rhs.idx_;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type uid_;
    size_type idx_;
    IncidentIterator(const Graph* graph, size_type uid, size_type idx)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), idx_(idx){
    } 

    // HW1 #3: YOUR CODE HERE
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

    /** Return the edge that the iterator points to*/
    value_type operator*() const{
      return Edge(graph_, uid_, graph_->adj_[uid_][idx_].pid_);
    }

   /** Operator to get to the next egde
    * by looking up in the adjacent matrix
    *
    * To make sure we do not iterate over every edge twice, we return edges
    * where the index of the first node is lower or equal than the index of the second node.
    *
    */
    edge_iterator& operator++(){
      ++idx_;
      while (uid_ < graph_->num_nodes()) {
        while (idx_ < graph_->adj_[uid_].size()) {
          if (uid_ < graph_->adj_[uid_][idx_].pid_)
            return *this;
          ++idx_;
        }
        ++uid_;
        idx_ = 0;
      }
      return *this;
    }

    /** Return if two edge iterator point to the same object or not*/
    bool operator==(const edge_iterator& rhs) const{
      return graph_ == rhs.graph_ and uid_ == rhs.uid_ and idx_ == rhs.idx_;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type uid_;
    size_type idx_;
    EdgeIterator(const Graph* graph, size_type uid, size_type idx)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), idx_(idx){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return the begin iterator of edges*/
  edge_iterator edge_begin() const {
    size_type uid = 0;
    while (uid < num_nodes() and adj_[uid].size() == 0) {
      ++uid;
    }
    return edge_iterator(this, uid, 0);
  }

  /** Return the end iterator of edges*/
  edge_iterator edge_end() const {
    return edge_iterator(this, num_nodes(),0);
  }

  /** remove a node to the graph and all the edges it belongs to
   * if the node is not existing, directly @return 0, the graph remain unchanged
   * if the node exists
   * @return 1
   * @post new num_nodes() == old num_nodes() - 1.
   * @post new num_edges() == old num_edges() - n.degree().
   * 
   * invalidate edge indexes, node indexes, node iterators,
   * edge iterators, and incident iterators
   *
   * Complexity: No more than O(numedges())
   */
  size_type remove_node(const Node& n){
    if (!has_node(n)) return 0;
    num_edges_ -= n.degree();
    size_type id = n.index();
    auto pred = [id](const edgeinfo & ef) {
      return ef.pid_ == id;
    };
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      Node n1 = e.node1();
      assert(n1 == n);
      Node n2 = e.node2();  // Adjacent node!
      adj_[n2.uid_].erase(std::remove_if(adj_[n2.uid_].begin(), 
        adj_[n2.uid_].end(), pred), adj_[n2.uid_].end());
    }

    adj_.erase(adj_.begin()+n.index());
    nodes_.erase(nodes_.begin()+n.index());

    for (size_type i = 0; i < num_nodes(); ++i) {
      for (size_type j = 0; j < adj_[i].size(); ++j) {
        if (adj_[i][j].pid_> id) adj_[i][j].pid_--;
      }
    }
    return 1;
  }

  /** similar to remove_node(const Node& n)
   * except @param[in] is an node iterator
   * and @return same as @param[in]
   */  
  node_iterator remove_node(node_iterator n_it){
    remove_node(*n_it);
    return n_it;
  }

  /** remove an edge to the graph
   * if the edge is not existing, directly @return 0, the graph remain unchanged
   * @pre @a a and @a b are distinct valid nodes of this graph
   * if the edge exists
   * @return 1
   * @post new num_edges() == old num_edges() - 1.
   *
   * invalidate edge indexes, node indexes, node iterators,
   * edge iterators, and incident iterators
   *
   * Complexity: No more than O(max(node.degreee()))
   */
  size_type remove_edge(const Node& a, const Node& b){
    if (!has_edge(a,b)) return 0;
    size_type aid = a.index();
    size_type bid = b.index();
    auto pred_a = [bid](const edgeinfo & ef) {
      return ef.pid_ == bid;
    };
    auto pred_b = [aid](const edgeinfo & ef) {
      return ef.pid_ == aid;
    };    
    adj_[a.uid_].erase(std::remove_if(adj_[a.uid_].begin(), 
      adj_[a.uid_].end(), pred_a), adj_[a.uid_].end());
    adj_[b.uid_].erase(std::remove_if(adj_[b.uid_].begin(), 
      adj_[b.uid_].end(), pred_b), adj_[b.uid_].end());
    --num_edges_;    
    return 1;
  }

  /** similar to remove_edge(const Node& a, const Node& b)
   * except @param[in] is an edge
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(),e.node2());
  }

  /** similar to remove_edge(const Node& a, const Node& b)
   * except @param[in] is an edge iterator
   * and @return same as @param[in]
   */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return e_it;
  }
 
 private:

  struct nodeinfo {
    Point position_;
    node_value_type value_;
    nodeinfo(const Point& position, const node_value_type& value) {
      position_ = position;
      value_ = value;
    }
  };

  struct edgeinfo {
    size_type pid_;
    edge_value_type value_;
    edgeinfo(const size_type& pid, const edge_value_type& value) {
      pid_ = pid;
      value_ = value;
    }
  };  

  std::vector<nodeinfo> nodes_;
  size_type num_edges_;
  std::vector<std::vector<edgeinfo> > adj_;

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
