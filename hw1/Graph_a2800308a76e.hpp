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

//using namespace std;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */

  // HW0: YOUR CODE HERE
  Graph() : nodes_(), values_(), edges_(), adj_() {}


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
      return ngraph_->nodes_[uid_];
      //return Point();
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
      return ngraph_->values_[uid_];
    }
    const node_value_type& value() const{
      return ngraph_->values_[uid_];
    }

    size_type degree() const {
      return ngraph_->adj_[uid_].size();
    }

    incident_iterator edge_begin() const{
      return incident_iterator(ngraph_, uid_, 0);
    }
    incident_iterator edge_end() const{
      return incident_iterator(ngraph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // if (ngraph_ == n.ngraph_ && uid_ == n.uid_)
      //   return true;

      // //(void) n;          // Quiet compiler warning
      // return false;
      return (ngraph_ == n.ngraph_ and uid_ == n.uid_);
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
      std::ptrdiff_t pd = ngraph_ - n.ngraph_;
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
    graph_type* ngraph_;
    size_type uid_;

    //private constructor
    Node(const Graph* graph, size_type uid)
        : ngraph_(const_cast<Graph*>(graph)), uid_(uid){
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
    nodes_.push_back(position);
    values_.push_back(value);
    adj_.push_back(std::vector<size_type>());
    return Node(this, num_nodes()-1);
  } 

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // if (n.ngraph_ == this)
    //   return true;
    // //(void) n;            // Quiet compiler warning
    // return false;
    if (size() == 0) return false;
    return (n.ngraph_ == this);
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
      return Node(egraph_, uid1_);
      // return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(egraph_, uid2_);
      //return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (egraph_ == e.egraph_ and uid1_ == e.uid1_ and uid2_ == e.uid2_) or
       (egraph_ == e.egraph_ and uid1_ == e.uid2_ and uid2_ == e.uid1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      std::ptrdiff_t pd = egraph_ - e.egraph_;
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
    graph_type* egraph_;
    size_type uid1_;
    size_type uid2_;

    //private constructor
    Edge(const Graph* graph, size_type uid1, size_type uid2)
        : egraph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
    }   
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
    //return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, edges_[i].first, edges_[i].second);
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
    size_type id = adj_[a.uid_].size() < adj_[b.uid_].size()? a.uid_ : b.uid_;
    for (size_type i = 0; i < adj_[id].size(); ++i) {
      if (adj_[id][i] == a.uid_ or adj_[id][i] == b.uid_) {
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
    // HW0: YOUR CODE HERE
    if(!has_edge(a,b)) {
      edges_.push_back(std::make_pair(a.uid_, b.uid_));
      adj_[a.uid_].push_back(b.uid_);
      adj_[b.uid_].push_back(a.uid_);
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
    edges_.clear();
    values_.clear();
    adj_.clear();
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
    
    value_type operator*() const{
      return graph_->node(idx_);
      // return Node(graph_,idx_);
    }
    node_iterator& operator++(){
      ++idx_;
      return *this;
    }
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

  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

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
    value_type operator*() const{
      return Edge(graph_, uid_, graph_->adj_[uid_][idx_]);
    }

    incident_iterator& operator++(){
      ++idx_;
      return *this;
    }

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

    value_type operator*() const{
      return graph_->edge(idx_);
    }
    edge_iterator& operator++(){
      ++idx_;
      return *this;
    }
    bool operator==(const edge_iterator& rhs) const{
      return graph_ == rhs.graph_ and idx_ == rhs.idx_;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    graph_type* graph_;
    size_type idx_;

    //private constructor
    EdgeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {
    }  
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  edge_iterator edge_end() const {
    return edge_iterator(this, num_edges());
  }


 private:

  std::vector<Point> nodes_;
  std::vector<node_value_type> values_;

  //here we keep track both the vector and adjacent matrix representation
  //so that both incident iterator and edge iterator are easy to do
  std::vector<std::pair<size_type, size_type> > edges_;
  std::vector<std::vector<size_type> > adj_;


  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
