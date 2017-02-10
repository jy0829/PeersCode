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
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

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

  /** Construct an empty graph.  */
  //This graph has zero edges, zero nodes. Initialize the right constants.
  Graph()
	: nodes_(), edges_(), adjacency() {}
    // HW0: YOUR CODE HERE

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
      return graph_ -> nodes_[uid_].first;
      // HW0: YOUR CODE HERE
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    //Return value of this node.
    node_value_type& value(){
      return graph_ -> nodes_[uid_].second;
    }
    //Same as above, but cast as a constant.
    const node_value_type& value() const{
      return graph_ -> nodes_[uid_].second;
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (uid_ == n.index() and n.graph_ == graph_);
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
      if (graph_ != n.graph_)
        return graph_ < n.graph_;
      else if (norm_2(position()) != norm_2(n.position()))
        return  norm_2(position()) < norm_2(n.position());
      else
        return &position() < &n.position();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
     friend class Graph;
     Graph* graph_; //The graph to which the node belongs.
     size_type uid_; //This node's uid.
     Node(const Graph* graph, size_type uid)
       : graph_(const_cast<Graph*>(graph)), uid_(uid){}
    //Private constructor for a node.
   public:

   //Return the degree of this node: number of edges attached.
   size_type degree() const{
     return graph_ -> adjacency[uid_].size();
   }
   //The beginning of this incident iterator.
   IncidentIterator edge_begin() const{
     return IncidentIterator(graph_, uid_, 0);
   }
   //The edge of this iterator: we would not like to be past this.
   IncidentIterator edge_end() const{
     return IncidentIterator(graph_, uid_, degree());
   }

}; //End Node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();
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
    //We first add the node to our nodes vector.
    nodes_.push_back(std::make_pair(position, val));
    //Note, we also need to update the adjacency matrix.
    //I.e. create a new vector inside of the matrix.
    std::vector<size_type> N;
    adjacency.push_back(N);
    return Node(this, num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (num_nodes() == 0){ //Check if graph is empty.
      return false;
    }
    else{ //If not empty.
      return n.graph_ == this;
    }
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

    //Our architecture doesn't let us access the second node
    //of the edge easily. We need to write a function to extract
    //it from the adjacency matrix.

    //This allows us to obtain the uid of the second node of this edge.
    size_type fetch_uid2() const{
      return graph_ -> adjacency[node_uid_][n_node_uid_];
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, fetch_uid2());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_ and e.node_uid_ == node_uid_ and e.fetch_uid2() == fetch_uid2())
        or (graph_ == e.graph_ and e.node_uid_ == fetch_uid2() and e.fetch_uid2() == node_uid_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //Check if the two graphs are the same.
      if(graph_ != e.graph_)
        return graph_ < e.graph_;
      else if (e.node_uid_ == node_uid_)
        return fetch_uid2() < e.fetch_uid2();
      else
        return e.node_uid < node_uid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_; //Graph containing this edge.
    size_type node_uid_; //uid of first node.
    size_type n_node_uid_; //place in vector of second node.
    Edge(const Graph* graph, size_type u1, size_type u2) :
      graph_(const_cast<Graph*>(graph)), node_uid_(u1),
      n_node_uid_(u2) {} //Constructor.


  }; //End Edge Class

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, edges_[i].first, edges_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (unsigned int i = 0; i < adjacency[a.node_uid_].size(); i++){
      if (b.node_uid_ == adjacency[a.node_uid_][i])
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
    for (unsigned int i = 0; i < adjacency[a.uid_].size(); i++){
      if(b.uid_ == adjacency[a.uid_][i]){
      //In this case, the edge already exists.
        return Edge(this, a.uid_, i);
      }
    }
    //Then you add both the uid of the first node and the index
    //of the second uid (not the uid itself!!).
    adjacency[a.uid_].push_back(b.uid_);
    adjacency[b.uid_].push_back(a.uid_);
    //Also add the relevent info in the edges vector.
    edges_.push_back(std::make_pair(a.uid_, adjacency[a.uid_].size() - 1));
    return Edge(this, a.uid_, adjacency[a.uid_].size() - 1);
    //The last argument is the index in that vector.
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    adjacency.clear();
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
    Node operator*() const{
      //Return the current node's uid.
      return Node(graph_,uid_);
    }
    NodeIterator& operator++(){
      //Change current node to the next one (in terms of uid).
      uid_++;
      return *this;
    }

    bool operator==(const NodeIterator& x) const{
      //Check for equality between this current node, and x's.
      return graph_ == x.graph_ and uid_ == x.uid_;
    }
   private:
    friend class Graph;
    Graph* graph_; //The uid of the node that the iterator is at.
    size_type uid_;
    //Constructor.
    NodeIterator(const Graph* graph, size_type uid) :
    graph_(const_cast<Graph*>(graph)), uid_(uid) {}
    // HW1 #2: YOUR CODE HERE
  }; //End Node Iterator Class

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  NodeIterator node_begin() const{
    //Return the first node of this graph.
    return NodeIterator(this, 0);
  }
  NodeIterator node_end() const{
    //This element does not actually exist: one past the end.
    return NodeIterator(this, num_nodes());
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
    Edge operator*() const{
    //Return the current edge.
      return Edge(graph_, node_uid_, n_node_uid_);
    }
    IncidentIterator& operator++(){
      //We simply increment the place in the vector of the next node.
      n_node_uid_++;
      return *this;
    }

    bool operator==(const IncidentIterator& x) const{
      return node_uid_ == x.node_uid_ and n_node_uid_ == x.n_node_uid_;
    }

   private:
    friend class Graph;
    size_type node_uid_; //uid of the current node
    size_type n_node_uid_;
    //This above is not actually an uid.
    Graph* graph_; //Pointer to the graph.

    IncidentIterator(const Graph* graph, size_type u1, size_type u2):
    graph_(const_cast<Graph*>(graph)), node_uid_(u1), n_node_uid_(u2) {}

    // HW1 #3: YOUR CODE HERE
  }; //End Incident Iterator class

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
    Edge operator*() const{
      //Access the edge using index and the edge() method.
      return graph_ -> edge(idx);
    }
    EdgeIterator& operator++(){
      //In this architecture, we only have to increment index.
      ++idx;
      return *this;
    }
    bool operator==(const EdgeIterator& x) const{
      //Two edge iterators are equal if they're on the same
      //graph and on the same index.
      return graph_ == x.graph_ and idx == x.idx;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type idx; //Index of the edge in the edges vectors - uid.
    Graph* graph_; //Pointer to the graph.
    EdgeIterator(const Graph* graph, size_type id):
      graph_(const_cast<Graph*>(graph)), idx(id) {}
    //Private constructor for edge class.
  }; //End Edge Iterator class

  //Beginning of your edge iterator.
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  }
  //End of the iterator - don't want to reach this stage: hence
  // the num_edges().
  edge_iterator edge_end() const{
   return EdgeIterator(this, this -> num_edges());
  }

  private:

  std::vector<std::pair<Point, V>> nodes_; //vector of nodes.
  std::vector<std::pair<size_type, size_type>> edges_; //vector of edges.
  std::vector<std::vector<size_type>> adjacency;
  //Adjacency matrix.

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

}; //End Graph class

#endif // CME212_GRAPH_HPP
