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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;
  struct internal_edge;

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

  /** Construct an empty graph. */
  Graph()
  : internal_nodes(), internal_edges(), next_node_uid(0), next_edge_uid(0) {
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
  class Node {
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
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch().uid;
    }





    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /**Returns the node value
      @return The value stored by the node
      */
    node_value_type& value(){
      return fetch().value;
    }

    /**Returns the node value
      @return The values stored by the node
    */
    const node_value_type& value() const {
      return fetch().value;
    }

    /**Sets value of node
      @pre value is a valid set_value
    */

    void set_value(node_value_type v){
      fetch().value = v;
    }
    /**Returns degree of node
    @return Degree of node
    */
    size_type degree() const{
      return fetch().adjacent_edges.size();
    }

    /**Returns initialized edge iterator
    @return Edge iterator initialized at edge 0
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(0, this);
    }

    /**Returns invalid edge iterator pointing to the position after the last edge
    @return Edge iterator initialized at edge degree()
    */
    incident_iterator edge_end() const{
      return IncidentIterator(degree(), this);
    }



    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((n.index()==index()) && (n.graph_pointer()==graph_pointer())){
        return true;
      } else {
        return false;
      }
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
      if (graph_<n.graph_pointer())
      {
        return true;
      } else{
        if (graph_ == n.graph_pointer())
        {
          if (index()<n.index()){
            return true;
          } else {
            return false;
          }
        } else{
          return false;
        }
      }

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // A pointer to the graph
    Graph<V>* graph_;
    // A unique identification number
    size_type uid_;

    /** Private Constructor */
    Node(const Graph<V>* graph, size_type uid)
        : graph_(const_cast<Graph<V>*>(graph)), uid_(uid) {
    }

    // Fetches the internal node corresponding to the node
    internal_node& fetch() const {
      return graph_->internal_nodes[uid_];

    }

    // Return the pointer to the graph
    const Graph<V>* graph_pointer() const {
      // HW0: YOUR CODE HERE
      return graph_;
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return next_node_uid;
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
  Node add_node(const Point& position,
    const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    internal_node new_internal_node =
      internal_node(position, value, next_node_uid);
    internal_nodes.push_back(new_internal_node);
    ++next_node_uid;
    return Node(this, next_node_uid-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_pointer() == this){
      if (n.index() < next_node_uid){
        return true;
      }
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch().uid;
    }

    void change_orientation(){
      orientation_ = !orientation_;
    }


    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if (orientation_){
        return  graph_->node(fetch().uid1);
      } else {
        return graph_->node(fetch().uid2);
      }
            // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (orientation_){
        return graph_->node(fetch().uid2);
      } else {
        return graph_->node(fetch().uid1);
      }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==( Edge& e) const {
      if (e.index()==index() && e.graph_pointer()==graph_){
        return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (const_cast<Graph<V>*>(graph_)<e.graph_pointer())
      {
        return true;
      } else{
        if (const_cast<Graph<V>*>(graph_) == e.graph_pointer())
        {
          if (index()<e.index()){
            return true;
          } else {
            return false;
          }
        } else{
          return false;
        }
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph<V>* graph_;
    size_type uid_;
    bool orientation_;

    /** Private Constructor */
    Edge(const Graph<V>* graph, size_type uid, bool orientation)
        : graph_(const_cast<Graph<V>*>(graph)), uid_(uid), orientation_(orientation) {
    }

    // Fetches the internal edge corresponding to the edge
    internal_edge& fetch() const {
      return (graph_->internal_edges)[uid_];
    }

    // Return the pointer to the graph
    const Graph<V>* graph_pointer() const{
      // HW0: YOUR CODE HERE
      return graph_;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return next_edge_uid;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i, true);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (unsigned int i = 0; i < next_edge_uid; ++i){
      const internal_edge* edge = &internal_edges[i];
      size_type uid1 = a.index();
      size_type uid2 = b.index();
      if ((edge->uid1 == uid1 && edge->uid2 == uid2) ||
        (edge->uid1 == uid2 && edge->uid2 == uid1)){
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

    internal_edges.push_back(internal_edge(next_edge_uid, a.fetch().uid, b.fetch().uid));
    a.fetch().adjacent_edges.push_back(next_edge_uid);
    b.fetch().adjacent_edges.push_back(next_edge_uid);



    ++next_edge_uid;
    return Edge(this, next_edge_uid-1, true);        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    internal_nodes.clear();
    internal_edges.clear();
    next_node_uid = 0;
    next_edge_uid = 0;
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
    NodeIterator(int ind_, const Graph<V>* graph_){
      ind = ind_;
      graph = graph_;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    // Returns index
    int index() const{
      return ind;
    }

    const Graph<V>* graph_pointer() const{
      return graph;
    }
    /**Dereferenced node object
    @return Dereferenced node of type Node
    */
    Node operator*() const{
      return graph->node(ind);
    }
    /**Returns next iterator
    @return Returns iterator that points to the node with the next uid
    */
    node_iterator& operator++(){
      ind++;
      return *this;
    }
    /** Determines whether two iterators are equal
    @return boolean equal to true iff both iterators point to the same node_iterator
    and belong to the same graph
    */
    bool operator==(const node_iterator& iter) const{
      if (ind == iter.index() && const_cast<Graph<V>*>(graph) == iter.graph_pointer()){
        return true;
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    int ind;
    const Graph<V>* graph;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /* Initialized node iterator
  @return Returns node iterator pointing to the node of uid 0
  */
  node_iterator node_begin() const{
    return NodeIterator(0, this);
  }
  /* Invalid initialized node iterator
  @return Returns node iterator pointing to the uid of a new node
  */
  node_iterator node_end() const{
    return NodeIterator(num_nodes(), this);
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



    IncidentIterator(int ind_, const Node* node_){
      ind = ind_;
      node = node_;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    // Returns index
    int index() const{
      return ind;
    }

    const Node* node_pointer() const{
      return node;
    }
    /* Dereferences edge iterator
    @return A valid edge with uid ind_
    */
    Edge operator*() const{
      int index = (node->fetch()).adjacent_edges[ind];

      Edge dereferenced_edge = node->graph_pointer()->edge(index);


      if (dereferenced_edge.node1().index()!=(node->index())){
        dereferenced_edge.change_orientation();
      }
      assert(dereferenced_edge.node1().index()==(node->index()) || dereferenced_edge.node2().index()==(node->index()));




      return dereferenced_edge;
    }

    /* Advances the iterator
    @return An iterator pointing to the next node (in uid order)
    */
    incident_iterator& operator++(){
      ind++;
      return *this;
    }
    /* Determines if two iterators are equal
    @return true iff both iterators belong to the same graph and point to the
    same edge
    */
    bool operator==(const incident_iterator& iter) const{
      if (ind == iter.index() && node == iter.node_pointer()){
        return true;
      } else {
        return false;
      }
    }

    private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    int ind;
    const Node* node;

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
    EdgeIterator(int ind_, const Graph<V>* graph_){
      ind = ind_;
      graph = graph_;
    }
    // Returns index
    int index() const{
      return ind;
    }

    const Graph<V>* graph_pointer() const{
      return graph;
    }
    /* Dereferences edge iterator
    @return Valid edge to which iterator is pointing
    */
    Edge operator*() const{
      return graph->edge(ind);
    }
    /* Advances iterator
    @return A valid edge with the next uid
    */
    edge_iterator& operator++(){
      ind++;
      return *this;
    }
    /* Determines if two iterators are equal
    @return true iff both iterators belong to the same graph and point to the
    same edge
    */
    bool operator==(const edge_iterator& iter) const{
      if (ind == iter.index() && const_cast<Graph<V>*>(graph) == iter.graph_pointer()){
        return true;
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    int ind;
    const Graph<V>* graph;

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  /*Initializes iterator
  @return iterator pointing to first edge of the graph
  */
  edge_iterator edge_begin() const{
    return edge_iterator(0, this);
  }
  /*Initializes iterator at end position
  @return iterator pointing to the position of a new edge of the graph
  */
  edge_iterator edge_end() const{
    return edge_iterator(next_edge_uid, this);
  }



 private:
   struct internal_node {
     internal_node( Point point_, node_value_type value_,size_type uid_){
       point = point_;
       value = value_;
       uid = uid_;
     }
     Point point;         // The position of the node
     node_value_type value;
     size_type uid;      // The unique identifcation for the node
     std::vector<size_type> adjacent_edges; //The edges adjacent to the node
   };

   struct internal_edge {
     internal_edge(size_type uid_, size_type uid1_, size_type uid2_){
       uid = uid_;
       uid1 = uid1_;
       uid2 = uid2_;
     }
     size_type uid;     //The unique identification for the edge
     size_type uid1;      // Pointer to node1
     size_type uid2;  // Pointer to node2
   };


  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Define a vector of nodes and edges, as well as the number of nodes, edges
  // and the next unique identifier
public:
  std::vector<internal_node> internal_nodes;
  std::vector<internal_edge> internal_edges;
  size_type next_node_uid;
  size_type next_edge_uid;


};

#endif // CME212_GRAPH_HPP
