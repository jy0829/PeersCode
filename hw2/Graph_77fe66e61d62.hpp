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
template <typename V, typename E>
class Graph {

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

  /** Type of value associated with node */
  using node_value_type = V;

  /** Type of value associated with edge */
  using edge_value_type = E;

  // Data structure to store node info
  struct node_info {
    Point p_;
    node_value_type nval_;
    std::vector<size_type> adj_;
    int idx_;
    node_info(Point p, node_value_type nval,
              std::vector<size_type> adj, int idx):
              p_(p), nval_(nval), adj_(adj), idx_(idx) {}
  };

 private:

  /** Vector containing all the nodes and info associated
    * with each node, indexed by node uid */
  std::vector<node_info> nodes_;

  // Vector containing the edge values in a parallel adjacency list
  std::vector<std::vector<edge_value_type>> edge_values;

  // Vector containing the node indices
  std::vector<size_type> i2u_;

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() :
    nodes_(),
    edge_values(),
    i2u_() {
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
      return graph_->nodes_[uid_].p_;
    }

    /** Return this node's modifiable position. */
    Point& position() {
      return graph_->nodes_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].idx_;
    }

    // Return modifiable index
    size_type index() {
      return graph_->nodes_[uid_].idx_;
    }

    // Return the node's uid, a number in the range [0,nodes_.size()]
    size_type uid() const {
      return uid_;
    }

    /** Return a reference to the value associated with the node */
    node_value_type& value() {
      return graph_->nodes_[uid_].nval_;
    }

    /** Return the value associated with the node */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].nval_;
    }

    /** Return the degree of a node, i.e. the number of edges
     *  incident to it */
    size_type degree() const {
      return graph_->nodes_[uid_].adj_.size();
    }

    /** Iterator over edges incident to a given node, returns
     *  iterator associated with the first edge incident to the node */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,uid_,graph_->nodes_[uid_].adj_.begin());
    }

    /** Iterator over edges incident to a given node, returns
     *  iterator associated with the last edge incident to the node */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,uid_,graph_->nodes_[uid_].adj_.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_==n.graph_ && uid_==n.uid_);
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
      if (graph_!=n.graph_) {
        return (graph_ < n.graph_);
      }
      else {
        return (uid_ < n.uid_);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions
    friend class Graph;
    // Pointer back to the graph container
    graph_type* graph_;
    // The node's index in the graph
    size_type uid_;
    // Private constructor of Node object
    Node(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
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
                const node_value_type& val = node_value_type()) {
    nodes_.push_back(node_info(position,val,std::vector<size_type>(),i2u_.size()));
    i2u_.push_back(nodes_.size()-1);
    edge_values.push_back(std::vector<edge_value_type>());
    return Node(this, nodes_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (nodes_[n.uid()].idx_ != -1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i2u_[i]);
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
      return Node(graph_,uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,uid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_!=e.graph_) {
        return false;
      }
      else {
        return ((e.node1().uid()==uid1_ && e.node2().uid()==uid2_)
               ||(e.node1().uid()==uid2_ && e.node2().uid()==uid1_));
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_!=e.graph_) {
        return (graph_ < e.graph_);
      }
      else {
        return (this->node1().uid()+this->node2().uid()
                < e.node1().uid()+e.node2().uid());
      }
    }

    // Return the value of the Edge
    edge_value_type& value() {
      size_type minimum = std::min(this->node1().uid(),this->node2().uid());
      size_type maximum = std::max(this->node1().uid(),this->node2().uid());
      size_type loc = 0;
      for (size_type i=0; i < graph_->nodes_[minimum].adj_.size(); i++) {
        if (graph_->nodes_[minimum].adj_[i] == maximum) {
          loc = i;
          break;
        }
      }
      return graph_->edge_values[minimum][loc];
    }

    // Return a constant value of this Edge
    const edge_value_type& value() const {
      return this->value();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer to the graph that the Edge belongs to
    graph_type* graph_;
    // UID of the first node of the Edge
    size_type uid1_;
    // UID of the second node of the Edge
    size_type uid2_;
    // Private constructor for Edge
    Edge(const graph_type* graph,size_type uid1,size_type uid2)
      : graph_(const_cast<graph_type*>(graph)),uid1_(uid1),uid2_(uid2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return std::distance(edge_begin(),edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(),i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (!(has_node(a)||has_node(b))) {
      return false;
    }
    else {
      size_type count = 0;
      for (size_type i = 0; i < nodes_[a.uid()].adj_.size(); i++) {
        if (b.uid() == nodes_[a.uid()].adj_[i]) {
          count += 1;
        }
      }
      return (count == 1);
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

  Edge add_edge(const Node& a, const Node& b,
                const edge_value_type& val = edge_value_type()) {
    size_type uid1 = a.uid();
    size_type uid2 = b.uid();
    if (!has_edge(a,b)) {
      nodes_[uid1].adj_.push_back(uid2);
      nodes_[uid2].adj_.push_back(uid1);
      edge_values[uid1].push_back(val);
      edge_values[uid2].push_back(val);
    }
    return Edge(this,uid1,uid2);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edge_values.clear();
    i2u_.clear();
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

    /** Dereference operation, returns the node the iterator
     *  is pointing to */
    Node operator*() const {
      return Node(graph_,graph_->i2u_[idx_]);
    }

    /** Increment operation, increments uid by 1, returns
     *  a reference to the pointer to the new node */
    node_iterator& operator++() {
      ++idx_;
      return *this;
    }

    /** Comparison operation, returns whether the iterator is
     *  currently pointing to the node we are interested in */
    bool operator==(const node_iterator& nit) const {
      return (graph_==nit.graph_ && idx_==nit.idx_);
    }

   private:
    friend class Graph;
    // Pointer to the graph the node iterator belongs to
    graph_type* graph_;
    // Node's index in i2u_
    size_type idx_;
    // Private constructor called by the NodeIterator class
    NodeIterator(const graph_type* graph, size_type idx)
      : graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }

  };

  /** Node iterator begin */
  node_iterator node_begin() const {
    return node_iterator(this,0);
  }

  /** Node iterator end */
  node_iterator node_end() const {
    return node_iterator(this,num_nodes());
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

    /** Dereference operator, returns the edge the pointer is
     *  pointing to */
    Edge operator*() const {
      return Edge(graph_,uid_,*it_);
    }

    /** Increment operator, increments the position of the
     *  iterator of adjacent nodes by one */
    IncidentIterator& operator++() {
      ++it_;
      return *this;
    }

    /** Comparison operator, compares whether two incident
     *  iterators are pointing to the same adjacent node */
    bool operator==(const IncidentIterator& iit) const {
      return (graph_==iit.graph_ && uid_==iit.uid_ && it_==iit.it_);
    }

   private:
    friend class Graph;
    // Pointer to the graph associated with the incident iterator
    graph_type* graph_;
    // ID of the node we are interested in iterating over its incident edges
    size_type uid_;
    // Index of the adjacent node in the adjacency list
    std::vector<size_type>::iterator it_;
    // Private constructor for IncidentIterator
    IncidentIterator(const graph_type* graph, size_type uid,
                     std::vector<size_type>::iterator it)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid), it_(it) {
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator>  {
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

    /** Define edge dereference operator */
    Edge operator*() const {
      return Edge(graph_,(*boss_).uid(),(*adj_edge_).node2().uid());
    }

    /** Increment operator, returns a reference to the EdgeIterator
     *  and increments the index the EdgeIterator refers to */
    EdgeIterator& operator++() {
      ++adj_edge_;
      fix();
      return *this;
    }

    /** Comparison operator */
    bool operator==(const EdgeIterator& eit) const {
      return (graph_==eit.graph_ && boss_==eit.boss_
              && adj_edge_==eit.adj_edge_);
    }

   private:
    friend class Graph;
    // Pointer to the graph the edge iterator belongs to
    graph_type* graph_;
    //Node iterator
    NodeIterator boss_ = (*graph_).node_begin();
    // Incident iterator of the adjacent nodes
    IncidentIterator adj_edge_ = (*boss_).edge_begin();

    // Fix function to decide which edges to iterate over
    void fix() {
      while (boss_ != (*graph_).node_end()) {
        while (adj_edge_ != (*boss_).edge_end()) {
          if ((*adj_edge_).node2().uid() > (*boss_).uid()) {
            ++adj_edge_;
          }
          else {
            break;
          }
        }
        if (adj_edge_ == (*boss_).edge_end()) {
          ++boss_;
          adj_edge_ = (*boss_).edge_begin();
        }
        else {
          break;
        }
      }
      if (boss_ == (*graph_).node_end()) {
        boss_ = (*graph_).node_begin();
        adj_edge_ = (*boss_).edge_begin();
      }
    }

    // Private constructor called by the EdgeIterator class
    EdgeIterator(const graph_type* graph, NodeIterator boss,
                 IncidentIterator adj_edge) :
                 graph_(const_cast<graph_type*>(graph)),
                 boss_(boss), adj_edge_(adj_edge) {
      fix();
    }
  };

  /** Define the start of the edge iterator */
  edge_iterator edge_begin() const {
    return edge_iterator(this,node_begin(),(*node_begin()).edge_begin());
  }

  /** Define the end of the edge iterator */
  edge_iterator edge_end() const {
    return edge_iterator(this,node_end(),(*node_end()).edge_end());
  }

  //

  /** The following 2 functions will remove a given node from the graph
    * @post the number of nodes in the graph is reduced by 1
    * @post the number of edges in the graph is reduced by the
    *       degree of the node
    * Complexity of operation is O(degree()^2 + num_nodes())
    * The node and all its incident edges are deleted and invalidated
    */
  size_type remove_node(const Node& n) {
    size_type uid1 = n.uid();
    size_type idx1 = n.index();
    // Remove all edges incident to the node
    for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit) {
      auto e = *eit;
      assert(e.node1() == n);
      size_type uid2 = e.node2().uid();
      for (size_type i = 0; i < nodes_[uid2].adj_.size(); i++) {
        if (uid1 == nodes_[uid2].adj_[i]) {
          std::iter_swap(nodes_[uid2].adj_.begin()+i,nodes_[uid2].adj_.end()-1);
          std::iter_swap(edge_values[uid2].begin()+i,edge_values[uid2].end()-1);
          nodes_[uid2].adj_.pop_back();
          edge_values[uid2].pop_back();
          break;
        }
      }
    }
    //Remove the node and assign new indices
    auto remove = find(i2u_.begin(),i2u_.end(),uid1);
    std::iter_swap(remove,i2u_.end()-1);
    i2u_.pop_back();
    nodes_[uid1].adj_.clear();
    edge_values[uid1].clear();
    nodes_[i2u_[idx1]].idx_ = idx1;
    nodes_[uid1].idx_ = -1;
    return 1;
  }

  // Remove the node through a node iterator
  // The node iterator will be incremented, thus invaldiating the old one
  node_iterator remove_node(node_iterator n_it) {
    auto temp = n_it;
    ++temp;
    remove_node(*n_it);
    return temp;
  }


  /** The folllowing three functions will remove and edge
    * @post num_edges() is reduced by 1
    * @post the degree of the two nodes of the edge is reduced by 1
    * Complexity of operation is O(degree())
    * The edge is deleted from the adjacency lists of both nodes
    */
  // Remove an edge given the two nodes of the edge
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!has_edge(n1,n2)) {
      return false;
    }
    size_type uid1 = n1.uid();
    size_type uid2 = n2.uid();

    // Remove node 1 from node 2's adjacency list
    for (size_type i = 0; i < nodes_[uid2].adj_.size(); i++) {
      if (uid1 == nodes_[uid2].adj_[i]) {
        std::iter_swap(nodes_[uid2].adj_.begin()+i, nodes_[uid2].adj_.end()-1);
        std::iter_swap(edge_values[uid2].begin()+i, edge_values[uid2].end()-1);
        nodes_[uid2].adj_.pop_back();
        edge_values[uid2].pop_back();
        break;
      }
    }

    // Remove node 2 from node 1's adjacency list
    for (size_type i = 0; i < nodes_[uid1].adj_.size(); i++) {
      if (uid2 == nodes_[uid1].adj_[i]) {
        std::iter_swap(nodes_[uid1].adj_.begin()+i, nodes_[uid1].adj_.end()-1);
        std::iter_swap(edge_values[uid1].begin()+i, edge_values[uid1].end()-1);
        nodes_[uid1].adj_.pop_back();
        edge_values[uid1].pop_back();
        break;
      }
    }

    return 1;
  }

  // Remove a given edge, calls the remove_edge function defined above
  size_type remove_edge(const Edge& e) {
    auto node1 = e.node1();
    auto node2 = e.node2();
    remove_edge(node1, node2);
    return 1;
  }

  /** Remove an edge given an edge iterator, calls the remove_edge function
    * defined above */
  edge_iterator remove_edge(edge_iterator e_it) {
    auto temp = e_it;
    ++temp;
    remove_edge(*e_it);
    return temp;
  }

};

#endif // CME212_GRAPH_HPP
