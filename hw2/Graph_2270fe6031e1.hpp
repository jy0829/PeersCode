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

  /** Synonym for Node Value type. */  
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Synonym for Edge Value type. */  
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
  Graph() : nodes_(), i2u_() {}

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
    Node() {}

    /** Return this node's position. */
    Point& position() {
      return graph_-> nodes_[uid_].p_;
    }
    
    /** Return this node's position, in const state. */
    const Point& position() const {
      return graph_-> nodes_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_ -> nodes_[uid_].idx_;
    }

    /** Return this node's value, an integer. */
    node_value_type& value() {
      return graph_-> nodes_[uid_].v_;
    }

    /** Return this node's value, an integer, in const state. */
    const node_value_type& value() const {
      return graph_-> nodes_[uid_].v_;
    }

    /** Return the number of edges incident to this node. */
    size_type degree() const {
      return graph_ -> nodes_[uid_].adj_.size();
    }

    /** Return the first incident iterator of this node */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    /** Return the last incident iterator of this node */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ and uid_ == n.uid_;
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
      return (graph_ < n.graph_) or (uid_ < n.uid_ and graph_ == n.graph_);
    }

   private:
    
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Data attributes of the Node class.
    Graph* graph_;
    size_type uid_;
    
    /** Construct a valid node. */
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
   
    nodeinfo new_node = nodeinfo(position, val, {}, num_nodes());
    nodes_.push_back(new_node);
    i2u_.push_back(nodes_.size() - 1);

    return Node(this, nodes_.size() - 1);    
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      return (this == n.graph_ and n.graph_ -> nodes_[n.uid_].idx_ != -1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return *std::next(node_begin(), i);  
  }

  /** Remove a node from the graph.
   * @return 0 if has_node(@a n) == false, @return 1 if has_node(@a n) ==  true. 
   *
   * @post If old has_node(@a n) == true then new has_node(@a n) == false 
   *                                          and @a n.index() == -1
   *                                          and new num_nodes() = old num_nodes() - 1.
   * @post Invalidates node indices: old node(@a i) might not equal new node(@a i). 
   *
   * @post Removes all incident edges to this node:
   *               If has_edge(@a n, @a b) == true then has_edge(@a n, @a b) == false for all
   *               nodes @a b adjacent to @a n.
   *
   * Complexity: O(num_nodes())
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)){
      return 0;
    }
    
    // Shift index of all valid nodes with index greater than @ n index down by 1
    for (size_type i = n.index()+1; i < num_nodes(); ++i) {
      nodes_[i2u_[i]].idx_ -= 1;
    }

    // Delete all edges incident to this node
    for (int i = n.degree()-1; i > -1; --i) {
      Node n2 = Node(this, nodes_[n.uid_].adj_[i].uid_);
      remove_edge(n,n2);
    }
      
    // Remove this node from index vector
    i2u_.erase(i2u_.begin() + n.index());

    // Invalidate this node
    nodes_[n.uid_].idx_ = -1;
      
    return 1;
  }

 /** Remove a node from the graph.
   * @pre node_begin() <= @a n_it < node_end()
   * @pre (* @a n_it) == a valid node, @a n.
   * @return node_begin() if has_node(@a n) == true else @return @a n_it.
   *
   * @post If old has_node(@a n) == true then new has_node(@a n) == false 
   *                                          and @a n.index() == -1
   *                                          and new num_nodes() = old num_nodes() - 1.
   * @post Invalidates node indices: old node(@a i) might not equal new node(@a i).
   *
   * @post Removes all incident edges to the node n:
   *               If has_edge(@a n, @a b) == true then has_edge(@a n, @a b) == false for all
   *               nodes @a b adjacent to @a n.
   *
   * Complexity: O(num_nodes())
   */
  node_iterator remove_node(node_iterator n_it) {
    if (remove_node((*n_it)) == 0) {
      return n_it;
    }
    return node_begin();
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, uid1_); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, uid2_);
    }

    /** Return this edge's value. */
    edge_value_type& value() {
      for (size_type i = 0; i < graph_->nodes_[uid1_].adj_.size(); ++i) {
	if (graph_->nodes_[uid1_].adj_[i].uid_ == uid2_) {
	  return graph_->nodes_[uid1_].adj_[i].v_;
	}
      }
      // Quiet comiler warning, above will return before reaching this point
      // if the edge exists.
      return graph_->nodes_[uid1_].adj_[0].v_;
    }

    /** Return this edge's value, in const state. */
    const edge_value_type& value() const {
      for (size_type i = 0; i < graph_->nodes_[uid1_].adj_.size(); ++i) {
	if (graph_->nodes_[uid1_].adj_[i].uid_ == uid2_) {
	  return graph_->nodes_[uid1_].adj_[i].v_;
	}
      }
      // Quiet comiler warning, above will return before reaching this point
      // if the edge exists.
      return graph_->nodes_[uid1_].adj_[0].v_;
    }

    
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_ and uid1_ == e.uid1_ and uid2_ == e.uid2_) or
	     (graph_ == e.graph_ and uid2_ == e.uid1_ and uid1_ == e.uid2_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return  (graph_ < e.graph_) or
	(graph_ == e.graph_ and uid1_ < e.uid1_) or
	(graph_ == e.graph_ and uid1_ == e.uid1_ and uid2_ < e.uid2_);
    }

   private:
    
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Data attributes of the edge class.
    Graph* graph_;
    size_type uid1_;
    size_type uid2_;

    /** Construct a valid edge. */
    Edge(const Graph* graph, size_type uid1, size_type uid2)
      : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {}
    
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
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto ei = a.edge_begin(); ei != a.edge_end(); ++ei) {
      if ((*ei).uid2_ == b.uid_) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    if (!has_edge(a,b)) {
      // Add the node to the nodes adjacency matrix
      edgeinfo new_edge1 = edgeinfo(val, b.uid_);
      edgeinfo new_edge2 = edgeinfo(val, a.uid_);
      nodes_[a.uid_].adj_.push_back(new_edge1);
      nodes_[b.uid_].adj_.push_back(new_edge2);
    }
    return Edge(this, a.uid_, b.uid_);
  }

 /** Remove an edge from the graph.
   * @return 0 if has_edge(@a a, @a b) == false, @return 1 if has_edge(@a a, @a b) == true
   *
   * @post If old has_edge(@a a, @a b) == true then new has_edge(@a a, @a b) == false
   *                                                and new num_edges() = old num_edges() -1.
   *
   * Complexity: O(max(@a a.degree(), @a b.degree()))
   */ 
  size_type remove_edge(const Node& n1, const Node& n2) {
    if(!has_edge(n1,n2)){
      return 0;
    }
    
    size_type aux1 = -1;
    for (size_type i = 0; i < n1.degree(); ++i) {
      if (nodes_[n1.uid_].adj_[i].uid_ == n2.uid_) {
        aux1 = i;
	break;
      }
    }

    size_type aux2 = -1;
    for (size_type i = 0; i < n2.degree(); ++i) {
      if (nodes_[n2.uid_].adj_[i].uid_ == n1.uid_) {
	aux2 = i;
	break;
      }
    }
    
    nodes_[n1.uid_].adj_.erase(nodes_[n1.uid_].adj_.begin() + aux1);
    nodes_[n2.uid_].adj_.erase(nodes_[n2.uid_].adj_.begin() + aux2);
    
    return 1;
  }

 /** Remove an edge from the graph.
   * @return 0 if has_edge(@a e.node1(), @a e.node2()) == false, return 1 if has_edge(@a e.node1(), @a e.node2()) == true
   *
   * @post If old has_edge(@a e.node1(), @a e.node2()) == true then new has_node(@a e.node1(), @a e.node2()) == false
   *                                                                and new num_edges() = old num_edges() -1.
   * Complexity: O(max(e.node1().degree(), e.node2().degree()))
   */ 
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

 /** Remove an edge from the graph.
   * @pre edge_begin() <= e_it < edge_end().
   * @pre (*e_it) == e, a valid edge.
   * @return e_it if has_edge(@a e.node1(), @a e.node2()) == false, return @begin_edge() if has_edge(@a e.node1(), @a e.node2()) == true
   *
   * @post If old has_edge(@a e.node1(), @a e.node2()) == true then new has_node(@a e.node1(), @a e.node2()) == false
   *                                                                and new num_edges() = old num_edges() -1.
   * Complexity: O(max(e.node1().degree(),e.node2().degree()))
   */   
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = (*e_it);
    if (remove_edge(e.node1(), e.node2()) == 0) {
      return e_it;
    }
    return edge_begin();
  }

  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
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
    NodeIterator() {}

    /** Return the node this node iterator is on */
    Node operator*() const {
      return Node(graph_, graph_->i2u_[idx_]);
    }

    /** Increment this node iterator by one, return this node iterator */
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /** Return if node iterator x is equal to this node iterator */
    bool operator==(const NodeIterator& x) const {
      return graph_ == x.graph_ and idx_ == x.idx_;
    }

   private:

    // Allow Graph to access NodeIterator's private member data and functions.
    friend class Graph;

    /** Construct a valid NodeIterator **/
    NodeIterator(const Graph* graph, size_type idx) : graph_(const_cast<Graph*>(graph)), idx_(idx) {}

    // Data attributes of the NodeIterator class
    Graph* graph_;
    size_type idx_;
  };

  /** Return the first node iterator of this graph */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return the last node iterator of this graph */
  NodeIterator node_end() const {
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

    /** Return the edge the incident iterator is on */
    Edge operator*() const {
      return Edge(graph_, uid1_, graph_->nodes_[uid1_].adj_[aux_uid2_].uid_);
    }

    /** Increment this incident iterator by one, return this incident iterator */
    IncidentIterator& operator++() {
      ++aux_uid2_;
      return *this;
    }
    
    /** Return if incident iterator iit is equal to this incident iterator */
    bool operator==(const IncidentIterator& iit) const {
      return graph_ == iit.graph_ and uid1_ == iit.uid1_ and aux_uid2_ == iit.aux_uid2_;  
    }

   private:
    friend class Graph;

    /** Construct a valid IncidentIterator **/
    IncidentIterator(const Graph* graph, size_type uid1, size_type aux_uid2) :
      graph_(const_cast<Graph*>(graph)), uid1_(uid1), aux_uid2_(aux_uid2) {}

    // Data attributes of this incident iterator
    Graph* graph_;
    size_type uid1_;
    size_type aux_uid2_;
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

    /** Return the edge this edge iterator is on. */
    Edge operator*() const {
      return Edge(graph_, graph_->i2u_[idx_], graph_->nodes_[graph_->i2u_[idx_]].adj_[aux_uid2_].uid_);
    }

    /** Increment this edge iterator by one, return this edge iterator */
    EdgeIterator& operator++() {
      ++aux_uid2_;
      fix();
      return *this;
    }

    /** Return if this edge iterator is equal to edge iterator x. */
    bool operator==(const EdgeIterator& x) const {
      return graph_ == x.graph_ and idx_ == x.idx_ and aux_uid2_ == x.aux_uid2_;
    }

   private:
    friend class Graph;

    /** Find the next valid edge */
    void fix() {
      while (true) {
        if (idx_ == graph_->num_nodes()) {
	  break;
	}
        else if (aux_uid2_ == graph_->nodes_[graph_->i2u_[idx_]].adj_.size()) {
          aux_uid2_ = 0;
	  ++idx_;
	}
        else if(graph_->i2u_[idx_] > graph_->nodes_[graph_->i2u_[idx_]].adj_[aux_uid2_].uid_) {
	  ++aux_uid2_;
	}
	else {
          break;
	}
      }
    }
    
    /** Construct a valid EdgeIterator **/
    EdgeIterator(const Graph* graph, size_type idx, size_type aux_uid2)
      : graph_(const_cast<Graph*>(graph)), idx_(idx), aux_uid2_(aux_uid2) {fix();}

    // Data attributes of this edge iterator
    Graph* graph_;
    size_type idx_;
    size_type aux_uid2_;
  };

  /** Return the first edge iterator of the graph */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }

  /** Return the last edge iterator of this graph */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_nodes(), 0);
  }


 private:

  // Data attributes of the Graph class.

  struct edgeinfo {
    edge_value_type v_;
    size_type uid_;

    edgeinfo(edge_value_type v, size_type uid) : v_(v), uid_(uid) {}
  };
  
  struct nodeinfo {
    Point p_;
    node_value_type v_;
    std::vector<edgeinfo> adj_;    
    int idx_;

    nodeinfo(Point p, node_value_type v, std::vector<edgeinfo> adj, size_type idx)
      : p_(p), v_(v), adj_(adj), idx_(idx) {} 
    
  };

  std::vector<nodeinfo>  nodes_;
  std::vector<size_type> i2u_;
  
};

#endif // CME212_GRAPH_HPP
