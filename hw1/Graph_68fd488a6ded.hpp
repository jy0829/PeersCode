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
  // Internal types are at the bottom

 public:

//
// PUBLIC TYPE DEFINITIONS
//

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  
  /** Synonym for V **/
  using node_value_type = V;

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
      : nodes_(), incident_matrix_() {
      	num_edges_ = 0;
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
  class Node : private totally_ordered<Node>{
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
      return graph_->nodes_[uid_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

   /** Return this node's value, a value of type node_value_type.
    * @return Value of the node.
    */
    node_value_type& value() {
    	return graph_->nodes_[uid_].second;
    }

   /** Return this node's value, a value of type node_value_type.
    * For when called on a constant type.
    * @return Value of the node.
    */
    const node_value_type& value() const {
    	return graph_->nodes_[uid_].second;
    }

   /** Return the degree of this node.
    * @return Number of edges incident to this node.
    */
    size_type degree() const {
    	return graph_->incident_matrix_[uid_].size();
    }

   /** Return incident_iterator pointing to the first edge incident to this node or the 
    * end iterator if there are no incident edges
    * @return Incident_iterator pointing to the first incident edge, or the end iterator if
    * 	there are no edges incident to this node.
    */
    incident_iterator edge_begin() const {
    	return incident_iterator(graph_, uid_, 0);
    }

   /** Return incident_iterator one past the last incident edge, if there are any,
    * of this node.
    * @return Incident_iterator one past the last incident edge of this node.
    * @post result must not be able to be dereferenced
    */
    incident_iterator edge_end() const {
    	return incident_iterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.uid_ == uid_ && n.graph_ == graph_) {
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
      if (graph_ == n.graph_ && uid_ < n.uid_) {
        return true;
      } else if (graph_ < n.graph_) {
      	return true;
      } else {
      	return false;
      }
    }

   private:
    // Pointer back to the node container
    Graph* graph_;
    // This node's unique identification number
    size_type uid_;
    // This node's value
    node_value_type value_;
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        	assert(graph != nullptr);
        	assert(uid < graph->num_nodes());
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value. Default is node_value_type holder.
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    std::pair<Point, node_value_type> new_pair (position, value);
    nodes_.push_back(new_pair);
    // Add new incident row
    std::vector<size_type> temprow;
    incident_matrix_.push_back(temprow);
    return Node(this, nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this and n.uid_ < size()) {
      return true;
    } else {
    	return false;
    }
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
      	return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      	return Node(graph_, uid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes in the same graph.
     */
    bool operator==(const Edge& e) const {
    	size_type tempmax = std::max(uid1_, uid2_);
    	size_type tempmin = std::min(uid1_, uid2_);
    	size_type tempmaxe = std::max(e.uid1_, e.uid2_);
    	size_type tempmine = std::min(e.uid1_, e.uid2_);
    	if (graph_ == e.graph_ and tempmax == tempmaxe and tempmin == tempmine) {
    		return true;
    	} else {
    		return false;
    	}
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function obeys trichotomy inside a graph.
     */
    bool operator<(const Edge& e) const {
      	size_type tempmax = std::max(uid1_, uid2_);
    	size_type tempmin = std::min(uid1_, uid2_);
    	size_type tempmaxe = std::max(e.uid1_, e.uid2_);
    	size_type tempmine = std::min(e.uid1_, e.uid2_);
      	if (graph_ < e.graph_) {
        	return true;
      	} else if (graph_ == e.graph_ && (tempmin < tempmine or (tempmin == tempmine and tempmax < tempmaxe))) {
       		return true;
      	} else {
      		return false;
      	}
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // This edge's first node unique identification number
    size_type uid1_;
    // This edge's index to the second node in the incident vector for the first node.
    size_type uid2_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type uid1, size_type uid2)
        : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
        	assert(graph != nullptr);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
  	EdgeIterator ei = edge_begin();
  	for (unsigned int j=0; j<i; ++j) {
  		++ei;
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
  	for (auto ii = a.edge_begin(); ii != a.edge_end(); ++ii) {
  		Edge e = *ii;
  		Node n2 = e.node2();
  		if (n2 == b) {
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
  	if (has_edge(a,b)) {
  		return Edge(this, a.uid_, b.uid_);
  	}

    // Not already in so add
    incident_matrix_[a.uid_].push_back(b.uid_);
    incident_matrix_[b.uid_].push_back(a.uid_);
    //std::pair<size_type,size_type> new_pair (a.uid_, b.uid_);
    //edges_.push_back(new_pair);
    ++num_edges_;
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    //edges_.clear();
    incident_matrix_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
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

   /** Dereference the node iterator.
    * @pre the NodeIterator is not a nullptr
    * @return the node the NodeIterator is pointing to
    */
    Node operator*() const {
    	return Node(graph_, pos_);
    }

   /** Increment the node iterator.
    * @pre the NodeIterator is not the end iterator
    * @return NodeIterator pointer to next node in node vector, or the end iterator
    */ 
    NodeIterator& operator++() {
    	pos_++;
    	return *this;
    }

   /** Test whether this NodeIterator and @a ni are equal.
    * 
    * Equal iterators have the same graph and are at the same position.
    */
    bool operator==(const NodeIterator& ni) const {
    	return (graph_ == ni.graph_ && pos_ == ni.pos_);
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Index of position in node vector
    size_type pos_;
    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type pos)
        : graph_(const_cast<Graph*>(graph)), pos_(pos) {
    }
  };

 /** Return node_iterator pointing to the first node of this graph or the 
  * end iterator if there are no nodes.
  * @return Node_iterator pointing to the first node of this graph, or the end iterator if
  * 	there are no nodes in the graph
  */
  node_iterator node_begin() const {
  	  return NodeIterator(this, 0);
  }

 /** Return node_iterator one past the last node of this graph, if there are any,
  * nodes in the graph.
  * @return Node_iterator one past the last node of this graph
  * @post result must not be able to be dereferenced
  */
  node_iterator node_end() const {
  	  return NodeIterator(this, nodes_.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

   /** Dereference the incident iterator.
    * @pre the IncidentIterator is not a nullptr
    * @return the edge the iterator is pointing to
    * @post result is an edge incident to the node with uid @a uid_
    */
    Edge operator*() const {
    	return Edge(graph_, uid_, graph_->incident_matrix_[uid_][pos_]);
    }

   /** Incremenet the incident iterator.
    * @pre the IncidentIterator is not the end iterator
    * @return IncidentIterator pointer to next node adjacent to node with uid @a uid_, or
    *	the end iterator
    */
    IncidentIterator& operator++() {
    	pos_++;
    	return *this;
    }

   /** Test whether this IncidentIterator and @a ii are equal.
    * 
    * Equal iterators have the same graph, are iterating over the same node, and are at the
    * same position.
    */
    bool operator==(const IncidentIterator& ii) const {
    	return (graph_ == ii.graph_ and uid_ == ii.uid_ and pos_ == ii.pos_);
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Uid of node iterating over
    size_type uid_;
    // Index of position in incidence vector
    size_type pos_;
    /** Private Constructor */
    IncidentIterator(const Graph* graph, size_type uid, size_type pos)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), pos_(pos) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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

   /** Dereference the edge iterator.
    * @pre the iterator is not the end iterator
    * @return the edge the iterator is pointing to
    * @post result is an edge of the graph iterating over
    */ 
    Edge operator*() const {
    	return Edge(graph_, node_pos_, graph_->incident_matrix_[node_pos_][edge_pos_]);
    }

   /** Incremenet the edge iterator.
    * @pre the EdgeIterator is not a nullptr
    * @return EdgeIterator pointer to next edge in graph, or the end iterator
    */ 
    EdgeIterator& operator++() {
    	++edge_pos_;
    	fix();
    	return *this;
    }
    
   /** Test whether this EdgeIterator and @a ei are equal.
    * 
    * Equal iterators have the same graph, are iterating over the same node and are at the
    * same position.
    */ 
    bool operator==(const EdgeIterator& ei) const {
    	if (graph_ == ei.graph_ and node_pos_ == ei.node_pos_ and edge_pos_ == ei.edge_pos_) {
    		return true;
    	} else {
    		return false;
    	}
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Index of node
    size_type node_pos_;
    // Index of position in incident edges to current node
    size_type edge_pos_;
    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type node_pos, size_type edge_pos)
        : graph_(const_cast<Graph*>(graph)), node_pos_(node_pos), edge_pos_(edge_pos) {
        	fix();
    }

    // Helper function to get to next edge with uid1 < uid2
    void fix() {
    	bool done = false;
    	while (!done) {
    		if (node_pos_ == graph_->num_nodes()) {
    			// at end of all edges
    			return;
    		} else if (edge_pos_ == graph_->incident_matrix_[node_pos_].size()) {
    			++node_pos_;
    			edge_pos_ = 0;
    		} else if (node_pos_ > graph_->incident_matrix_[node_pos_][edge_pos_]) {
    			++edge_pos_;
    		} else {
    			done = true;
    		}
    	}
    }
  };

 /** Return edge_iterator pointing to the first edge of this graph or the 
  * end iterator if there are no edges.
  * @return Incident_iterator pointing to the first edge, or the end iterator if
  * 	there are no edges
  */ 
  edge_iterator edge_begin() const {
  	return EdgeIterator(this, 0, 0);
  }
  
 /** Return edge_iterator one past the last edge of this graph.
  * @return edge_iterator one past the last edge of this graph
  * @post result must not be able to be dereferenced
  */ 
  edge_iterator edge_end() const {
  	return EdgeIterator(this, num_nodes(), 0);
  }

 private:

  // vector of node positions and values 
  std::vector<std::pair<Point, node_value_type>> nodes_;
  // number of edges counter
  size_type num_edges_;
  // vector of vectors storing incident edges
  std::vector<std::vector<size_type>> incident_matrix_;

};

#endif // CME212_GRAPH_HPP
