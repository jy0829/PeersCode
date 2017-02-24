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
 private:
  // Internal types are at the bottom

 public:

//
// PUBLIC TYPE DEFINITIONS
//

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

  /** Type of indexes and sizes.
  Return type of Graph::Node::index(), Graph::num_nodes(),
  Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /** Predeclaration of NodeInfo type. */
  struct NodeInfo;
  /** Synonym for NodeInfo */
  using node_info = NodeInfo;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  
  /** Synonym for V **/
  using node_value_type = V;

  /** Synonym for unsigned int */
  using uid_type = size_type;

  /** Synonym for int */
  using idx_type = int;

  /** Predeclaration of EdgeInfo */
  struct EdgeInfo;
  /** Synonym for EdgeInfo */
  using edge_info = EdgeInfo;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Synonym for E **/
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



  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
      : nodes_(), i2u_() {
      	num_edges_ = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::NodeInfo
   * @brief Class representing the graph's node representations
   */
  struct NodeInfo {
  	Point p_;
  	node_value_type value_;
  	std::vector<edge_info> adj_;
  	idx_type idx_;
  	NodeInfo(Point p,node_value_type value,std::vector<edge_info> adj,idx_type idx) {
  		p_ = p;
  		value_ = value;
  		adj_ = adj;
  		idx_ = idx;
  	}
  };

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
      return graph_->nodes_[uid_].p_;
    }

    /** Return a reference to this node's position */
    Point& position() {
    	return graph_->nodes_[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[uid_].idx_;
    }

   /** Return this node's value, a value of type node_value_type.
    * @return Value of the node.
    */
    node_value_type& value() {
    	return graph_->nodes_[uid_].value_;
    }

   /** Return this node's value, a value of type node_value_type.
    * For when called on a constant type.
    * @return Value of the node.
    */
    const node_value_type& value() const {
    	return graph_->nodes_[uid_].value_;
    }

   /** Return the degree of this node.
    * @return Number of edges incident to this node.
    */
    size_type degree() const {
    	return graph_->nodes_[uid_].adj_.size();
    }

   /** Return the uid of this node.
    * @return The uid of this node
    */
    size_type uid() const {
    	return uid_;
    } 

   /** Return incident_iterator pointing to the first edge incident to this node or the 
    * end iterator if there are no incident edges
    * @return Incident_iterator pointing to the first incident edge, or the end iterator if
    * 	there are no edges incident to this node.
    */
    incident_iterator edge_begin() const {
    	return incident_iterator(graph_, uid_, 0);
    }
    incident_iterator begin() const {
    	return edge_begin();
    }

   /** Return incident_iterator one past the last incident edge, if there are any,
    * of this node.
    * @return Incident_iterator one past the last incident edge of this node.
    * @post result must not be able to be dereferenced
    */
    incident_iterator edge_end() const {
    	return incident_iterator(graph_, uid_, degree());
    }
    incident_iterator end() const {
    	return edge_end();
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
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        	assert(graph != nullptr);
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
   * @param[in] value The new node's value. Default is node_value_type holder.
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    std::vector<edge_info> new_adj;
    node_info new_node_info(position, value, new_adj, i2u_.size());
    nodes_.push_back(new_node_info);
    // Update i2u
    i2u_.push_back(nodes_.size()-1);
    return Node(this, nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this and nodes_[n.uid_].idx_ != -1 and nodes_[n.uid_].idx_ < (int)size()) {
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
  Node node(idx_type i) const {
    return Node(this, i2u_[i]);
  }

  /** Return the index of node with uid @a uid.
   *
   */
  idx_type uid2idx(size_type uid) const {
  	return nodes_[uid].idx_;
  }

  //
  // EDGES
  //

  /** @class Graph::EdgeInfo
   * @brief Class representing the edge representation
   */
  struct EdgeInfo {
  	uid_type uid_other_;
  	edge_value_type value_;
  	EdgeInfo(uid_type uid_other,edge_value_type value) {
  		uid_other_ = uid_other;
  		value_ = value;
  	}
  };

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
    	size_type tempmaxe = std::max(e.uid1_, uid2_);
    	size_type tempmine = std::min(e.uid1_, uid2_);
      	if (graph_ < e.graph_) {
        	return true;
      	} else if (graph_ == e.graph_ && (tempmin < tempmine or (tempmin == tempmine and tempmax < tempmaxe))) {
       		return true;
      	} else {
      		return false;
      	}
    }

    /** Return the edges length
     */
    double length() const {
    	return norm(node1().position() - node2().position());
    }

    /** Return the edges value
    */
    edge_value_type& value() {
    	/*Edge opposite(graph_, node2().uid_, node1().uid_);
    	graph_->nodes_[uid2_].adj_[opposite.inc_id()].value_ = graph_->nodes_[uid1_].adj_[inc_id()].value_;*/

    	return graph_->nodes_[uid1_].adj_[inc_id()].value_;
    }

    const edge_value_type& value() const {
    	return graph_->nodes_[uid1_].adj_[inc_id()].value_;
    }

    /*void set_opp_value(edge_value_type value) {
    	graph_->nodes_[uid1_].adj_[inc_id()].value_ = value;
    }*/

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // This edge's first node unique identification number.
    size_type uid1_;
    // This edge's second node unique identification number.
    size_type uid2_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type uid1, size_type uid2)
        : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
        	assert(graph != nullptr);
    }
    /* OLD Helper function to get second node's uid */
    /*size_type second_uid() const {
    	return graph_->nodes_[uid1_].adj_[inc_id_].uid_other_;
    }*/
    /* Helper function to get inc_id */
    size_type inc_id() const {
    	int ctr = 0;
    	for (auto e : node1()) {
    		if (e.node2() == node2()) {
    			return ctr;
    		} else {
    			ctr++;
    		}
    	}
    	return node1().degree();
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
  	for (auto e : a) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()){
  	// Check if already in the graph
  	if (has_edge(a,b)) {
  		return Edge(this, a.uid_, b.uid_);
  	}

    // Not already in so add
    edge_info new_edge_info(b.uid_,value);
    nodes_[a.uid_].adj_.push_back(new_edge_info);

    edge_info new_edge_info2(a.uid_,value);
    nodes_[b.uid_].adj_.push_back(new_edge_info2);

    ++num_edges_;
    return Edge(this, a.uid_, b.uid_);
  }

  // REMOVE METHODS

  /** Remove node @a n from the graph, if it is present.
   * @param[in] n Node to be removed.
   * @return 1 if node succesfully removed, 0 if node not in graph.
   * 
   * @post new num_nodes() = old num_nodes() - 1, if 1 returned
   * @post new num_edges = old num_edges - n.degree(), if 1 returned
   * @post If the last Node is being removed all elements index are not affected.
   			Otherwise the Node with position num_nodes()-1 is moved to position old n.index(), and all other elements are not affected.
   *
   * @note Invalidates Node that @a n_it is pointing to. All other nodes remain validated.
   * @note Invalidates all Edges incident to the Node that @a n_it is pointing to. All other 		edges remain validated. 
   * @note Invalidates NodeIterator node_end()-1 and any NodeIterator pointing to @a n. All 		other NodeIterators remain valid.
   * @note Invalidates all EdgeIterators and IncidentIterators incident to @a n, or incident to a Node adjacent to @a n. All other EdgeIterators and IncidentIterators remain valid.
   * Complexity: O(1) (assuming max degree squared of a node is O(1)).
   */
  size_type remove_node(const Node& n) {
  	if (has_node(n)) {
  		incident_iterator start = n.edge_begin();
  		while (start != n.end()) {
  			start = remove_edge(start);
  		}
  		nodes_[i2u_.back()].idx_ = nodes_[n.uid_].idx_;
  		i2u_[nodes_[n.uid_].idx_] = i2u_.back();
  		i2u_.pop_back();
  		nodes_[n.uid_].idx_ = -1;
  		return 1;
  	}
  	return 0;
  }

  /** Remove from the graph the node that NodeIterator @a n_it is pointing to.
   * @param[in] n_it NodeIterator pointing to node to be removed.
   * @return If removing the last Node return the end NodeIterator
   			Otherwise return a NodeIterator pointing to the new location of the Node that was last in the graph.
   * 
   * @post new num_nodes() = old num_nodes() - 1, if 1 returned
   * @post new num_edges = old num_edges - n.degree(), if 1 returned
   * @post If the last Node is being removed all elements index are not affected.
   			Otherwise the Node with index num_nodes()-1 is moved to index old n.index(), and all other elements are not affected.
   *
   * @note Invalidates Node that @a n_it is pointing to. All other nodes remain validated.
   * @note Invalidates all Edges incident to the Node that @a n_it is pointing to. All other 		edges remain validated. 
   * @note Invalidates NodeIterator node_end()-1 and any NodeIterator pointing to @a n. All 		other NodeIterators remain valid.
   * @note Invalidates all EdgeIterators and IncidentIterators incident to @a n, or incident to a Node adjacent to @a n. All other EdgeIterators and IncidentIterators remain valid.
   * Complexity: O(1) (assuming max degree squared of a node is O(1)).
   */
  node_iterator remove_node(node_iterator n_it) {
 	remove_node(*n_it);
 	return n_it;
  }

  /** Remove edge with nodes @a n1 and @a n2 from graph.
   * @param[in] n1 First node of the edge to be removed from graph.
   * @param[in] n2 Second node of the edge to be removed from graph.
   * @return 1 if succesfully removed, 0 if not in graph.
   *
   * @post new num_edges() = old num_edges() - 1
   * @post new n1.degree() = old n1.degree() - 1
   * @post new n2.degree() = old n2.degree() - 1
   *
   * @note All Nodes remain valid.
   * @note Invalidates the Edge being removed. All other Edges remain valid.
   * @note All NodeIterators remain valid.
   * @note Invalidates all EdgeIterators and IncidentIterators pointing to an edge incident to 			either @a n1 or @a n2. All other EdgeIterators and IncidentIterators remain valid.
   * Complexity: O(1) (assuming max degree of edges if O(1)).
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
  	// check has edge
  	if (has_edge(n1, n2)) {
  		size_type inc_id1 = Edge(this, n1.uid_, n2.uid_).inc_id();
  		size_type inc_id2 = Edge(this, n2.uid_, n1.uid_).inc_id();

  		// swap back element for back
  		nodes_[n1.uid_].adj_[inc_id1] = nodes_[n1.uid_].adj_.back();
  		nodes_[n2.uid_].adj_[inc_id2] = nodes_[n2.uid_].adj_.back();

  		// pop back
  		nodes_[n1.uid_].adj_.pop_back();
  		nodes_[n2.uid_].adj_.pop_back();
  		num_edges_--;
  		return 1;
  	}
  	return 0;
  }

  /** Remove edge @a e from graph.
   * @param[in] e Edge to be removed from graph.
   * @return 1 if succesfully removed, 0 if not in graph.
   *
   * @post new num_edges() = old num_edges() - 1
   * @post new e.node1().degree() = old e.node1().degree() - 1
   * @post new e.node2().degree() = old e.node2().degree() - 1
   *
   * @note All Nodes remain valid.
   * @note Invalidates the Edge being removed. All other Edges remain valid.
   * @note All NodeIterators remain valid.
   * @note Invalidates all EdgeIterators and IncidentIterators pointing to an edge incident to 			either @a n1 or @a n2. All other EdgeIterators and IncidentIterators remain valid.
   * Complexity: O(1) (assuming max degree of edges if O(1)).
   */
  size_type remove_edge(const Edge& e) {
  	return remove_edge(e.node1(), e.node2());
  }

  /** Remove edge e that IncidentIterator @a i_it is pointing to from graph.
   * @param[in] i_it IncidentIterator pointing to Edge to be removed from graph.
   * @return If @a i_it is pointing to last Edge of the adjacency vector for e.node1(),return 			the end IncidentIterator for e.node1().
   			Otherwise return a IncidentIterator pointing to the new location of the Edge that was last in the adjacency vector for e.node1().
   *
   * @post new num_edges() = old num_edges() - 1
   * @post new node1().degree() = old node1().degree() - 1
   * @post new node2().degree() = old node2().degree() - 1
   *
   * @note All Nodes remain valid.
   * @note Invalidates the Edge being removed. All other Edges remain valid.
   * @note All NodeIterators remain valid.
   * @note Invalidates all EdgeIterators and IncidentIterators pointing to an edge incident to 			either @a n1 or @a n2. All other EdgeIterators and IncidentIterators remain valid.
   * Complexity: O(1) (assuming max degree of edges if O(1)).
   */
  incident_iterator remove_edge(incident_iterator i_it) {
  	remove_edge(*i_it);
  	return i_it;
  }

  /** Remove edge e that EdgeIterator @a e_it is pointing to from graph.
   * @param[in] e_it EdgeIterator pointing to Edge to be removed from graph.
   * @return If @a e_it is edge_end()-1 return edge_end().
			Else if @a e_it points to the same edge as node1().edge_end()-1 or the edge that node1().edge_end()-1 points to cannot be pointed to by a valid EdgeIterator, return the next valid EdgeIterator after @a e_it.
			Otherwise return the EdgeIterator node1().edge_end()-1.
   *
   * @post new num_edges() = old num_edges() - 1
   * @post new node1().degree() = old node1().degree() - 1
   * @post new node2().degree() = old node2().degree() - 1
   *
   * @note All Nodes remain valid.
   * @note Invalidates the Edge being removed. All other Edges remain valid.
   * @note All NodeIterators remain valid.
   * @note Invalidates all EdgeIterators and IncidentIterators pointing to an edge incident to 			either @a n1 or @a n2. All other EdgeIterators and IncidentIterators remain valid.
   * Complexity: O(1) (assuming max degree of edges if O(1)).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
  	remove_edge(*e_it);
  	e_it.fix();
  	return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    i2u_.clear();
    num_edges_ = 0;
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
    	return Node(graph_, graph_->i2u_[pos_]);
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
  	  return NodeIterator(this, i2u_.size());
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
    	return Edge(graph_, uid_, graph_->nodes_[uid_].adj_[pos_].uid_other_);
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
    	return Edge(graph_, graph_->i2u_[node_pos_], graph_->nodes_[graph_->i2u_[node_pos_]].adj_[edge_pos_].uid_other_);
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
    		} else if (edge_pos_ == graph_->nodes_[graph_->i2u_[node_pos_]].adj_.size()) {
    			// at end of incident edges for this node
    			++node_pos_;
    			edge_pos_ = 0;
    		} else if (graph_->i2u_[node_pos_] > graph_->nodes_[graph_->i2u_[node_pos_]].adj_[edge_pos_].uid_other_) {
    			// Not good, try next one
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
  std::vector<node_info> nodes_;
  // number of edges counter
  size_type num_edges_;
  // vector mapping uid to idx
  std::vector<uid_type> i2u_;
  // vector of vectors storing incident edges
  //std::vector<std::vector<size_type>> incident_matrix_;

};

// TO USE AUTO WE DEFINE THE FOLLOWING
template<typename V, typename E>
struct NodeRange {
	Graph<V, E>& g_;
	using GraphType = Graph<V, E>;
	using NodeIterator = typename GraphType::node_iterator;
	NodeIterator begin() {
		return g_.node_begin();
	}
	NodeIterator end() {
		return g_.node_end();
	}
};

template<typename V,typename E>
NodeRange<V, E> nodesof(Graph<V, E>& g) {
	return NodeRange<V, E> {g};
}

template<typename V, typename E>
struct EdgeRange {
	Graph<V, E> &g_;
	using GraphType = Graph<V, E>;
	using EdgeIterator = typename GraphType::edge_iterator;
	EdgeIterator begin() {
		return g_.edge_begin();
	}
	EdgeIterator end() {
		return g_.edge_end();
	}
};

template<typename V, typename E>
EdgeRange<V,E> edgesof(Graph<V, E>& g) {
	return EdgeRange<V, E> {g};
}


#endif // CME212_GRAPH_HPP
