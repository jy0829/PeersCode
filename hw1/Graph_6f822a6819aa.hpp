#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP
#define NDEBUG

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

	class InternalNode;

 public:

	//
	// PUBLIC TYPE DEFINITIONS
	//

	/** Type of this graph. */
	using graph_type = Graph;

	/** Synonym for the type of node value */
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
	
	//-----------------------------//
	// CONSTRUCTORS AND DESTRUCTOR //
	//-----------------------------//

	/** Construct an empty graph. */
	Graph() {
	}

	/** Default destructor */
	~Graph() = default;

	//-------//
  	// NODES //
  	//-------//

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
			// Do nothing to constuct an invalid node.
		}

		/** Return this node's position. */
		const Point& position() const {
			return fetch().point;
		}

		/** Return this node's index, a number in the range [0, graph_size). */
		size_type index() const {
			return uid_;
		}

		/** Test whether this node and @a n are equal.
		 *
		 * Equal nodes have the same graph and the same index.
		 */
		bool operator==(const Node& n) const {
			if (n.graph_== graph_ && n.uid_ == uid_) {
				return true;
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
			return (n.uid_ < uid_);
		}

		// Return the node's value
		node_value_type& value() {
			return fetch().value;
		}

		// Return the node's value (const version)
		const node_value_type& value() const {
			return fetch().value;
		}

		// Returns the number of neighbors belonging to the node
		size_type degree() const {
			return fetch().neighbors.size();
		}

		// Returns the beginning incident iterator 
		incident_iterator edge_begin() const {
			return IncidentIterator(graph_, uid_, 0);
		}

		// Returns the ending incident iterator
		incident_iterator edge_end() const {
			return IncidentIterator(graph_, uid_, degree());
		}

	 private:

		// Allow Graph to access Node's private member data and functions.
		friend class Graph;
		
		graph_type* graph_;
		size_type uid_;

		// Private constructor
		Node(const graph_type* graph, size_type uid)
			: graph_(const_cast<graph_type*>(graph)), uid_(uid) {

		}

		/* Helper method to return the appropriate element
		 * It loops through the elements until the element with the correct
		 * index is found.
		 */
		InternalNode& fetch() const {
			assert(!(uid_ > graph_->size()));
			return graph_->nodes[uid_]; 
		}

	};

	/** Return the number of nodes in the graph.
	 *
	 * Complexity: O(1).
	 */
	size_type size() const {
		return nodes.size();
	}

	/** Synonym for size(). */
	size_type num_nodes() const {
		return size();
	}

	/** Add a node to the graph, returning the added node.
	 * @param[in] position The new node's position
	 * @param[in] value_in The value associated with this node
	 * @post new num_nodes() == old num_nodes() + 1
	 * @post result_node.index() == old num_nodes()
	 *
	 * Complexity: O(1) amortized operations.
	 */
	Node add_node(const Point& position, const node_value_type& value_in = node_value_type()) {
		InternalNode node = InternalNode(position, nodes.size(), value_in);
		nodes.push_back(node);
		return Node(this, nodes.size()-1);
	}

	/** Determine if a Node belongs to this Graph
	 * @return True if @a n is currently a Node of this Graph
	 *
	 * Complexity: O(1).
	 */
	bool has_node(const Node& n) const {
		// HW0: YOUR CODE HERE
		InternalNode node;
		node.uid = n.index();
		// Try to find node in vector of nodes
		auto p = std::find(nodes.begin(), nodes.end(), node);
		if (p == nodes.end()) {
			return false;
		}
		return true;
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

	//-------//
	// EDGES //
	//-------//

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
			// Do nothing to construct an invalid edge.
		}

		/** Return a node of this Edge */
		Node node1() const {
			return graph_->node(uid1_);
		}

		/** Return the other node of this Edge */
		Node node2() const {
			return graph_->node(uid2_);
		}

		/** Test whether this edge and @a e are equal.
		 *
		 * Equal edges represent the same undirected edge between two nodes.
		 */
		bool operator==(const Edge& e) const {
			if (graph_ == e.graph_ && uid1_ == e.uid1_ && uid2_ == e.uid2_)
				return true;
			if (graph_ == e.graph_ && uid1_ == e.uid2_ && uid2_ == e.uid1_)
				return true;
			return false;
		}

		/** Test whether this edge is less than @a e in a global order.
		 *
		 * This ordering function is useful for STL containers such as
		 * std::map<>. It need not have any interpretive meaning.
		 */
		bool operator<(const Edge& e) const {
			if (graph_ == e.graph_)
				return uid1_ < e.uid1_;
			return graph_ < e.graph_;
		}

	 private:
		// Allow Graph to access Edge's private member data and functions.
		friend class Graph;

		// Private instance variables
		graph_type* graph_;
		size_type uid1_;
		size_type uid2_;

		// Private constructor
		Edge(const graph_type* graph, size_type uid1, size_type uid2)
			: graph_(const_cast<graph_type*> (graph)), uid1_(uid1), uid2_(uid2) {
				assert(graph_ != nullptr);
				assert(uid1_ != uid2_);
			}
	};

	/** Return the total number of edges in the graph.
	 *
	 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	 */
	size_type num_edges() const {
		return numedges_;
	}

	/** Return the edge with index @a i.
	 * @pre 0 <= @a i < num_edges()
	 *
	 * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
	 */
	Edge edge(size_type i) const {
		EdgeIterator ei = edge_begin();
		while (i > 0) {
			--i;
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
		// Try to find b in the vector containing the neighbors of a
		size_type a_uid = a.index();
		size_type b_uid = b.index();
		auto p = std::find(nodes[a_uid].neighbors.begin(), nodes[a_uid].neighbors.end(), b_uid);
		if (p == nodes[a_uid].neighbors.end())
			return false;
		return true;
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
		// Check if there is already an edge between these two
		size_type a_uid = a.index();
		size_type b_uid = b.index();
		if (has_edge(a, b)) {
			if (a_uid < b_uid)
				return Edge(this, a_uid, b_uid);
			else
				return Edge(this, b_uid, a_uid);
		}
		// If it doesn't already exist, add a new edge
		nodes[a_uid].neighbors.push_back(b_uid);
		nodes[b_uid].neighbors.push_back(a_uid);
		numedges_++;
		if (a_uid < b_uid)
			return Edge(this, a_uid, b_uid);
		else
			return Edge(this, b_uid, a_uid);
 	}

	/** Remove all nodes and edges from this graph.
	 * @post num_nodes() == 0 && num_edges() == 0
	 *
	 * Invalidates all outstanding Node and Edge objects.
	 */
	void clear() {
		numedges_ = 0;
		nodes.clear();
	}

	//---------------//
	// Node Iterator //
	//---------------//

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

		/** Return the node that the iterator points to */
		Node operator*() const {
			return graph_->node(uid_);
		}

		/** Return the node iterator that the next node
		 * @pre uid_ < graph->size()
		 * @post result->index() = uid_ + 1
		 */
		NodeIterator& operator++() {
			++uid_;
			return *this;
		}
		
		/** Return true if this node iterator is the same as
		 * @a ni. They are the same if they are part of the same
		 * graph and they point to the same object.
		 */
		bool operator==(const NodeIterator& ni) const {
			if (graph_ == ni.graph_ && uid_ == ni.uid_)
				return true;
			return false;
		}

	 private:
		friend class Graph;

		// Private instance variables
		graph_type* graph_;
		size_type uid_;

		// Private constructor
		NodeIterator(const graph_type* graph, size_type uid)
		: graph_(const_cast<graph_type*>(graph)), uid_(uid) {
		}

	};

	/** Return the node iterator that points to the first node in the graph */
	node_iterator node_begin() const {
		return NodeIterator(this, 0);
	}

	/** Return a node iterator that points to just beyond the last node in the graph
	 * This node iterator should not be dereferenced because it does not point 
	 * to an actual node.
	 */
	node_iterator node_end() const {
		return NodeIterator(this, nodes.size() + 1);
	}

	//-------------------//
	// Incident Iterator //
	//-------------------//

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

		/** Construct and return the edge that the iterator points to */
		Edge operator*() const {
			return Edge(graph_, uid_, graph_->nodes[uid_].neighbors[edge_idx_]);
		}
		
		/** Return the incident iterator that points to the next edge */
		IncidentIterator& operator++() {
			if (edge_idx_ == graph_->nodes[uid_].neighbors.size()) {
				return *this;
			} else {
				edge_idx_++;
				return *this;
			} 
		}
		
		/** Test whether two incident iterators are part of the same graph,
		 * point to the same edge, and are part of the same node
		 */
		bool operator==(const IncidentIterator& ii) const {
			return (graph_ == ii.graph_ && uid_ == ii.uid_ && edge_idx_ == ii.edge_idx_);
		}

	 private:
		friend class Graph;
		
		graph_type* graph_;
		size_type uid_;
		size_type edge_idx_;

		// Private constructor
		IncidentIterator(const graph_type* graph, size_type uid, size_type edge_idx)
		: graph_(const_cast<graph_type*>(graph)), uid_(uid), edge_idx_(edge_idx) {
		}

		};

	//---------------//
	// Edge Iterator //
	//---------------//

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

		/** Return the edge constructed by the current state of the edge iterator */
		Edge operator*() const {
			return Edge(graph_, node_, graph_->nodes[node_].neighbors[neighbor_uid_]);
		}

		/** Return the edge iterator that points to the next edge.
		 * @pre has_edge((*this).node1(), (*this).node2())
		 * @post (*this).node1() < (*this).node1()
		 */
		EdgeIterator& operator++() {
			do {
				if (neighbor_uid_ < (graph_->nodes[node_].neighbors.size() - 1)) {
					neighbor_uid_++;
				} else {
					neighbor_uid_ = 0;
					node_++;
					while (node_ != graph_->nodes.size() &&
						graph_->nodes[node_].neighbors.size() == 0) {
						node_++;
					}
				}
			} while (node_ != graph_->nodes.size() && 
				graph_->nodes[node_].neighbors[neighbor_uid_] < node_);

			return *this;
		}

		/** Test if two edge iterators are the same
		 * The are the same if they are part of the same graph,
		 * reference the same node, and reference the same neighbor
		 */
		bool operator==(const EdgeIterator& ei) const {
			return (graph_ == ei.graph_ && node_ == ei.node_ && neighbor_uid_ == ei.neighbor_uid_);
		}

	 private:
		friend class Graph;

		graph_type* graph_;
		size_type node_;
		size_type neighbor_uid_;

		// Private constructor
		EdgeIterator(const graph_type* graph, size_type node, size_type neighbor_uid)
			: graph_(const_cast<graph_type*>(graph)), node_(node), neighbor_uid_(neighbor_uid) {
			}
	};

	/** Return the edge iterator that points to the graph's first edge */
	edge_iterator edge_begin() const {
		unsigned int i = 0;
		while (i < nodes.size()) {
			if (nodes[i].neighbors.size() != 0) {
				return EdgeIterator(this, i, 0);
			}
			i++;
		}
		return EdgeIterator(this, size(), 0);
	}

	/** Return the edge iterator that points to just beyond the last edge 
	 * of the graph.  This edge iterator should not be dereferenced.
	 */
	edge_iterator edge_end() const {
		return EdgeIterator(this, nodes.size(), 0);
	}

 private:  

	class InternalNode {
	 public:

	 	InternalNode(Point point, size_type uid, node_value_type value)
	 	: point(point), uid(uid), value(value) {
	 	}

		Point point;
		size_type uid;
		node_value_type value;
		std::vector<size_type> neighbors;

		// Need to define the == operator for nodes
		// for use in has_node()
		bool operator==(const InternalNode& n) const {
			if (uid == n.uid) {
				return true;
			}
			return false;
		}
	};

	size_type numedges_;
	std::vector<InternalNode> nodes;

};

#endif // CME212_GRAPH_HPP
