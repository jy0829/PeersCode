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


/** min and max inline function to compare indices */
inline unsigned min(unsigned x, unsigned y) {
	return x < y ? x : y; 
}
inline unsigned max(unsigned x, unsigned y) { 
	return x > y ? x : y; 
}

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph{
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
  using graph_type = Graph<V,E>;

	/** Type of the node value */
	using node_value_type = V;

	/** Type of the edge value */
	using edge_value_type = E;

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
  Graph() {
   // HW0: YOUR CODE HERE
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
			return graph_->nodes_[uid_].position_;
    }

    /** Return this node's position. */
    Point& position() {
			return graph_->nodes_[uid_].position_;
		}

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].i_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

		/** Returns a reference on the value of a valid node. */
		node_value_type& value(){
			return graph_->nodes_[uid_].v_;
		}

		/** Returns a constant reference on the value of a valid node. */
		const node_value_type& value() const{
			return graph_->nodes_[uid_].v_;
		}

		/** Returns the degree of a valid node. */
		size_type degree() const{
			return graph_->nodes_[uid_].incident_nodes_.size();
		}

		/** Returns an incident iterator set to the first adjacent edge of       * node object. */
		incident_iterator edge_begin() const{
			return IncidentIterator(graph_, uid_, 0);
		}

		/** Returns an incident iterator set to the last adjacent edge of       * node object. */
		incident_iterator edge_end() const{
			return IncidentIterator(graph_, uid_, graph_->nodes_[uid_].incident_nodes_.size());
		}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const{
      // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
			if(graph_ == n.graph_)
				return uid_ < n.uid_;
			else
				return graph_ < n.graph_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(graph_type *graph_, size_type uid_):
					graph_(graph_), uid_(uid_) {}
    
    graph_type *graph_;
    size_type uid_;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
		nodes_.push_back(nodeinfo(i2u_.size(), position, value));
		i2u_.push_back(nodes_.size()-1);
    return Node(this, nodes_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.graph_ == this and nodes_[n.uid_].i_ < i2u_.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(const_cast<graph_type*>(this), i2u_[i]);
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
			return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
		 return Node(graph_, uid2_);
    }

		double length() const {
			return norm(graph_->nodes_[uid1_].position_ - graph_->nodes_[uid2_].position_);
		}

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph_==e.graph_ 
				and min(uid1_, uid2_)==min(e.uid1_, e.uid2_)
				and max(uid1_, uid2_)==max(e.uid1_, e.uid2_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
			if(graph_== e.graph_){

				size_type m=min(uid1_, uid2_), M=max(uid1_, uid2_),
				  em=min(e.uid1_, e.uid2_), eM=max(e.uid1_, e.uid2_);

				if(m == em)
					return M < eM;
				else
					return m < em;
			}
			else
				return graph_ < e.graph_;
    }

		edge_value_type& value() {
				size_type uidm=min(uid1_, uid2_), uidM=max(uid1_, uid2_);

				for(size_type j=0 ; j<graph_->nodes_[uidm].incident_nodes_.size() ; ++j)
					if(graph_->nodes_[uidm].incident_nodes_[j].uid_ == uidM)
						return graph_->nodes_[uidm].incident_nodes_[j].v_;	
				assert(false);
		}

		const edge_value_type& value() const {
			size_type uidm=min(uid1_, uid2_), uidM=max(uid1_, uid2_);

			for(size_type j=0 ; j<graph_->nodes_[uidm].incident_nodes_.size() ; ++j)
				if(graph_->nodes_[uidm].incident_nodes_[j].uid_ == uidM)
					return graph_->nodes_[uidm].incident_nodes_[j].v_;	
		}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  Edge(graph_type *graph_, size_type uid1_, size_type uid2_) :
    graph_(graph_), uid1_(uid1_), uid2_(uid2_) {}

    graph_type *graph_;
    size_type uid1_, uid2_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
		size_type s=0;
		for(size_type i=0 ; i<num_nodes() ; ++i)
			s+=nodes_[i2u_[i]].incident_nodes_.size();

		assert(s%2==0);
		return s/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
		return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for(size_type i=0 ; i<nodes_[a.uid_].incident_nodes_.size() ; ++i)
			if(nodes_[a.uid_].incident_nodes_[i].uid_ == b.uid_)
				return true;
			
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
    if(has_edge(a, b))
			return Edge(this, a.uid_, b.uid_);

		nodes_[a.uid_].incident_nodes_.push_back(edgeinfo(b.uid_));
		nodes_[b.uid_].incident_nodes_.push_back(edgeinfo(a.uid_));

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
    i2u_.clear();
  }

	/** Removes Edge(@a node1, @a node2) from the graph.
		* @pre Edge(@a node1, @a node2) is a valid edge.
		* @post new num_edges() = old num_edges()-1
		* @return new num_edges()
		* Complexity O(no
		*/
	size_type remove_edge(const Node& node1, const Node& node2){

		for(size_type j=0 ; j<nodes_[node1.uid_].incident_nodes_.size() ; ++j)
			if(nodes_[node1.uid_].incident_nodes_[j].uid_ == node2.uid_){
				nodes_[node1.uid_].incident_nodes_.erase(
						nodes_[node1.uid_].incident_nodes_.begin()+j);
				break;
			}
		
		for(size_type j=0 ; j<nodes_[node2.uid_].incident_nodes_.size() ; ++j)
			if(nodes_[node2.uid_].incident_nodes_[j].uid_ == node1.uid_){
				nodes_[node2.uid_].incident_nodes_.erase(
						nodes_[node2.uid_].incident_nodes_.begin()+j);
				break;
			}
		return num_edges();
	}

	/** Removes @a e from the graph.
		* @pre @a e is a valid edge.
		* @post new num_edges() = old num_edges()-1
		* @return new num_edges()
		*/
	size_type remove_edge(const Edge& e){
		return remove_edge(e.node1(), e.node2());
	}

	/** Removes @a *e_it from the graph.
		* @pre @a e_it is a valid edge iterator and e_it != edge_end()
		* @post new num_edges() = old num_edges()-1
		* @return an edge_iterator pointing at the new successor 
		* of the former edge. If the former edge was the last edge,
		* returns edge_end()
		*/
	edge_iterator remove_edge(edge_iterator e_it){
		remove_edge(e_it.node1(), e_it.node2());
		return e_it;
	}

	/** Removes @a n and all incident edges from this graph.
		* @pre n is a valid node.
		* @post new size() = old size()-1
		* @post All edges containing @a n are removed
		* @return new size()
		* Complexity: O(n.degree())
		*/
	size_type remove_node(const Node& n){
		while(nodes_[n.uid_].incident_nodes_.size()){
			remove_edge(n,
				 	Node(this, nodes_[n.uid_].incident_nodes_.back().uid_));
		}
		
		nodes_[i2u_[num_nodes()-1]].i_=nodes_[n.uid_].i_;
		i2u_[nodes_[n.uid_].i_]=i2u_[num_nodes()-1];
		i2u_.pop_back();
		return num_nodes();
	}

	/** Removes @a *n_it and all incident edges from this graph.
		* @pre n_it != node_end()
		* @post new size() = old size()-1
		* @post All edges containing @a *n_it are removed
		* @return an iterator pointing at the new successor of the former node
		* If the former node was the last node, returns node_end()
		* Complexity: O(n.degree())
		*/
	node_iterator remove_node(node_iterator n_it){
		remove_node(*n_it);
		return n_it;
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

		/** Returns the node currently referred to by the iterator
		  * @pre 0 <= j_ < num_nodes() 
		  */
		Node operator*() const{
			return graph_->node(j_);
		}

		/** Increments the iterator
			* @pre 0<= j_ < num_node()
			*/
		NodeIterator& operator++(){
			++j_;
			return *this;
		}

		/** Test equality for two node iterators
			* Two node iterators are equal if they have same @a graph_ and same @a j_
			*/
		bool operator==(const node_iterator& nit) const{
			return graph_==nit.graph_ and j_==nit.j_;
		}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
		NodeIterator(graph_type *graph_, size_type j_) :
			graph_(graph_), j_(j_) {}

		graph_type *graph_;
		size_type j_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

	/** Returns a node iterator set at the first node of the graph. */
	node_iterator node_begin() const{
		return NodeIterator(const_cast<graph_type*>(this), 0);
	}

	/** Returns a node iterator set at the last node of the graph. */
	node_iterator node_end() const{
		return NodeIterator(const_cast<graph_type*>(this), i2u_.size());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

		/** Returns the edge to which the iterator is referrring to.
			* @pre 0 <= uid_ < num_nodes()
			* @pre 0 <= j_ < @a edges_[uid].size()
			*/
		Edge operator*() const{
	  	return Edge(graph_, uid_, graph_->nodes_[uid_].incident_nodes_[j_].uid_);
		}

		/** Increments iterator to refer to the next adjacent edges of node node(uid_)
			* @pre 0 <= uid_ < num_nodes()
			* @pre 0 <= j_ < @a edges_[uid].size()
			*/
		incident_iterator& operator++(){
			++j_;
			return *this;
		}

		/** Test two iterators for equality
			* Two incident iterators are said equal if they have same:
			* @a graph_, @a, @a uid_ and @a j_
			*/
		bool operator==(const incident_iterator& iit) const{
			return graph_==iit.graph_ and uid_==iit.uid_ and j_==iit.j_;
		}

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
		IncidentIterator(graph_type* graph_, size_type uid_, size_type j_)
		 	:	graph_(graph_), uid_(uid_), j_(j_) {}

		graph_type* graph_;
		size_type uid_;
		size_type j_;
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

		/** Returns the edge referred to by the edge iterator
			* @pre 0 <= uid_ < num_nodes()
			* @pre 0 <= j_ < @a graph_->edges_[uid_].size()
			*/
		Edge operator*() const{
			return Edge(graph_, graph_->i2u_[i_],
					graph_->nodes_[graph_->i2u_[i_]].incident_nodes_[j_].uid_);
		}

		/** Increment the edge iterator to refer to the next edge
			* @pre 0 <= i_ < num_nodes()
			* @pre 0 <= j_ < @a graph_->nodes_[graph_->i2u_[i_]].incident_nodes_.size()
			*/
		EdgeIterator& operator++(){
			if(++j_==graph_->nodes_[graph_->i2u_[i_]].incident_nodes_.size()){
				j_=0;
				++i_;
			}

			if(i_==graph_->size()
					or graph_->i2u_[i_] < graph_->nodes_[graph_->i2u_[i_]].incident_nodes_[j_].uid_)
				return *this;
			else
				return ++(*this);
		}

		/** Tests two iterators for equality
			* Two iterators are said equal if they have the same:
			* @a graph_, @a uid_ and @a j_
			*/
		bool operator==(const EdgeIterator& eit) const{
			return graph_==eit.graph_ and i_==eit.i_ and j_==eit.j_;
		}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
		EdgeIterator(graph_type* graph_, size_type i_, size_type j_) :
			graph_(graph_), i_(i_), j_(j_) {}

		graph_type* graph_;
		size_type i_, j_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

	/** Returns an edge iterator set to the first edge of the graph. */
  edge_iterator edge_begin() const{
		size_type i=-1;
		while(++i < i2u_.size() and nodes_[i].incident_nodes_.size()==0);
		return EdgeIterator(const_cast<graph_type*>(this), i, 0);
	}

	/** Returns an edge iterator set to the end of the graph. */
	edge_iterator edge_end() const{
		return EdgeIterator(const_cast<graph_type*>(this),
				i2u_.size(), 0);
	}

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
	struct edgeinfo {
		edgeinfo(size_type uid_) : uid_(uid_) {}

		size_type uid_;
		edge_value_type v_;
	};

	struct nodeinfo {
		nodeinfo (size_type i_, Point position_, node_value_type v_) :
			i_(i_), position_(position_), v_(v_) {}

		size_type i_;
		Point position_;
		node_value_type v_;
		std::vector<edgeinfo> incident_nodes_;
	};


	std::vector<nodeinfo> nodes_;
	std::vector<size_type> i2u_;
};

#endif // CME212_GRAPH_HPP
