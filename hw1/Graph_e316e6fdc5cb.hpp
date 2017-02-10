#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
To run:
./viewer tiny.nodes tiny.tets
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  using node_value_type = V;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes_(), nodeVals_(), edgeList_(), adjList_() {
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
      return graph_->nodes_[idx_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
    }

    
    /** Return this node's value by reference, of type @a node_value_type
    * Serves as a getter of the node's value
    * Can also serve as a setter for the node's value
    */
    node_value_type& value() {
	return graph_->nodeVals_[index()];
    }

    /** Return this node's value as a const reference, of type @a node_value_type
    * This function is a strict getter of the node's value
    */
    const node_value_type& value() const {
	return (const node_value_type&) graph_->nodeVals_[index()];
    }

    // HW1: YOUR CODE HERE
    /** Returns the degree of the current node
     * also known as the number of edges incident to this node
     */
    size_type degree() const {
	return graph_->adjList_[index()].size();
    }
    /** Returns a new incident edge iterator with:
     * 1) the current graph pointer
     * 2) the index of the node in consideration (the one who's 
     *         edges we are iterating through)
     * 3) the beginning of the iterator through that node's adjacency set
     *        (this set contains all of the nodes it is attached to)
     */
    incident_iterator edge_begin() const {
	return IncidentIterator(graph_, index(), graph_->adjList_[index()].begin());
    }
    /** Returns a new incident edge iterator with:
     * 1) the current graph pointer
     * 2) the index of the node in consideration (the one who's 
     *         edges we are iterating through)
     * 3) the end of the iterator of that node's adjacency set
     *        (this set contains all of the nodes it is attached to)
     */
    incident_iterator edge_end() const {
	return IncidentIterator(graph_, index(), graph_->adjList_[index()].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
	return ((graph_ == n.graph_) && (idx_ == n.index()));
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
	if (graph_ < n.graph_) {
		return true;
	} else if ((graph_ == n.graph_) && (idx_ < n.idx_)) {
		return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
	Graph* graph_;
	size_type idx_;

	Node(const Graph* graph, size_type idx) : graph_(const_cast<Graph*>(graph)), idx_(idx) {}
  };

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
    nodes_.push_back(position);
    nodeVals_.push_back(val);

    // add an empty vector into the adjacency list
    // so it can later store edges
    std::set<size_type> newSet;
    adjList_.push_back(newSet);
    
    Node n = Node(this, (size()-1));
    return n;    
  }



  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_ == this);

    //return(n.index() < size());
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

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, n1ID_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, n2ID_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_) {
	// the edges are not equal because they are not from the same graph
	return false;
      }
      if ((this->node1() == e.node1()) && (this->node2() == e.node2())) {
	return true;
      } else if ((this->node1() == e.node2()) && (this->node2() == e.node1())) {
	return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
	if (this->graph_ < e.graph_) {
	// the edges are not equal because they are not from the same graph
	return true;
      } else if ((graph_ == e.graph_)) {
	Node thisMin = std::min(this->node1(), this->node2());
	Node eMin = std::min(e.node1(), e.node2());
	Node thisMax = std::max(this->node1(), this->node2());
	Node eMax = std::max(e.node1(), e.node2());

      	if (thisMin < eMin) {
		return true;
      	} else if ((thisMin == eMin) && (thisMax < eMax)) {
		return true;
      	}
      }
      return false;

    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects


    Graph* graph_;
    size_type n1ID_;
    size_type n2ID_;

    Edge(const Graph* graph, size_type n1ID, size_type n2ID) : graph_(const_cast<Graph*>(graph)), n1ID_(n1ID), n2ID_(n2ID) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edgeList_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, edgeList_[i].first, edgeList_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    std::set<size_type>::iterator it = adjList_[a.index()].find(b.index());
    return (it != adjList_[a.index()].end());
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
    if (has_edge(a, b) == true) {
	// return that edge
 	return Edge(this, a.index(), b.index());
    }
    // otherwise, create a new one
    // put it into the edgeList_
    std::pair<size_type, size_type> thisPair (a.index(), b.index());
    edgeList_.push_back(thisPair);

    // and also put it into the adjList_ for both nodes
    adjList_[a.index()].insert(b.index());
    adjList_[b.index()].insert(a.index());
    

    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    nodeVals_.clear();
    edgeList_.clear();
    adjList_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
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
    /** Dereferences the node iterator's current position and returns
     * the corresponding node object
     */
    Node operator*() const {
	return Node(graph_, nodeID_);
    }
    /** Increments the position of the node iterator
     */
    NodeIterator& operator++() {
	++nodeID_;
	return *this;
    }
    /** Test whether this node iterator and @a ni are equal.
     *
     * Equal node iterators represent the same graph pointers and the same nodeID's.
     */
    bool operator==(const NodeIterator& ni) const {
	return ((graph_ == ni.graph_) && (nodeID_ == ni.nodeID_));
    }

    /** Test whether this node iterator and @a ni are not equal.
     *
     * Unequal node iterators represent either:
     * 1) different graph pointers
     * 2) different nodeID's.
     */
    bool operator!=(const NodeIterator& ni) const {
	return !(*this == ni);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
	Graph* graph_;
	size_type nodeID_;

	// construct a valid node iterator
	NodeIterator(const Graph* graph, size_type nodeID) : graph_(const_cast<Graph*>(graph)), nodeID_(nodeID) {}
  };

  // HW1 #2: YOUR CODE HERE
  /** Returns a new node iterator with the current graph pointer
   * as well as the beginning point for the iterator, i.e. at node of index 0.
   * This node of index 0 is the first node ever added to the Graph.
   */
  node_iterator node_begin() const {
	return node_iterator(this, 0);
  }
  /** Returns a new node iterator with the current graph pointer
   * as well as the end, at an index after the last node that exists
   * such that the node_[index] does not exist
   */
  node_iterator node_end() const {
	return node_iterator(this, this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
    /** Dereferences the incident iterator's current position and returns
     * the corresponding edge object
     */
    Edge operator*() const {
	return (Edge(graph_, node1ID_,*setIT_));
    }

    /** Increments the position of the incident iterator
     */
    IncidentIterator& operator++() {
	++setIT_;
	return *this;
    }

    /** Test whether this incident iterator and @a iit are not equal.
     *
     * Equal incident iterators have the same:
     * 1) graph pointers
     * 2) node1ID's.
     * 3) pointer to a set iterator.
     */
    bool operator==(const IncidentIterator& iit) const {
	return ((graph_ == iit.graph_) && (node1ID_ == iit.node1ID_) && (setIT_ == iit.setIT_));
    }

    /** Test whether this incident iterator and @a iit are not equal.
     *
     * Unequal incident iterators represent either:
     * 1) different graph pointers
     * 2) different node1ID's.
     * 3) different set iterator pointers.
     */
    bool operator!=(const IncidentIterator& iit) const {
	return !(*this == iit);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

	Graph* graph_;
	size_type node1ID_;
	std::set<size_type>::iterator setIT_;

	// construct a valid incident iterator
	IncidentIterator(const Graph* graph, size_type node1ID, std::set<size_type>::iterator setIT) : graph_(const_cast<Graph*>(graph)), node1ID_(node1ID), setIT_(setIT) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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
    /** Dereferences the edge iterator's current position and returns
     * the corresponding edge object
     */
    Edge operator*() const {
	return Edge(graph_, graph_->edgeList_[edgeID_].first, graph_->edgeList_[edgeID_].second);
    }

    /** Increments the position of the edge iterator
     */
    EdgeIterator& operator++() {
	++edgeID_;
	return *this;
    }

    /** Test whether this edge iterator and @a EI are equal.
     *
     * Equal edge iterators all have the same:
     * 1) graph pointers.
     * 2) edgeID's.
     */
    bool operator==(const EdgeIterator& EI) const {
	return ((graph_ == EI.graph_) && (edgeID_ == EI.edgeID_));
    }

    /** Test whether this edge iterator and @a EI are not equal.
     *
     * Unequal edge iterators represent either:
     * 1) different graph pointers
     * 2) different edgeID's
     */
    bool operator!=(const EdgeIterator& EI) const {
	return !(*this == EI);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

	Graph* graph_;
	size_type edgeID_;

	// construct a valid edge iterator
	EdgeIterator(const Graph* graph, size_type edgeID) : graph_(const_cast<Graph*>(graph)), edgeID_(edgeID) {}
	
  };

  // HW1 #5: YOUR CODE HERE
  /** Returns a new edge iterator with the current graph pointer
   * as well as the beginning point for the iterator, i.e. an edge of index 0.
   * This edge of index 0 is the first edge ever added to the Graph.
   */
  edge_iterator edge_begin() const {
	return EdgeIterator(this, 0);
  }
  /** Returns a new edge iterator with the current graph pointer
   * as well as the end, at an index after the last edge that exists
   * such that the edgeList_[index] does not exist
   */
  edge_iterator edge_end() const {
	return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  
  std::vector<Point> nodes_;
  std::vector<node_value_type> nodeVals_;

  std::vector<std::pair<size_type, size_type>> edgeList_;
  std::vector<std::set<size_type>> adjList_;


};

#endif // CME212_GRAPH_HPP
