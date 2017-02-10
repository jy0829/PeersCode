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

  /** Type of node value*/
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->node_container[index_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /** Return this node's value, the type of value matches node_value_type*/
    node_value_type& value() {
    	return graph_->node_container[index_].value;
    }
    
    /** Return this node's value as const, the type of value matches node_value_type*/
    const node_value_type& value() const {
    	return graph_->node_container[index_].value;
    }

    /** Return this node's value as const, the type of value matches node_value_type*/
    size_type degree() const {
    	return graph_->node_connection[index_].size();
    }

    incident_iterator edge_begin() const {
    	return IncidentIterator(graph_, index_,0);
    }
    incident_iterator edge_end() const {
    	return IncidentIterator(graph_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((n.index() == this->index_) && (n.graph_ == this->graph_))
         return true;
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
      // HW0: YOUR CODE HERE
      //In the same graph
      if (this->graph_ == n.graph_) {
      	if (n.index() < index_) {        	
		return true;
	} else {		
		return false;
        }
      }
      //In the different graphs.
      if (this->graph_ < n.graph_)
      	return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    /** Pointer to Graph Container */
    Graph* graph_;

    /** Node Index */
    size_type index_;
    

    /** Private Constructor */
    Node(const Graph* graph)
        : graph_(const_cast<Graph*>(graph)){}

    Node(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index){} 
 
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_container.size();
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
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    internal_node newnode;
    newnode.position = position;
    newnode.value = node_value;
    node_container.push_back(newnode);       
    return Node(this,node_container.size()-1);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this->graph_ == n.graph_)
	return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);        
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_index); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_index);    
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ == e.graph_) {
	if ((this->node1() == e.node1()) && (this->node2() == e.node2()))
         	return true;
      	if ((this->node1() == e.node2()) && (this->node2() == e.node1()))
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
      if (this->graph_ == e.graph_) {
      	if (this->node1_index < e.node1_index) {
		return true;
	}
        if ((this->node1_index == e.node1_index) && (this->node2_index < e.node2_index))
		return true;
      }	
      
      if (this->graph_ < e.graph_)
      	return true;

      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    /** Pointer to Graph Container */
    Graph* graph_;

    /** Edge Index */
    size_type index_;
   

    /** Node Index */
    size_type node1_index;
    size_type node2_index;

    /** Private Constructor */ 
    Edge(const Graph* graph, size_type index, size_type n1, size_type n2)
        : graph_(const_cast<Graph*>(graph)), index_(index), node1_index(n1), node2_index(n2){} 
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_container.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this,i,edge_container[i].node1_index, edge_container[i].node2_index);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (edge_container.size() != 0) {
        for (size_type i=0;i<edge_container.size();i++ ) {
            if ((edge_container[i].node1_index == a.index()) && (edge_container[i].node2_index == b.index()))
                return true;
            if ((edge_container[i].node1_index == b.index()) && (edge_container[i].node2_index == a.index()))
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
    internal_edge new_edge;
    new_edge.node1_index = a.index();
    new_edge.node2_index = b.index();
    
    //If edge is not in the graph.
    if (!has_edge(a,b)) {
    	edge_container.push_back(new_edge);
        node_connection[a.index()].push_back(b.index());
        node_connection[b.index()].push_back(a.index());
    	return Edge(this,edge_container.size()-1,a.index(),b.index());
    }

    //If edge exists in the graph.
    for (size_type i=0;i<edge_container.size();i++ ) {
            if ((edge_container[i].node1_index == a.index()) && (edge_container[i].node2_index == b.index()))
                return Edge(this,i,a.index(), b.index());
            if ((edge_container[i].node1_index == b.index()) && (edge_container[i].node2_index == a.index()))
                return Edge(this,i,b.index(),a.index());
    }        
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_container.clear();
    edge_container.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator(){
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const{
    	return Node(graph_,n);
    }
    
    NodeIterator& operator++() {
	n++;
	return *this;
    }
    bool operator==(const NodeIterator& Iter2) const {
    	if ((this->graph_ == Iter2.graph_) && (this->n == Iter2.n))
		return true;
	return false;
    }

	
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
     
    /**Pointer to Node*/
    size_type n;

    /** Pointer to Graph Container */
    Graph* graph_;
    
    /** Private Constructor */  
    NodeIterator(const Graph* graph, size_type new_n):graph_(const_cast<Graph*>(graph)), n(new_n) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const {
  	return node_iterator(this,0);
  }
  node_iterator node_end() const {
  	return node_iterator(this, num_nodes());
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const {
    	return Edge(this->graph_, node1_index, this->graph_.node_connection[node1_index][node2_connection]);
    }

    IncidentIterator& operator++() {
 
	//Last of node1's connections. 	
 	if (node2_connection == this->graph_.node_connection[node1_index].size())
		return *this;
        
	//Iterate next node node1 connected to.
        node2_connection++;
	return *this;
    }

    bool operator==(const IncidentIterator& Iter2) const {
    	if ((this->graph_ == Iter2.graph_) && (this->node1_index == Iter2.node1_index) && (this->node2_connection == Iter2.node2_connection))
		return true;
    	return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type node1_index;
    size_type node2_connection;
    Graph* graph_;

    /** Private Constructor */  
    IncidentIterator(const Graph* graph, size_type n1, size_type n2):graph_(const_cast<Graph*>(graph)), node1_index(n1), node2_connection(n2) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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
    //bool operator==(const EdgeIterator& Iter2) const
    
    Edge operator*() const {
		return Edge(this->graph_, node1_index, this->graph_.node_connection[node1_index][node2_connection]);
    	
    }
    
    EdgeIterator& operator++() {
	while (node1_index > graph_->node_connection[node1_index][node2_connection]) {
		if (node2_connection < graph_->node_connection[node1_index].size()-1) {
			node2_connection++;
		} else {
			if (node1_index == graph_->node_connection.size()-1) {
					node2_connection = graph_->node_connection[node1_index].size();
					break;		
			}
			node1_index++;
			node2_connection = 0;
		}		        
	}    	
	return *this;
    }

    bool operator==(const EdgeIterator& Iter2) const {
    	if ((this->graph_ == Iter2.graph_) && (this->node1_index == Iter2.node1_index) && (this->node2_connection == Iter2.node2_connection))
		return true;
    	return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type node1_index;
    size_type node2_connection;

    /** Private Constructor */  
    EdgeIterator(const Graph* graph, size_type n1, size_type n2):graph_(const_cast<Graph*>(graph)), node1_index(n1), node2_connection(n2) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  /** Internal struct for Nodes */
  struct internal_node {
     Point position;
     node_value_type value;
     size_type degree;
     internal_node(): value(node_value_type()), degree(size_type()) {}
  };
  /** Internal struct for Edges */
  struct internal_edge {
     size_type node1_index;
     size_type node2_index;
  };

  /** Vector to store nodes that connected */
  std::vector<std::vector<size_type>> node_connection;

  /** Vector to store Nodes and Edges */
  std::vector<internal_node> node_container;
  std::vector<internal_edge> edge_container;

  


};


#endif // CME212_GRAPH_HPP
