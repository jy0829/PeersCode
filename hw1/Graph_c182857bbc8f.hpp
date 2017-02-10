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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

  struct node_pair;
  struct node_element;

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
  /** Synonym for Node Value. */
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
      return (graph->nodes[idx]).point;  
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Get the value of the node. */
    node_value_type& value() {
      return (graph->nodes[idx]).node_value;
    }

    /** const version of value() */
    const node_value_type& value() const {
      return (graph->nodes[idx]).node_value;
    }

    /** Get the number of edges connected to the node. */
    size_type degree() const {
      int count = 0;
      for (int i = 0; i < graph->num_edge(); i++) {
        if (*this == graph->edges[i].node1() || *this == graph->edges[i].node2()){
          count++;
        }
      }
      return count;
    }

    /** Get the begin incident iterator of the node. */
    incident_iterator edge_begin() const {
      size_type i = 0;
      while (i < graph->num_edges() && *this != (*graph).edge(i).node1()
             && *this != (*graph).edge(i).node2()) {
        i++;
      }
      return IncidentIterator(graph, i, idx);
    }
  
    /** Get the end incident iterator of the node. */
    incident_iterator edge_end() const {
      size_type i = graph->num_edges();
      while (i > 0 && *this != (*graph).edge(i-1).node1() 
             && *this != (*graph).edge(i-1).node2()) {
        i--;
      }   
      return IncidentIterator(graph, i, idx);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const node_type& n) const {
      return (this->graph == n.graph && this->idx == n.idx);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const node_type& n) const {
      if (this->graph != n.graph) {
        return (this->graph < n.graph);
      } else {
	return (this->idx < n.idx);
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    //Pointer back to the graph.
    graph_type* graph;
    //Index of the node in the container
    size_type idx;

    //Construct valid node object 
    Node(const graph_type* g, size_type i)
        : graph(const_cast<Graph*>(g)), idx(i) {
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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& v = node_value_type()) {
    nodes.push_back(node_element(position, v));
    return Node(this, nodes.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const node_type& n) const {
    return (this == n.graph && nodes.size() > n.idx); 
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    assert(i < nodes.size());
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
    node_type node1() const {
      return (*(this->graph)).node(graph->edges[idx].node1);
    }

    /** Return the other node of this Edge */
    node_type node2() const {
      return (*(this->graph)).node(graph->edges[idx].node2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const edge_type& e) const {
      return ((this->graph == e.graph) 
              && (((this->graph->edges[idx].node1 == e.graph->edges[idx].node1)
              && (this->graph->edges[idx].node2 == e.graph->edges[idx].node2))
              || ((this->graph->edges[idx].node1 == e.graph->edges[idx].node2) 
              && (this->graph->edges[idx].node2 == e.graph->edges[idx].node2))));  
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const edge_type& e) const {
      if (this->graph != e.graph) {
        return (this->graph < e.graph);
      } else {
	return (this->idx < e.idx);
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    //Pointer back to graph
    graph_type* graph;
    //Index of the edge in the container
    size_type idx;

    //Construct valid edge
    Edge(const graph_type* g, size_type i)
        : graph(const_cast<Graph*>(g)), idx(i) {
    } 
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type edge(size_type i) const {
    assert(i < edges.size());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const node_type& a, const node_type& b) const {
    for (auto pair : edges) {
      if (((pair.node1 == a.idx) && (pair.node2 == b.idx))
          || ((pair.node1 == b.idx) && (pair.node2 == a.idx))) {
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
  edge_type add_edge(const node_type& a, const node_type& b) {
    for (size_type i = 0; i < edges.size(); i++) {
      if ((edges[i].node1 == a.idx && edges[i].node2 == b.idx)
          || (edges[i].node1 == b.idx && edges[i].node2 == a.idx)) {
        return Edge(this, i);
      }
    } 
    node_pair new_pair(a.idx, b.idx);
    edges.push_back(new_pair);
    return edge(edges.size() - 1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** @return the Node the iterator is pointing to. */ 
    Node operator*() const {
      return (*graph).node(idx); 
    }  
  
    /** Move the iterator to the next Node in the container. */
    NodeIterator& operator++() {
      (this->idx)++;      
      return *this;
    }
  
    /** Compare this Node Iterator to @a ni 
     *  Equal Iterators have the same graph pointer and the same index.
     */
    bool operator==(const NodeIterator& ni) const {
      return (this->graph == ni.graph && this->idx == ni.idx);
    } 


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    graph_type* graph;
    size_type idx; 
 
    NodeIterator(const graph_type* g, size_type i)
        : graph(const_cast<Graph*>(g)), idx(i) {
    }
   
  };

  // HW1 #2: YOUR/ CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Get the begin Node Iterator */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Get the end Node Iterator */
  node_iterator node_end() const {
    return NodeIterator(this, nodes.size());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    
    /** Get the Edge the incident iterator is pointing to. */
    Edge operator*() const {
      if ((*graph).edge(e_idx).node1() != (*graph).node(n_idx)) {
        graph->edges[e_idx].swap_order();  
      }
      return (*graph).edge(e_idx);      
    }

    /** Move the Incident Iterator to point the next incident edge. */  
    IncidentIterator& operator++() {
      if (IncidentIterator(graph, e_idx, n_idx) != (*graph).node(n_idx).edge_end()) {
        e_idx++;
        while (IncidentIterator(graph, e_idx, n_idx) != (*graph).node(n_idx).edge_end()
               && (*graph).edge(e_idx).node1() != (*graph).node(n_idx)
               && (*graph).edge(e_idx).node2() != (*graph).node(n_idx)) {
          e_idx++; 
        }
      }
      return *this;
    }

    /** Compare this Incident Iterator with @a iter  
     *  Equal Iterators have the same graph pointer, e_idx and n_idx
     */
    bool operator==(const IncidentIterator& iter) const {
      return (this->graph == iter.graph && this->e_idx == iter.e_idx
              && this->n_idx == iter.n_idx);
    }

   private:
    friend class Graph;
    
    graph_type* graph;
    size_type e_idx;
    size_type n_idx; 
 
    IncidentIterator(const graph_type* g, size_type e_i, size_type n_i )
        : graph(const_cast<Graph*>(g)), e_idx(e_i), n_idx(n_i) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:totally_ordered<EdgeIterator> {
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

    /** @return the Edge the iterator is pointing to.  */
    Edge operator*() const {
      return (*graph).edge(idx);
    }
 
    /** Move the Edge Iterator to point the next Edge in the container. */
    EdgeIterator& operator++() {
      (this->idx)++;
      return *this;
    }

    /** Compare this Edge iterator with @a ei. 
     *  Equal iterators have the same graph pointer and idx.
     */
    bool operator==(const EdgeIterator& ei) const {
      return (this->graph == ei.graph && this->idx == ei.idx);
    } 
   
 
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph;
    size_type idx; 
 
    EdgeIterator(const graph_type* g, size_type i)
        : graph(const_cast<Graph*>(g)), idx(i) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Get the begin edge iterator. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Get the end edge iterator. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }


 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  
  struct node_pair {
    size_type node1;
    size_type node2;

    node_pair(size_type node1, size_type node2) {
      this->node1 = node1;
      this->node2 = node2;
    }
  
    void swap_order() {
      auto temp = this->node1;
      this->node1 = node2;  
      this->node2 = temp;
    }
  };

  struct node_element {
    Point point;   
    node_value_type node_value;

    node_element(Point point, node_value_type node_value) {
      this->point = point;
      this->node_value = node_value;
    }
  };
  
  std::vector<node_element> nodes;
  std::vector<node_pair> edges;
};

#endif // CME212_GRAPH_HPP
