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

    
    class internal_node; // vector 
    class internal_edge; // 
    

 public:
    using node_value_type = V;

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

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
      num_edges_= 0;
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
      return fetch().point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return the value for the node */
    node_value_type& value() {
        return fetch().value_;
    }

    /** Return the value for the node */
    const node_value_type& value() const {
        return fetch().value_;
    }

    /** Return the number of edges for the given node */
    size_type degree() const{
        return fetch().neighbor_nodes.size();
    }

    /** Return an iterator that points to the first edge for the node */
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_, uid_, 0);
    }

    /** Return an iterator that points to the last edge for the node */
    incident_iterator edge_end() const {
        return IncidentIterator(graph_, uid_, degree() );
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        if (n.graph_ == graph_ && n.uid_ == uid_) {
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

        if (graph_ == n.graph_ && uid_ < n.uid_)
            return true;

        if (graph_ < n.graph_)
            return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
      graph_type* graph_;
      size_type uid_;

      /** Private Constructor */
      Node(const graph_type* graph, size_type uid) 
          : graph_(const_cast<Graph*>(graph)), uid_(uid) {
      }
      // Helper method to return the appropriate node
      internal_node& fetch() const {
        if (uid_ < graph_->size() ) {
            return graph_->node_vec[uid_];
        }
        //return graph_->node_vec[0];
        assert(false);
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_vec.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type() ) {
      internal_node node_new;
      node_new.value_ = val;
      node_new.point_ = position;
      node_vec.push_back(node_new);   
    return Node(this,node_vec.size()-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //internal_node temp_node;
    //temp_node.point_ = node_vec[n.uid_].point_;
    //temp_node.uid_ = n.index();
      
    //auto is_bool = std::find(node_vec.begin(), node_vec.end(), temp_node); 
    //if (is_bool != node_vec.end() ) {
    //   return true;
    //}
    // return false;
    return (n < size());
  }
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);        // Invalid node
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const { 
      return graph_->node(node2_uid);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ == e.graph_) {
        if (node1_uid == e.node1_uid && node2_uid == e.node2_uid) {
          return true;
        }
        if (node1_uid == e.node2_uid && node2_uid == e.node1_uid) {
          return true;
        }
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        if (graph_ == e.graph_) {
            return (node1_uid < e.node1_uid);
        }
        return (graph_ < e.graph_);
        
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
      
    Graph* graph_;
    size_type node1_uid;
    size_type node2_uid;
    
    /** Private Constructor */
    Edge(const Graph* graph , size_type uid1, size_type uid2) 
        : graph_(const_cast<Graph*>(graph)), node1_uid(uid1), node2_uid(uid2) {
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
    return Edge(this, edge_vec[i].node1, edge_vec[i].node2);    
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

      internal_node n = node_vec[a.index()];
  
      auto is_true = std::find(n.neighbor_nodes.begin(), n.neighbor_nodes.end(), b.index());     
      if (is_true != n.neighbor_nodes.end() ) {
          return true;
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
    size_type node1 = a.index();
    size_type node2 = b.index(); 

    if (has_edge(a,b) ){
        return Edge(this, node1, node2);
    }
    node_vec[node1].neighbor_nodes.push_back(node2);
    node_vec[node2].neighbor_nodes.push_back(node1);
    ++num_edges_;

    internal_edge e;
    e.node1 = node1;
    e.node2 = node2;
    edge_vec.push_back(e);
    return Edge(this, node1, node2);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    num_edges_ = 0;
    node_vec.clear();
    edge_vec.clear();

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

    /** Return the node indicated by the NodeIterator */
    Node operator*() const {
        return graph_->node(uid_);
    }

    /** Move to the next NodeIterator */
    NodeIterator& operator++(){
        if (uid_ < graph_->size()) {
            uid_++;
        }
        return *this;
    }

    /** Operator overloud that detemines if two nodes are equal. Returns bool. */
    bool operator==(const NodeIterator& n) const {
        return ( (n.graph_ == graph_ && n.uid_ == uid_) );
    }

    private:
    friend class Graph;
    graph_type* graph_;
    size_type uid_;

    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type uid)
    : graph_(const_cast<Graph*>(graph)), uid_(uid) {

    }
  };

  /** Return the node iterator that points to the first node of the graph */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return the node iterator that points to the last node of the graph */
  NodeIterator node_end() const {
    return NodeIterator(this, node_vec.size() );
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

    /** Return the edge indicated by the IncidentIterator */
    Edge operator*() const {
        return Edge(graph_, node_uid_, graph_->node_vec[node_uid_].neighbor_nodes[edge_uid_] );
    }

    /** Return the IncidentIterator that points to the next edge */
    IncidentIterator& operator++() {
        if (edge_uid_ < graph_->node_vec[node_uid_].neighbor_nodes.size()) {
            edge_uid_++;
        }
        return *this;
    }

    /** Operator overload that determines if two edges are equal. Returns bool. */
    bool operator==(const IncidentIterator& iit) const {
        return (graph_ == iit.graph_ && node_uid_ == iit.node_uid_ && iit.edge_uid_ == edge_uid_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_uid_;    // uid of node invoked
    size_type edge_uid_;    // uid of neighbor_nodes that contains edge nodes

    /** Private Constructor */
    IncidentIterator(const Graph* graph, size_type node, size_type edge)
    : graph_(const_cast<Graph*>(graph)), node_uid_(node), edge_uid_(edge) {
    }
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

    /** Return the edge indicated by the EdgeIterator */
    Edge operator*() const {
        return Edge(graph_, node1_, node2_);
    }

    /** Returns the EdgeIterator that points to the next edge */
    EdgeIterator& operator++() {
        if (edge_uid_ < graph_->edge_vec.size()){
            ++edge_uid_;
            node1_ = graph_->edge_vec[edge_uid_].node1;
            node2_ = graph_->edge_vec[edge_uid_].node2;
        }
        return *this;
    }


    /** Operator overload that determines if two edges are equal. Returns bool. */
    bool operator==(const EdgeIterator& eit) const {
        return (graph_ == eit.graph_ && node1_ ==  eit.node1_ && node2_ == eit.node2_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node1_;
    size_type node2_;
    size_type edge_uid_;

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type node1, size_type node2, size_type edge_uid)
    : graph_(const_cast<Graph*>(graph)), node1_(node1), node2_(node2), edge_uid_(edge_uid) {
    }

  };

  /** Return the edge iterator that points to the first edge of the graph */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edge_vec[0].node1, edge_vec[0].node2, 0);
  }

  /** Return the edge iterator that points to the last edge of the graph */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_vec[num_edges_-1].node1, edge_vec[num_edges_-1].node2, num_edges_-1);
  }

    private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.  
  class internal_edge {
     public:
     size_type node1;
     size_type node2;
     //size_type uid_;
     
      bool operator==(const internal_edge& e) const {
          if ( (node1 == e.node1 && node2 == e.node2) || (node2 == e.node1 && node1 == e.node2) ) {
              return true;
          }
          return false;
      } 
  };
    
  class internal_node {
      public:
      Point point_;
      node_value_type value_;
      std::vector<size_type> neighbor_nodes;

  };
    
  std::vector<internal_node> node_vec;
  std::vector<internal_edge> edge_vec;
  size_type num_edges_;
              
};

#endif // CME212_GRAPH_HPP
