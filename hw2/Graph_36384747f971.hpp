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
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

 private:

  //node internal struct
  struct node_info;

  //edge interal struct
  struct edge_info;

  //track total number of edges
  size_type total_edges = 0;

  //node array for graph
  std::vector<node_info> nodes_vec;

  //adjaceny array for graph
  std::vector<std::vector<edge_info>> adj_vec;

  //node indices vec
  std::vector<size_type> idx2nid;

 public:
 
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  using node_value_type = V;

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


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes_vec(), adj_vec() {
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
      //no additiontional code needed to construct an invalid node
    }

    /** Return this node's position. */
    Point& position() { 
     //get node struct from graph and extract position   
      return n_graph_->nodes_vec[node_id_].P;
    }

    /** Return this node's position. */
    const Point& position() const { 
     //get node struct from graph and extract position   
      return n_graph_->nodes_vec[node_id_].P;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //check indices vector
      for (size_type i = 0; i < n_graph_->idx2nid.size(); ++i) {
        if (node_id_ == n_graph_->idx2nid[i]) {
          return i;
        }
      }
      return 0;
    }

    /**
    *value() can be used to add additional information about the node
    *@return @a val by reference for Node of type V
    * 
    *@pre the function must be called on a valid node
    */
    node_value_type& value() {
      return n_graph_->nodes_vec[node_id_].val;
    }

    /**
    *value() const can be used to add additional information about the node
    *@return @a val of type const by reference for node of type V
    *
    *@pre The graph pointer must point to a valid graph object.
    */

    const node_value_type& value() const {
      return n_graph_->nodes_vec[node_id_].val;
    }

    /**
    *degree() determines the total number of neighbors of a given node
    *@return @a d of type size_type
    *
    *@pre the function must be called on a valid node
    */

    size_type degree() const {
      size_type d = n_graph_->adj_vec[node_id_].size();
      return d;
    }

    /**
    *edge_begin() constructs @a EB, an Incident_Iterator that points to the first
    *neighbor of a root node 
    *@return @a EB of type incident_iterator
    *
    *@pre the function must be called on a valid node
    */

    incident_iterator edge_begin() const {
      //use constructor
      incident_iterator I(node_id_, 0, const_cast<Graph*>(n_graph_));
      return I;   
    }

    /**
    *edge_end() constructs @a EE, an Incident_Iterator that points one past the 
    *last neighbor of a root node 
    *@return @a EE of type incident_iterator
    *
    *@pre the function must be called on a valid node
    */

    incident_iterator edge_end() const {
      //use constructor
      incident_iterator I(node_id_, degree(), const_cast<Graph*>(n_graph_));
      return I; 
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // check if they have the same graph and index
      if (n.node_id_ == node_id_ && n.n_graph_ == n_graph_) {
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
      //if they are in the same graph, order based on node id
      if ( n.n_graph_ == n_graph_) {
        if (node_id_ < n.node_id_) {
          return true;
        }
      }
      /*if different graphs, global order will depend on graph location 
      in memory, i.e. if mem location of graph of node n is larger than the
      mem location of this graph, node n is larger in global order, regardless 
      of node id */
      else if (n_graph_ <  n.n_graph_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    //pointer back to the Graph (8 bits)
    Graph* n_graph_;

    //node id (8 bits)
    size_type node_id_;

  }; //end Node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    //get length of vector holding unique node info
    size_type total_nodes = idx2nid.size();
    return total_nodes;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] a A value associated with the node
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& a = node_value_type()) {
    
    //create and populate struct for node being added
    node_info current_node;
    current_node.P = position;
    current_node.val = a;
   
    //add this node struct to the graph
    nodes_vec.push_back(current_node);

    //add to indices vector
    idx2nid.push_back(nodes_vec.size()-1);

    //add placeholder to adj vec to partition space for this node
    adj_vec.emplace_back(std::vector<edge_info>());

    //create node to be returned and populate
    Node N;
    N.n_graph_ = const_cast<Graph*>(this);
    //get node id
    N.node_id_ = nodes_vec.size()-1;
    
    return N;     
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //make sure graph is the same, evaluate if index is larger than 
    //indices vector size
    if (n.index() < idx2nid.size() && n.n_graph_ == const_cast<Graph*>(this)) {
      return true;
    }
    return false;
  }

  /** Remove a node from the graph and return an unsigned integer. 
   * @param[in] n the node that is to be removed from the graph.
   * @return Returns an unsigned integer (in this case, zero).
   * @pre The node @a n must be valid and have a unique id.
   * @post The graph has one less node, and there is one less entry in the 
   * indices vector. 
   *
   *Complexity: On the order of the degree of the node.
   */
  size_type remove_node(const Node& n) {
    idx2nid.erase(idx2nid.begin()+n.index());
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      size_type n_1 = (*i).n_1_id_;
      size_type n_2 = (*i).n_2_id_;
      for (size_type j = 0; j < (*i).node2().degree(); ++j) {
        if (adj_vec[n_2].size() > j) {
          if (adj_vec[n_2][j].N_2_id == n_1) {
            adj_vec[n_2].erase(adj_vec[n_2].begin()+j);
            total_edges = total_edges-1;
            break;
          }
        }
      }
    }
    return 0;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    //node to be returned
    Node node_at_i;
    //assign node id
    node_at_i.node_id_ = idx2nid[i];
    //assign graph
    node_at_i.n_graph_ = const_cast<Graph*>(this);

    return node_at_i;   
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
      //no additional code needed to construct an invalid edge
    }

    /** Return a node of this Edge */
    Node node1() const {
      //construct node 1
      Node n1;
      n1.node_id_ = n_1_id_;
      n1.n_graph_ = const_cast<Graph*>(e_graph_);
      return n1; 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //construct node 2
      Node n2;
      n2.node_id_ = n_2_id_;
      n2.n_graph_ = const_cast<Graph*>(e_graph_);
      return n2;   
    }

    //edge length calculation
    double length() const {
      return norm(node1().position()-node2().position());
    }

/**
    *evalue() can be used to add additional information about the edge
    *@return @a val by reference for Edge of type E
    * 
    *@pre the function must be called on a valid node
    */
    edge_value_type& evalue() {
      for (size_type i = 0; i < e_graph_->adj_vec[node1().node_id_].size(); i++) {
        if (node2().node_id_ == (e_graph_->adj_vec[node1().node_id_][i]).N_2_id) {
          return e_graph_->adj_vec[node1().node_id_][i].e_val;
        }
        
      }
      return e_graph_->adj_vec[node1().index()][0].e_val;       
        
    }

    /**
    *evalue() const can be used to add additional information about the edge
    *@return @a val of type const by reference for edge of type E
    *
    *@pre The graph pointer must point to a valid graph object.
    */

    const edge_value_type& evalue() const {
      for (size_type i = 0; i < e_graph_->adj_vec[node1().node_id_].size(); i++) {
        if (node2().node_id_ == (e_graph_->adj_vec[node1().node_id_][i]).N_2_id) {
          return e_graph_->adj_vec[node1().index()][i].e_val;
        }

      }
        return e_graph_->adj_vec[node1().index()][0].e_val;
      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //check if (a,b)=(a,b) or (a,b)=(b,a) using min and max node indices
      size_type min_node = std::min(e.n_1_id_, e.n_2_id_);
      size_type max_node = std::max(e.n_1_id_, e.n_2_id_);
      if (std::min(n_1_id_, n_2_id_) == min_node && std::max(n_1_id_, n_2_id_) == max_node && e.e_graph_ == e_graph_) {
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
      //for edges in the same graph, will be determined by smallest node index. 
      //If they have the same min node index, will then be decided by smallest
      //second node.
      if (e_graph_ == e.e_graph_) {
        if (std::min(n_1_id_, n_2_id_) < std::min(e.n_1_id_, e.n_2_id_)) {
          return true;
        }
        else if (std::min(n_1_id_, n_2_id_) == std::min(e.n_1_id_, e.n_2_id_)) {
          if (std::max(n_1_id_, n_2_id_) < std::max(e.n_1_id_, e.n_2_id_)) {
            return true;
          }
        }
      }
      /* if they aren't in the same graph, use smallest graph memory location
      (similar to node global order */
      else if (e_graph_ < e.e_graph_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    //pointer to graph (8 bits)
    Graph* e_graph_;

    //node 1 id (8 bits)
    size_type n_1_id_;

    //node 2 id (8 bits)
    size_type n_2_id_;

  }; //end Edge class

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return total_edges;
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
    //see if node pair exists in adjacency matrix
    for (size_type i = 0; i < adj_vec[a.node_id_].size(); i++) {
      if (b.node_id_ == (adj_vec[a.node_id_][i]).N_2_id) {
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
    //check if edge currently exists
    if (has_edge(a, b)) {
      //create edge from these nodes
      Edge e_return;
      e_return.n_1_id_ = a.node_id_;
      e_return.n_2_id_ = b.node_id_;

      e_return.e_graph_ = const_cast<Graph*>(this);
      return e_return;
    }
    //if not, build edges struct and add it to the graph
    total_edges = total_edges + 1;
    edge_info a_edge;
    a_edge.N_2_id = b.node_id_;

    edge_info b_edge;
    b_edge.N_2_id = a.node_id_;

    //add node-edge info to the adjacency vector
    adj_vec[a.node_id_].push_back(a_edge);
    adj_vec[b.node_id_].push_back(b_edge);

    //create edge from these nodes
    Edge e_new;
    e_new.n_1_id_ = a.node_id_;
    e_new.n_2_id_ = b.node_id_;
    e_new.e_graph_ = const_cast<Graph*>(this);

    return e_new;        
  }

  /** Remove an edge from the graph and return an unsigned integer.
  * 
  * @param[in] a,b Nodes of the edge to be removed from the graph.
  * @return Returns an unsigned integer (zero, in this case).
  * @pre The nodes @a a and @a b must be valid.
  * @post The graph will have one less edge. The degree of nodes @a a and @a b
  * have been decreased by one.
  *
  *Complexity is of the order of the largest degree of the two nodes.
  */
  size_type remove_edge(const Node& a, const Node& b) {
    //remove for node a = node 1
    for (size_type i = 0; i < a.degree(); ++i) {
      if (adj_vec[a.node_id_].size() > i) {
        if (adj_vec[a.node_id_][i].N_2_id == b.node_id_) {
          adj_vec[a.node_id_].erase(adj_vec[a.node_id_].begin()+i);
          //update total number of edges
          total_edges = total_edges-1;
        }
      }
    }
    //remove for node b = node 1
    for (size_type i = 0; i < b.degree(); ++i) {
      if (adj_vec[b.node_id_].size() > i) {
        if (adj_vec[b.node_id_][i].N_2_id == a.node_id_) {
          adj_vec[b.node_id_].erase(adj_vec[b.node_id_].begin()+i);
        }
      }
    }

    return 0;
  }

  /** Remove an edge from the graph and return an unsigned integer.
  * 
  * @param[in] e Edge to be removed from the graph.
  * @return Returns an unsigned integer (zero, in this case).
  * @pre The edge @a e must be valid.
  * @post The graph will have one less edge. The degree of nodes @a a and @a b
  * have been decreased by one.
  *
  *Complexity is of the order of the largest degree of the two nodes (comes 
  *from calling remove_edge(const Node& a, const Node& b)).
  */

  size_type remove_edge(const Edge& e) {
  //edge nodes
  Node a = e.node1();
  Node b = e.node2();
  //utilize other remove edge function
  return remove_edge(a, b);
  }

  /** Remove an edge from the graph viad edge iterator and return the edge 
  * iterator of the next edge.
  * 
  * @param[in] e_it Edge iterator to be used to remove an edge from the graph.
  * @return Returns the edge iterator immediately after @a e_it.
  * @pre The edge iterator @a e_it must be valid, i.e. it cannot be at the last
  *element in the list of edges.
  * @post The graph will have one less edge. The degree of nodes @a a and @a b
  * have been decreased by one.
  *
  *Complexity is of the order of the largest degree of the two nodes (comes 
  *from calling remove_edge(const Edge& e)).
  */

  edge_iterator remove_edge (edge_iterator e_it) {
    Edge current = *e_it;
    Edge next = ++e_it;
    remove_edge(current);
    return next;
  } 

 
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // clear the vectors
    nodes_vec.clear();
    adj_vec.clear();
    idx2nid.clear();
    total_edges = 0;
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

    /**
    *operator*() returns the Node that the NodeIterater is currently pointing to.
    *@return returns the Node @a N
    *
    *@pre iterator != graph.node_end() : The iterator cannot be at the end of 
    *  the nodes 
    */
    Node operator*() const {

      return ni_graph_->node(p_) ;
    }

    /**
    *operator++() increments the iterator by one position
    *@return returns a reference to the incremented NodeIterator.
    *
    *@pre iterator != graph.node_end() : The iterator cannot be at the end of 
    *  the nodes  
    */
    NodeIterator& operator++() {
     
      this->p_ = p_+1;
      return *this;
    }

    /**
    *operator==() compares the current NodeIterator to a provided NodeIterator
    *and determines whether they are equal 
    *
    *@param[in] @a NI the input node iterator for comparison
    *@return returns @a b, a boolean value 
    */
    bool operator==(const NodeIterator& n) const {
      if (n.ni_graph_ == ni_graph_  && n.p_ == p_) {
        return true;
      }
      return false;
    }

    /**
    *operator!=() compares the current NodeIterator to a provided NodeIterator
    *and determines if they are not equal (uses operator==() for evaluation)
    *
    *@param[in] @a NI the input node iterator for comparison
    *@return returns @a b, a boolean value 
    */

    bool operator!=(const NodeIterator& n) const {
      return !(*this == n);
    }

  /** Remove a node from the graph using a node iterator and return 
   * a node iterator for the next element in the array. 
   * @param[in] n_it the node iterator that will be used to remove a node from
   * the graph.
   * @return Returns a node iterator (the iterator for the next node in sequence).
   * @pre The node iterator @a n_it must be valid, i.e. it cannot be pointing
   * to the end of the array.
   * @post The graph has one less node, and there is one less entry in the 
   * indices vector. 
   *
   *Complexity: On the order of the degree of the node, since it calls remove_node(Node& n).
   */

    NodeIterator remove_node(node_iterator n_it) {
      Node n = *n_it;
      NodeIterator next_n = ++n_it;
      ni_graph_->remove_node(n);
      return next_n; 
    }

   private:

    friend class Graph;

    size_type p_;
    Graph* ni_graph_;
  
    NodeIterator(size_type p1, const Graph* g) : p_(p1), ni_graph_(const_cast<Graph *>(g)) {
    }
  };

  /**
  *node_begin() returns the pointer to the beginning of the stored node structure
   *@return returns @NB, the node_iterator that points to the node id with zero
  */
  node_iterator node_begin() const {

    return NodeIterator(0,this);
  }

  /**
  *node_end() returns the pointer to one past the end of the stored node structure
  *@returns returns @NE, the node_iterator that points to one position after
  *  the last stored node
  */
  node_iterator node_end() const {
  
    return NodeIterator(idx2nid.size(),this);
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

    /**
    *operator*() returns the Edge that the IncidentIterater is currently pointing to.
    *@return returns the Edge @a E
    *
    *@pre IncidentIterator != node.edge_end() : The iterator cannot be at the end of 
    *  the edges of the node
    */

    Edge operator*() const {
      Edge Ed; //invalid edge
      //populate the edge
      Ed.n_1_id_ = n1i_;
      Ed.n_2_id_ = ii_graph_->adj_vec[n1i_][ii_].N_2_id;
      Ed.e_graph_ = ii_graph_;

      return Ed;
    }

    /**
    *operator++() increments the iterator by one position
    *@return returns a reference to the incremented IncidentIterator.
    *
    *@pre iterator != node.edge_end() : The iterator cannot be at the end of 
    *  the edges of the node  
    */

    IncidentIterator& operator++() {
      this->ii_ = ii_+1;
      return *this;
    }

    /**
    *operator==() compares the current IncidentIterator to a provided IncidentIterator
    *and determines whether they are equal
    *
    *@param[in] @a II the input incident iterator for comparison
    *@return returns @a b, a boolean value 
    */

    bool operator==(const IncidentIterator& i) const {
      if (i.ii_graph_ == ii_graph_  && i.ii_ == ii_ && i.n1i_ == n1i_) {
        return true;
      }
      return false;
    }

    /**
    *operator!=() compares the current IncidentIterator to a provided IncidentIterator
    *and determines if they are not equal (uses operator==() for evaluation)
    *
    *@param[in] @a II the input incident iterator for comparison
    *@return returns @a b, a boolean value 
    */

    bool operator!=(const IncidentIterator& i) const {
      return !(*this == i);
    }

   private:
    friend class Graph;

   //outermost iterator (based on node id)
   size_type n1i_;
   //inner iterator (iterate across connectivity for a specific node)
   size_type ii_;
   //graph pointer
   Graph* ii_graph_;

    IncidentIterator(size_type n1, size_type i1, const Graph* g) : n1i_(n1), 
    ii_(i1), ii_graph_(const_cast<Graph *>(g)) {
    }
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

    /**
    *operator*() returns the Edge that the EdgeIterater is currently pointing to.
    *@return returns the Edge @a E
    *
    *@pre EdgeIterator != graph.edge_end() : The iterator cannot be at the end of 
    *  the edges of the graph
    */

    Edge operator*() const {
      Edge Ed;
      Ed.e_graph_ = ei_graph_;
      Ed.n_1_id_ = (*ni_).index();
      Ed.n_2_id_ = (*iii_).node2().index();

      return Ed;
    }

    
    /**
    *operator++() increments the iterator by one position
    *@return returns a reference to the incremented EdgeIterator.
    *
    *@pre iterator != graph.edge_end() : The iterator cannot be at the end of 
    *  the edges of the graph  
    */

    EdgeIterator& operator++() {
      ++iii_;
      FML_II();
      return *this;
   }

    /**
    *operator==() compares the current EdgeIterator to a provided EdgeIterator
    *and determines whether they are equal
    *
    *@param[in] @a EI the input edge iterator for comparison
    *@return returns @a b, a boolean value 
    */

    bool operator==(const EdgeIterator& e) const {
      return (e.ei_graph_ == e.ei_graph_ && e.ni_ == ni_ && e.iii_ == iii_);
    }

    /**
    *operator!=() compares the current EdgeIterator to a provided EdgeIterator
    *and determines if they are not equal 
    *
    *@param[in] @a EI the input edge iterator for comparison
    *@return returns @a b, a boolean value 
    */

    /**
    *operator!=() compares the current EdgeIterator to a provided EdgeIterator
    *and determines if they are not equal 
    *
    *@param[in] @a EI the input edge iterator for comparison
    *@return returns @a b, a boolean value 
    */

    bool operator!=(const EdgeIterator& e) const {
      return !(e.ei_graph_ && e.ei_graph_ && e.ni_ == ni_ && e.iii_ == iii_);
    }
    
   private:
    friend class Graph;

    NodeIterator ni_ = (*ei_graph_).node_begin();
    IncidentIterator iii_ = (*ni_).edge_begin();
    Graph* ei_graph_;

    void FML_II() {
      while (ni_ != (*ei_graph_).node_end()) { //not at last node
        while (iii_ != (*ni_).edge_end()) { //not at end of incidents
            //guarantee specific node ordering to prevent duplicate edges
            if ((*ni_) < (*iii_).node2()) {
            ++iii_;
          }
          else {
            break;
          }
        } //end inner while
        if (iii_ == (*ni_).edge_end()) { //update new positions
          ++ni_;
          iii_ = (*ni_).edge_begin();
        }
        else {
          break;
        }
      } //end outer while

      if (ni_ == (*ei_graph_).node_end() ) {
        ni_ = (*ei_graph_).node_begin();
        iii_ = (*ni_).edge_begin();
      }
    } //end function


    EdgeIterator(NodeIterator n1, IncidentIterator i1, const Graph* g) : 
    ni_(n1), iii_(i1), ei_graph_(const_cast<Graph *>(g)) {
      FML_II();
    }

  };

  /**
  *edge_begin() returns the pointer to the first edge of a graph
   *@return returns @EB, an edge_iterator 
  */

  edge_iterator edge_begin() const {
    return EdgeIterator(node_begin(), (*node_begin()).edge_begin(), this);
  }

  /**
  *edge_end() returns the pointer to one past the last edge of the graph
  *@returns returns @EE, an edge_iterator 
  */

  edge_iterator edge_end() const {
    return EdgeIterator(node_end(), (*node_begin()).edge_end(), this);
  }

 private:

  //interal type for node
  struct node_info {

    Point P; //Cartesian location of the node
    node_value_type val;

  };

  //interal type for edge
  struct edge_info {
    
    size_type N_2_id; //id of node 2 of the edge
    edge_value_type e_val; 
  };

};


#endif // CME212_GRAPH_HPP
