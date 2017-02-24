
#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <tuple>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
// Defining template
template <typename V, typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Define an internal vector to store the nodes
  

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V,E>;
  typedef V node_value_type ;
  typedef E edge_value_type ;

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
  using uid_type = unsigned;
  using idx_type = signed;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : internal_nodes(), idx2uid_() {
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
    } // Construction of invalid node

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->internal_nodes[uid_].position_;
    }

    // Modifiable Node position
    Point& position(){
      // HW2: YOUR CODE HERE
      return graph_->internal_nodes[uid_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). 
    * @pre this is a valid node
    */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->internal_nodes[uid_].idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /*Return value associated with a node 
    */
    node_value_type& value() {
      return graph_->internal_nodes[uid_].value_;
    }

    const node_value_type& value() const {
      return graph_->internal_nodes[uid_].value_;
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->graph_ == n.graph_)
        if (this->uid_ == n.uid_)
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
      // Comparing nodes in same graph 
      if (this->graph_ == n.graph_)
        if (this->uid_ < n.uid_)
          return true;

      // Comparing nodes in different graphs
      return (this->graph_ < n.graph_);
    }

    /** Return the degree of this node
    */
    size_type degree() const {
      return graph_->internal_nodes[uid_].adjacency_.size();
    }

    /** Return an incident iterator that points past the last edge incident on 
    a node
    */
    IncidentIterator edge_end() const {
      return IncidentIterator(graph_,uid_,
                              graph_->internal_nodes[uid_].adjacency_.size());
    }

     /** Return an incident iterator to the first edge incident on a node
    */
    IncidentIterator edge_begin() const {
      // If degree is not zero
      if (degree()!= 0) return IncidentIterator(graph_,uid_,0);
      return edge_end(); // Else return edge_end()
    }

    

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_; // Pointer to Graph class
    size_type uid_; // ID number    
    Node(const Graph* graph, size_type uid)
       : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    } // Valid Constructor

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return idx2uid_.size();
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
  // Default constructor of node_value_type is used to assign a default value
  Node add_node(const Point& position,const node_value_type& value =node_value_type()){
    // HW0: YOUR CODE HERE
    // Adding nodes to the internal_nodes vector
    internal_nodes.push_back(NodeInfo(position,value,
                              idx2uid_.size(), 
                              std::vector<AdjacencyInfo> ())); 
    idx2uid_.push_back(internal_nodes.size()-1); // index to uid mapping
    return Node(this,internal_nodes.size()-1);   // Return node to the user
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_ == this)
      if (n.graph_->internal_nodes[n.uid_].idx_!=-1)
        if (n.graph_->internal_nodes[n.uid_].idx_<(idx_type)n.graph_->num_nodes())
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
    // HW0: YOUR CODE HERE
    return Node(this,idx2uid_[i]);  // Node with index i
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
      return Node(this->graph_,this->node1_);      // Node 1
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_,this->node2_);     // Node 2
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const { 
      return (std::make_tuple(std::min(node1(),node2()),
                              std::max(node1(),node2())) ==
              std::make_tuple(std::min(e.node1(),e.node2()),
                              std::max(e.node1(),e.node2()))); 
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Comparing edges in same graph || different graphs
      return ((std::make_tuple(this->graph_,this->uid_)==
              std::make_tuple(e.graph_,e.uid_)) || 
              (this->graph_ < e.graph_));
    }

    // Get length of edges
    double length() const{
      return norm(node1().position()-node2().position());
    }

    // Access value of edges
    edge_value_type& value(){
      size_type i = 0;
      while (graph_->internal_nodes[node1().uid_].adjacency_[i].node_id_ !=
             node2().uid_) {
        ++i;
      }      
      return graph_->internal_nodes[node1().uid_].adjacency_[i].value_;
    };

    const edge_value_type& value() const{
      size_type i = 0;
      while (graph_->internal_nodes[node1().uid_].adjacency_[i].node_id_ !=
             node2().uid_) {
        ++i;
      }      
      return graph_->internal_nodes[node1().uid_].adjacency_[i].value_;
    };

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    size_type uid_,node1_,node2_; // edge id, node ids
    Edge(const Graph* graph, size_type uid, size_type node1, size_type node2)
       : graph_(const_cast<Graph*>(graph)), uid_(uid),node1_(node1),node2_(node2){
    } // Valid Constructor

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return std::distance(edge_begin(),edge_end());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // Using edge iterator to find the edge corresponding to index i
    for (auto ei=edge_begin();ei!=edge_end();++ei){
      if ((*ei).uid_==i)
        return *ei;
    }
    return Edge();       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Run through the shorter of the two adjacency lists
    if (a.degree()<=b.degree()){
      for (size_type i = 0; i< internal_nodes[a.uid_].adjacency_.size(); i++){
        if (a.graph_->internal_nodes[a.uid_].adjacency_[i].node_id_ == b.uid_)
          return true;
      }
      return false; 
    } 
    else{
      for (size_type i = 0; i<b.degree(); i++){
        if (b.graph_->internal_nodes[b.uid_].adjacency_[i].node_id_ == a.uid_)
          return true;
      }
      return false; 
    }      
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
  Edge add_edge(const Node& a, const Node& b, 
                const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // Return existing edge
    for (size_type i=0; i<a.degree();++i)
      if (a.graph_->internal_nodes[a.uid_].adjacency_[i].node_id_ == b.uid_)
        return Edge(this,a.graph_->internal_nodes[a.uid_].adjacency_[i].edge_id_,
                    a.uid_,b.uid_);       

    // Pushing back into the adjacency list maintained by each node
    internal_nodes[a.uid_].adjacency_.push_back(AdjacencyInfo(b.uid_,num_edges(),
                                                              value));
    internal_nodes[b.uid_].adjacency_.push_back(AdjacencyInfo(a.uid_,num_edges()-1,
                                                              value));
    return Edge(this,num_edges()-1,a.uid_,b.uid_);          
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    internal_nodes.clear();
    idx2uid_.clear();
    return;
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

    /** Deference node iterator 
    *@pre NodeIterator points to a valid Node object of this graph
    *@post returns Node object, such that
    *Node.graph_ = NodeIterator.graph_;
    *Node.uid_ = NodeIterator.uid_;
    */
    Node operator*() const {
      // Dereferencing operator
      return Node(graph_,graph_->idx2uid_[idx_]);
    } 

    /** Increment operator
    *@pre NodeIterator points to a valid Node object of this graph
    *@post returns NodeIterator, such that
    *NodeIterator_return.graph_ = NodeIterator.graph_;
    *NodeIterator_return.uid_ = NodeIterator.uid_ +1;
    */
    NodeIterator& operator ++(){
      // Increment operator
      ++idx_;
      return *this;
    }

    /** Test whether this NodeIterator is equal to @a node_iterator in a global 
    * order.
    */    
    bool operator ==(const NodeIterator& node_iterator) const {
      return ((graph_ == node_iterator.graph_)
        && (idx_ == node_iterator.idx_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_; // Pointer to Graph class
    size_type idx_; // ID number    
    NodeIterator(const Graph* graph, size_type idx)
       : graph_(const_cast<Graph*>(graph)), idx_(idx) {
    } // Valid Constructor

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns a NodeIterator to the start of the container containing nodes
  */
  NodeIterator node_begin() const {
    return NodeIterator(this,0);
  }

  /** Returns a NodeIterator to one past the end of the container
  * containing nodes
  */
  NodeIterator node_end() const {
    return NodeIterator(this,num_nodes());
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

    /** Deference IncidentIterator
    *@pre IncidentIterator points to a valid Edge object of this graph
    *@post returns Edge object
    */
    Edge operator *() const {
      return Edge(graph_,graph_->internal_nodes[node_id_].adjacency_[id_].edge_id_, 
        node_id_, graph_->internal_nodes[node_id_].adjacency_[id_].node_id_);
    }

    /** Increment operator
    *@pre IncidentIterator points to a valid Edge object of this graph
    *@post returns IncidentIterator
    */
    IncidentIterator& operator ++() {
      ++id_;
      return *this;
    }

    /** Test whether this IncidentIterator is equal to @a iit in a global order.
    */
    bool operator ==( const IncidentIterator & iit ) const {
      return ((graph_ == iit.graph_) && (node_id_ == iit.node_id_)
        && (id_ == iit.id_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    uid_type node_id_, id_; // node id, item number in adjacency list of node_id
    IncidentIterator(const Graph* graph, size_type node_id, size_type id)
       : graph_(const_cast<Graph*>(graph)), node_id_(node_id), id_(id) {
    } // Valid Constructor
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
    EdgeIterator()  {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Deference node iterator 
    *@pre EdgeIterator points to a valid Edge object of this graph
    *@post returns edge object
    */
    Edge operator*() const {
      // Dereferencing operator
      auto node_id = graph_->idx2uid_[node1_]; // Get uid_ of node
      return Edge(graph_,graph_->internal_nodes[node_id].adjacency_[id_].edge_id_,
                  graph_->idx2uid_[node1_],
                  graph_->internal_nodes[graph_->idx2uid_[node1_]].
                                                      adjacency_[id_].node_id_);
    } 

    /** Increment operator
    *@pre EdgeIterator points to a valid Edge object of this graph
    *@post returns EdgeIterator that points to the next valid edge of the graph
    */
    EdgeIterator& operator ++(){
      ++id_;
      // Go through nodes of the graph
      while ((size_type) node1_< graph_->num_nodes()){
        // Go throuch each incident "unique" edge
        while (id_< Node(graph_,graph_->idx2uid_[node1_]).degree()){
          // Check if we get a new edge
          auto adj_node = graph_->internal_nodes[graph_->idx2uid_[node1_]].
                                                      adjacency_[id_].node_id_;
          if (graph_->idx2uid_[node1_] < adj_node) { // To avoid repeating edges
            return *this;
          }
        ++id_;
        }
        // If we reach end of adjacency list of a particular node, 
        // increase index and start again
        ++node1_;
        id_=0;
      }      
      return *this;   
    }

    /** Test whether this EdgeIterator is equal to @a edge_iterator in a global order.
    */
    bool operator ==(const EdgeIterator& edge_iterator) const {
      return ((graph_ == edge_iterator.graph_)
        && (id_ == edge_iterator.id_)
        && (node1_ == edge_iterator.node1_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_; // Pointer to Graph class
    idx_type node1_; // Index of One node of a edge  
    uid_type id_; // ID number on adjacency list of node1
      
    EdgeIterator(const Graph* graph, idx_type node1, uid_type id)
       : graph_(const_cast<Graph*>(graph)), node1_(node1), id_(id){
    }; // Valid Constructor

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const


  /** Returns a EdgeIterator to one past the end of the container
  * containing edges
  */
  EdgeIterator edge_end () const {
    return EdgeIterator(this,num_nodes(),0);
  }

   /** Returns a EdgeIterator to the start of the container containing edges
  */
  EdgeIterator edge_begin () const {    
    for (size_type node_idx=0; node_idx < num_nodes(); ++node_idx){
      if (Node(this, idx2uid_[node_idx]).degree()>0)
        return EdgeIterator(this,node_idx,0);
    }
    return edge_end();
  }

  /** Remove a node and all its incident edges from the graph.
   * @pre @a n is a valid node of this graph
   *
   * @post has_node(@a n) == false
   * @post new num_nodes() == old num_nodes()-1
   * @post new has_edge(@a n, @a n2) == false for all n2 such that, 
   * old has_edge(@a n, @a n2) == true 
   * @post new num_edges() == old num_edges()- old @a n.degree()
   * 
   * @post new @a n1.index() == old @a n1.index()-1 for all n1 such that,
   * old @a n1.index() > @a n.index()
   * @result 0
   *
   * Complexity: O(n.degree()+num_nodes()) 
  */
  size_type remove_node(const Node& n){
    // Remove all incident edges
    size_type i = 0;
    while (i < n.degree()){
      remove_edge(n, Node(this, internal_nodes[n.uid_].adjacency_[i].node_id_));
    }

    // Decrement indices of all other nodes after 'n' by 1
    auto idx = (size_type) n.index();
    for (auto i=idx+1; i<num_nodes(); ++i){
      internal_nodes[idx2uid_[i]].idx_ -= 1;
    }

    // Remove n form the idx2uid and change its index in internal_nodes data
    internal_nodes[idx2uid_[idx]].idx_=-1;
    idx2uid_.erase(idx2uid_.begin()+idx);
    return 0;
  }

/** Remove a node and all its incident edges from the graph, using a node iterator
  * @pre @a n_it points to a valid node object of this graph
   * @result NodeIterator that points to the node @a n1, such that
   * old @a n1.index() == @a (*n_it).index()+1
  */
  node_iterator remove_node(node_iterator n_it){
    auto n = *n_it;
    remove_node(n);
    return n_it;
  }

  /** Remove an edge from the graph
   * @param n1,n2 Two ends of the edge
   * 
   * @pre @a n1, n2 are valid nodes of the graph
   * @post has_edge(@a n1, @a n2) == false
   * @post new num_edges() == old num_edges()- 1
   * 
   * @post new @a e.uid_ == old @a e.uid_-1 for all edge @a e such that,
   * old @a e.uid_ > @a e_remove.uid_ where e_remove.node1() == n1, 
   * e_remove.node2() == n2
   * @post new n1.degree() == old n1.degree-1
   * @post new n2.degree() == old n2.degree-1
   * @post Node n1 and edge info is deleted from the adjacency list of n2 
   * and vice versa
   * @result 0
   *
   * Complexity: O(n.degree()+num_nodes()) 
  */
  size_type remove_edge(const Node& n1, const Node& n2){
    // Remove from adjacency list of n1
    
    size_type edge_id = 0;
    if (has_edge(n1,n2)){
      for (size_type i=0; i< n1.degree(); ++i){
        if (internal_nodes[n1.uid_].adjacency_[i].node_id_== n2.uid_){
          edge_id = internal_nodes[n1.uid_].adjacency_[i].edge_id_;
          internal_nodes[n1.uid_].adjacency_.erase(internal_nodes[n1.uid_].
                                                    adjacency_.begin()+i);
          break;
        }
      }
      // Remove from adjacency list of n2
      for (size_type i=0; i< internal_nodes[n2.uid_].adjacency_.size(); ++i){
        if (internal_nodes[n2.uid_].adjacency_[i].node_id_==n1.uid_){
          internal_nodes[n2.uid_].adjacency_.erase(internal_nodes[n2.uid_].
                                                    adjacency_.begin()+i);
          break;
        }
      }

      // Renumber all the edges greater than the removed edge
      for (auto ni= node_begin();ni!= node_end();++ni){
        for (size_type i=0; i< internal_nodes[(*ni).uid_].adjacency_.size(); ++i){
          if (internal_nodes[(*ni).uid_].adjacency_[i].edge_id_ > edge_id){
            internal_nodes[(*ni).uid_].adjacency_[i].edge_id_ -= 1;
          }          
        }
      }
    }
    (void) edge_id;
    return 0;
  }

  // Remove edge from a graph
  /* @a param e Edge of the graph
  * @result 0
  */
  size_type remove_edge(const Edge& e){
    remove_edge(e.node1(),e.node2());
    return 0;
  }

  /* Remove edge using an edge iterator
  * @param e_it Edge iterator
  * @result EdgeIterator that points to the edge @a e1, such that
  * old @a e1.uid_ == @a (*e_it).uid_+1
  */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return e_it;
  }


 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Details of adjacency list
  struct AdjacencyInfo{
    AdjacencyInfo() {};
    AdjacencyInfo(uid_type node_id, uid_type edge_id, edge_value_type value)
                  : node_id_(node_id),edge_id_(edge_id), value_(value){};

    uid_type node_id_;
    uid_type edge_id_;
    edge_value_type value_;
  };

  // Details of a particular node
  struct NodeInfo{
    NodeInfo() {};
    NodeInfo(Point position, node_value_type value, idx_type idx,
             std::vector<AdjacencyInfo> adjacency)
      : position_(position), value_(value),idx_(idx),adjacency_(adjacency) {};

    Point position_; // Position
    node_value_type value_; // Value
    idx_type idx_; // Index
    std::vector<AdjacencyInfo> adjacency_;  // Adjacency list  
  };

  std::vector<NodeInfo> internal_nodes;
  std::vector<uid_type> idx2uid_;
};

#endif // CME212_GRAPH_HPP