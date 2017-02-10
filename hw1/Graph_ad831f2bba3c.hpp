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
  /** Synonym for Node Value */
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
  class Node:private totally_ordered<Node> {
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
      return graph_->NodeInfo[uid_].NodeLocation;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /**Return the value stored in the graph for this proxy node
    *@pre Proxy node references a valid node in graph_
    *@return a mutable value for this non-constant node. 
    * Complexity O(1)  
    */
    node_value_type& value(){
      return graph_->NodeInfo[uid_].NodeValue;
    }
    /**Return the value stored in the graph for this constant proxy node
    *@pre Proxy node references a valid node in graph_
    *@return a constant value for this constant node.
    * Complexity O(1) 
    */
    const node_value_type& value() const{
      return graph_->NodeInfo[uid_].NodeValue;
    }
    /** Return the number of edges connect to the node
    *@pre Proxy node references a valid node in graph_
    *@return the number of edges connected to this node.
    * complexity O(1)
    */
    size_type degree() const{
      return graph_->adjList[uid_].size(); 
    } 
    /** Return iterator to the first node the input node is connected to via an edge
    *@return an iterator to the first node the input node is connected to via an edge 
    *@post if this is an isolated node .edge_begin()==.edge_end() 
    * Complexity O(1) 
    */
    incident_iterator edge_begin() const{
       return IncidentIterator(graph_,uid_,0);
    }
    /*Return an iterator to 1 past the last node the input node is connected to via an edge 
    *@return an iterator to 1 past the last node the input node is connected to via an edge 
    *@post if this is an isolated node .edge_begin()==.edge_end() 
    * Complexity O(1) 
    */
    incident_iterator edge_end() const{
       return IncidentIterator(graph_,uid_,graph_->adjList[uid_].size());
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      //compare 2 graphs - same uid and node on same graph
      return ((n.graph_ ==  graph_) && (n.uid_ == uid_));
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
      //There are 2 cases - if on different graphs and if on same graph
      if (n.graph_ != graph_) return (n.graph_ < graph_);
      else return (n.uid_ < uid_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // A node is defined by its pointer to the graph, and its id
    Graph* graph_;
    // This element's unique identification number
    size_type uid_;
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
             : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return NodeInfo.size();
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
    //when we add a node, we need to add to the - vector of nodes, but also
    // to the edges via the adjoined List
    
    NodeInfo.push_back(NodeStruct(position,val));
    adjList.push_back(std::vector<size_type>());
    return Node(this, (size()-1));
  }
/*
  Node add_node(const Point& position) {
    //when we add a node, we need to add to the - vector of nodes, but also
    // to the edges via the adjoined List
    NodeLocations.push_back(position);
    adjList.push_back(std::vector<size_type>());
    return Node(this, (size()-1));
  }
*/

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_ == this && (n.uid_<=size()));
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
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
      //conditions: belong to same graph and see if the nodes id's match (check both directions)

      return (((uid1_== e.uid1_) && (uid2_ == e.uid2_) && (graph_ == e.graph_)) ||
              ((uid1_== e.uid2_) && (uid2_ == e.uid1_) && (graph_ == e.graph_)));

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //IDEA - look, in this order, at both graphs and nodes
      //If different graphs, look directly at the graphs id's
      // If same graphs, check first by minumum, and if same minimum, by maximum
      //If different graph, simiply look at their graphs ids
      if (graph_ != e.graph_) {return graph_ < e.graph_;}
      // if on same graph, check the min and max between 2 nodes
      else 
         {
            if (std::min(e.uid1_, e.uid2_) != std::min(uid1_,uid2_))
                {return std::min(e.uid1_, e.uid2_) < std::min(uid1_,uid2_);}
            else
                {return std::max(e.uid1_, e.uid2_) < std::max(uid1_,uid2_);}
         }

    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // An Edge is defined by its pointer to a graph, and its 2 nodes(id's of nodes)
    Graph* graph_;
    // This element's unique identification number
    size_type uid1_;
    size_type uid2_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type uid1, size_type uid2)
           : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
    }

      //
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
  Edge edge(size_type i) const {
    size_type id1 = edges[i].first;
    size_type id2 = edges[i].second;
    return Edge(this, id1, id2);

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //check via the adjointed list, where the adjointed list is like a vector
    // of vector; to remember later in terms of x[i][j], in terms of the search

    // Note could use std vector find which implements the same thing
    for (size_type i=0; i<adjList[a.uid_].size(); i++) {
        if (adjList[a.uid_][i]==b.uid_) return true;
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
    //Case 1: if the add we want to add already exists
    if (has_edge(a, b)) {return Edge(this, a.uid_, b.uid_);}
    // Case 2: if the edge does not exist
    else {
          // Add to the adjointed list
          adjList[a.uid_].push_back(b.uid_);
          adjList[b.uid_].push_back(a.uid_);
          // Add also to the pairs, in a sorted way
          if (a.uid_ < b.uid_) edges.push_back(std::make_pair(a.uid_, b.uid_));
          else edges.push_back(std::make_pair(b.uid_, a.uid_));

          return Edge(this, a.uid_, b.uid_);
          }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges.clear();
    NodeInfo.clear();
    adjList.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private totally_ordered<NodeIterator> {
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
    /* Dereference Node iterator 
    * @pre *this is in the range [.node_begin(),.node_end()) 
    * @return the node this iterator points to
    * Complexity O(1)
    */
    Node operator*() const{
       return Node(graph_,uid_);
    }
    /* Increment the node iterator 
    * @pre *this is in the range [.node_begin(),.node_end()) 
    * @pre *this is in the range [.node_begin()+1,.node_end()] 
    * @return the iterator pointing to the next node 
    * Complexity O(1)
    */
    NodeIterator& operator++(){
       // increment our vector iterator
       ++uid_;
       return *this;
    } 
    /* Equality comparsion between node iterators
    *@param[in] OtherIt a nodeiterator to the "other" node 
    *@return if the two iterators are equal 
    * Complexity O(1)
    */
    bool operator==(const NodeIterator& OtherIt) const{
       return (graph_==OtherIt.graph_ &&uid_==OtherIt.uid_);
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph * graph_;
    size_type uid_;
    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type uidin)
           : graph_(const_cast<Graph*>(graph)), uid_(uidin) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
    /* Return a node_iterator to the first node 
    * @return a node iterator to the first node if g.size()>0
    * else node_begin()==node_end() 
    * Complexity O(1)
    */
  node_iterator node_begin() const{
     node_iterator beginit(this,0);
     return beginit;
  }
    /* Return a boundary to the node iterators 
    * @return an iterator to an invalid node that is signifies the end of valid nodes
    * Complexity O(1)
    */
  node_iterator node_end() const{
     node_iterator endit(this,NodeInfo.size());
     return endit; 
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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
    /* Dereference incident iterator
    * @pre *this is in the range [.edge_begin(),.edge_end()) 
    * @return The edge this iterator points to
    * Complexity O(1)
    */
    Edge operator*() const{
       return Edge(graph_,uid1_,graph_->adjList[uid1_][uid2_]); 
    }
    /* Increments incident operator
    * @pre *this is in the range [.edge_begin(),.edge_end()) 
    * @pre *this is in the range [.edge_begin()+1,.edge_end()] 
    * @return this iterator when it points to the next edge or (*.edge_end())
    * Complexity O(1)
    */
    IncidentIterator& operator++(){
       ++uid2_; 
       return *this; 
    }
    /* Equality check between incident iterators
    * @param[in] otherIt the "other" iterator aka the RHS of the ==
    * @return The true if this iterator and the other iterator point to the same element
    * Complexity O(1)
    */
    bool operator==(const IncidentIterator& otherIt) const{
       return (graph_==otherIt.graph_ && uid1_==otherIt.uid1_ && uid2_==otherIt.uid2_);
    }

   private:
    friend class Graph;
    Graph * graph_;
    size_type uid1_;
    size_type uid2_;
    // HW1 #3: YOUR CODE HERE
    /** Private Constructor */
    IncidentIterator(const Graph* graph, size_type uid1,size_type uid2)
           : graph_(const_cast<Graph*>(graph)), uid1_(uid1),uid2_(uid2) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private totally_ordered<EdgeIterator> {
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
    /* Dereference edge iterator
    * @pre *this is in the range [.edge_begin(),.edge_end()) 
    * @return The edge this iterator points to
    * Complexity O(1)
    */
    Edge operator*() const{
       return Edge(graph_,(graph_->edges[uid_]).first,(graph_->edges[uid_]).second);
    }
    /* Increments edge operator
    * @pre *this is in the range [.edge_begin(),.edge_end()) 
    * @pre *this is in the range [.edge_begin()+1,.edge_end()] 
    * @return this iterator when it points to the next edge or (*.edge_end())
    * Complexity O(1)
    */
    EdgeIterator& operator++(){
       ++uid_;
       return *this; 
    }
    /* Equality check between edge iterators
    * @param[in] otherIt the "other" iterator aka the RHS of the ==
    * @return true if this iterator and the other iterator point to the same element
    * Complexity O(1)
    */
    bool operator==(const EdgeIterator& otherIt) const{
       return (graph_==otherIt.graph_&&uid_==otherIt.uid_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type uid_;
    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type uidin)
           : graph_(const_cast<Graph*>(graph)), uid_(uidin) {
    }
     
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return an iterator to the first edge in the graph
  *@return an iterator to the first edge in the graph 
  * if the graph has no edges .edge_begin()==.edge_end() 
  * Complexity O(1) 
  */
  edge_iterator edge_begin() const{
     return EdgeIterator(this,0);
  }
  /** Return an iterator to 1 past the last edge in the graph
  *@return an iterator to 1 past the last node the input node is connected to via an edge 
  * if the graph has no edges .edge_begin()==.edge_end() 
  * Complexity O(1) 
  */
  edge_iterator edge_end() const{
     return EdgeIterator(this,edges.size());
  }

 private:
    struct NodeStruct{
      Point NodeLocation; //stores the nodes(point)
      node_value_type NodeValue; //stores the nodes(value)
      NodeStruct(Point pt, node_value_type val):
         NodeLocation(pt),NodeValue(val){}
    };
    std::vector<NodeStruct> NodeInfo; //vector to store Node position and value
    std::vector<std::vector<size_type>> adjList; //Store the nodes in adjointed list format
    //aka for each node, we all know all nodes to which it is connected
    std::vector<std::pair<size_type, size_type>> edges; //store simple pairs of nodes, aka
    //"simple edges" - (i,j)
};

#endif // CME212_GRAPH_HPP
