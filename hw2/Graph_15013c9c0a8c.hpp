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
template <typename V,typename E>
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
  /** Synonym for Edge Value */
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

    /** Return this node's position. */
    Point& position() {
      return graph_->NodeInfo[uid_].NodeLocation;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(graph_->NodeInfo[uid_].index_>=0);
      assert(graph_->NodeInfo[uid_].index_<graph_->size());
      return graph_->NodeInfo[uid_].index_;
    }
    
    /** Return this node's uid, a number not necessarily in the range [0, graph_size). */
    size_type identifier() const {
      return uid_;
    }

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
    Node(const Graph* graph, size_type idx)
             : graph_(const_cast<Graph*>(graph)), uid_(graph->index2Uid[idx]) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return index2Uid.size();
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
    
    index2Uid.push_back(NodeInfo.size());
    NodeInfo.push_back(NodeStruct(position,val,index2Uid.size()-1));
    adjList.push_back(std::vector<std::pair<size_type,edge_value_type> >());
    return Node(this, (size()-1));
  }

  /** Remove a node from the graph and it's incident edges.
   * @param[in]: const Node n 
   * @post: new .size()=old .size()-return value
   * @post: no proxy nodes are invalidated
   * @post: current iterator and end iterator are invalidated 
   * @return: Number of nodes deleted
   * Complexity O(1) because it is O(max(degree)^2) where degree is limited above by a constant
   */
  size_type remove_node(const Node& n) {

    //std::cout<<"Deleting n uid_: "<<n.identifier()<<"\n"; 
    if(!has_node(n)) return 0; //Can't delete if not in the graph
    auto pred=[&](std::pair<size_type,edge_value_type> element){return element.first==n.identifier();};
    //for loop iterates over incident nodes n.degree()=O(1)
    for(auto it=adjList[n.identifier()].begin();it!=adjList[n.identifier()].end();++it){
       //find_if iterates through the incident node, incident nodes looking for original node n
       // n.degree()=O(1)
       auto eraseIt= find_if(adjList[(*it).first].begin(),adjList[(*it).first].end(),pred);
       //find if will always find our value since we checked that the graph has the node 
       //O(1) by same argument as above
       adjList[(*it).first].erase(eraseIt);
    }
    //O(1)
    adjList[index2Uid[n.index()]].clear();
   // std::cout<<"Clearing row: "<<adjList[index2Uid[n.index()]].size()<<"\n";
    //O(1)
    NodeInfo[index2Uid[size()-1]].index_=n.index(); 
    NodeInfo[index2Uid[n.index()]].NodeValue=node_value_type(); 
    NodeInfo[index2Uid[n.index()]].NodeLocation=Point(); 
    //std::cout<<"Setting end indx: "<<NodeInfo[index2Uid[size()-1]].index_<<"\n";
    //O(1)
    index2Uid[n.index()]=index2Uid[size()-1];
   // std::cout<<"Pop off indx2uid: "<<index2Uid.size()<<"\n";
    //O(1)
    index2Uid.pop_back(); 
    //std::cout<<"Pop off indx2uid: "<<index2Uid.size()<<"\n";
    return 1;
  }
  
  /** Remove a node from the graph and it's incident edges.
   * @param[in]: node_iterator it 
   * @pre: n_it is a valid iterator in range [.begin(),.end())
   * @post: new .size()=old .size()-1 if
   * @post: no proxy nodes are invalidated
   * @post: current iterator and end iterator are invalidated 
   * @return: iterator to next item to continue iteration over nodes 
   * Complexity O(1) because it is O(max(degree)^2) where degree is limited above by a constant
   */
  node_iterator remove_node(node_iterator n_it) {
    auto test=remove_node(*n_it);//removes value of n_it
    assert(test==1);// if n_it is a valid iterator then the node it points
    // to should be in the graph and able of being removed
    return n_it;//returns same iterator which now points to a different value 
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
    return (n.graph_ == this && (index2Uid[n.index()]==n.identifier()));
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

    /**Return the value stored in the graph for this proxy edge
    *@pre Proxy edge references a valid edge in graph_
    *@return a mutable value for this non-constant edge. 
    * Complexity O(1)  
    */
    edge_value_type& value(){
      size_type firstIndex; 
      size_type searchIndex; 
      firstIndex=std::min(uid1_,uid2_);
      searchIndex=std::max(uid1_,uid2_);
      auto pred=[&](std::pair<size_type,edge_value_type> element){return element.first==searchIndex;};
      auto foundIt= find_if(graph_->adjList[firstIndex].begin(),graph_->adjList[firstIndex].end(),pred);
      assert(foundIt!=(graph_->adjList[firstIndex]).end());//Valid edge is in adjList
      return (*foundIt).second;
    }

    /**Return the value stored in the graph for this proxy edge
    *@pre Proxy edge references a valid edge in graph_
    *@return a constant value for this constant edge. 
    * Complexity O(1)  
    */
    const edge_value_type& value() const{
      size_type firstIndex; 
      size_type searchIndex; 
      firstIndex=std::min(uid1_,uid2_);
      searchIndex=std::max(uid1_,uid2_);
      auto pred=[&](std::pair<size_type,edge_value_type> element){return element.first==searchIndex;};
      auto foundIt= find_if((graph_->adjList[firstIndex]).begin(),(graph_->adjList[firstIndex]).end(),pred);
      assert(foundIt!=(graph_->adjList[firstIndex]).end());//Valid edge is in adjList
      return (*foundIt).second;
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
    size_type size_=0; 
    for(auto it=adjList.begin();it!=adjList.end();++it){
       size_+=(*it).size(); 
    } 
    assert(size_%2==0);
    return size_/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(),i);

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
        if (adjList[a.uid_][i].first==b.uid_) return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val= edge_value_type()) {
    //Case 1: if the add we want to add already exists
    if (has_edge(a, b)) {return Edge(this, a.uid_, b.uid_);}
    // Case 2: if the edge does not exist
    else {
          // Add to the adjointed list
          adjList[a.uid_].push_back(std::make_pair(b.uid_,val));
          adjList[b.uid_].push_back(std::make_pair(a.uid_,val));
          return Edge(this, a.uid_, b.uid_);
          }
  }

  /** Remove a node from the graph and it's incident edges.
   * @param[in]: const Node n1, const Node n2 both valid nodes of the graph 
   * @post: new .num_edges=old .num_edges-return value
   * @post: no proxy edges are invalidated
   * @post: all edge_iterators [it,end()) are invalidated 
   * @return: Number of edges deleted (ex: deleted 1 edge that we stored twice)
  */
  size_type remove_edge(const Node& n1, const Node& n2){
     size_type response=0; 
     //To Do: Implement this instead using find_if (not wrong as is)
     for(auto it=adjList[n1.identifier()].begin();it!=adjList[n1.identifier()].end();++it){
        if((*it).first==n2.identifier()){
           adjList[n1.identifier()].erase(it);
           response=1;
           break;
        }
     }
     //To Do: Implement this instead using find_if (not wrong as is)
     for(auto it=adjList[n2.identifier()].begin();it!=adjList[n2.identifier()].end();++it){
        if((*it).first==n1.identifier()){
           adjList[n2.identifier()].erase(it);
           break;
        }
     }
     return response;//This will return 0 if n1,n2 is not a valid edge
     //No need to call has_edge with the above logic
  }

  /** Remove a node from the graph and it's incident edges.
   * @param[in]: const Edge e1  
   * @post: new .num_edges=old .num_edges-return value
   * @post: no proxy edges are invalidated
   * @post: all edge_iterators [it,end()) are invalidated 
   * @return: Number of edges deleted (ex: deleted 1 edge that we stored twice)
  */
  size_type remove_edge(const Edge& e1){
     return remove_edge(e1.node1(),e1.node2());
  }
  
  /** Remove a node from the graph and it's incident edges.
   * @param[in]: valid edge_iterator e_it   
   * @pre: e_it in range [.begin(),.end())  
   * @post: no proxy edges are invalidated
   * @post: all edge_iterators [it,end()) are invalidated  
   * @return: Number of edges deleted (ex: deleted 1 edge that we stored twice)
  */
  edge_iterator remove_edge(const edge_iterator e_it){
     Edge e1=*e_it;
     auto test=remove_edge(e1.node1(),e1.node2()); 
     assert(test==1);//valid edge_iterator must point to edge in graph
     return e_it; //Return the same iterator
     //which now points to a different edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    NodeInfo.clear();
    adjList.clear();
    index2Uid.clear();
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
       return Node(graph_,graph_->index2Uid[idx_]);
    }
    /* Increment the node iterator 
    * @pre *this is in the range [.node_begin(),.node_end()) 
    * @pre *this is in the range [.node_begin()+1,.node_end()] 
    * @return the iterator pointing to the next node 
    * Complexity O(1)
    */
    NodeIterator& operator++(){
       // increment our vector iterator
       ++idx_;
       return *this;
    } 
    /* Equality comparsion between node iterators
    *@param[in] OtherIt a nodeiterator to the "other" node 
    *@return if the two iterators are equal 
    * Complexity O(1)
    */
    bool operator==(const NodeIterator& OtherIt) const{
       return (graph_==OtherIt.graph_ &&idx_==OtherIt.idx_);
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph * graph_;
    size_type idx_;
    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type idxin)
           : graph_(const_cast<Graph*>(graph)), idx_(idxin) {
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
     node_iterator endit(this,index2Uid.size());
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
       return Edge(graph_,uid1_,graph_->adjList[uid1_][uid2_].first); 
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
       return Edge(graph_,graph_->index2Uid[node_],graph_->adjList[graph_->index2Uid[node_]][incident_].first);
    }
    /* Increments edge operator
    * @pre *this is in the range [.edge_begin(),.edge_end()) 
    * @pre *this is in the range [.edge_begin()+1,.edge_end()] 
    * @return this iterator when it points to the next edge or (*.edge_end())
    * Complexity O(1) on average over range [.edge_begin(),.edge_begin()), worst case complexity for single call O(num_edges)
    */
    EdgeIterator& operator++(){
       ++incident_;
       if(incident_==graph_->adjList[graph_->index2Uid[node_]].size()){
       incident_=0;
       ++node_;
       }
       while(node_<graph_->size() && graph_->adjList[graph_->index2Uid[node_]][incident_].first>graph_->index2Uid[node_]){
          ++incident_;
          if(incident_==graph_->adjList[graph_->index2Uid[node_]].size()){ 
             incident_=0;
             ++node_;
          }
       }
       return *this; 
    }
    /* Equality check between edge iterators
    * @param[in] otherIt the "other" iterator aka the RHS of the ==
    * @return true if this iterator and the other iterator point to the same element
    * Complexity O(1)
    */
    bool operator==(const EdgeIterator& otherIt) const{
       return (graph_==otherIt.graph_&&node_==otherIt.node_ &&incident_==otherIt.incident_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type node_;
    size_type incident_;
    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type uidin,size_type uidin2)
           : graph_(const_cast<Graph*>(graph)), node_(uidin),incident_(uidin2) {
      while(node_<graph_->index2Uid.size() && graph_->adjList[graph_->index2Uid[node_]][incident_].first>graph_->index2Uid[node_]){
          ++incident_;
          if(incident_==graph_->adjList[graph_->index2Uid[node_]].size()){ 
             incident_=0;
             ++node_;
          }
       }

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
      size_type node_=0;
      size_type incident_=0; 
      /* 
      while(adjList[index2Uid[node_]][incident_].first>index2Uid[node_] && node_<index2Uid.size()){
          ++incident_;
          if(incident_==adjList[index2Uid[node_]].size()){ 
             incident_=0;
             ++node_;
          }
       }
      */
     return EdgeIterator(this,node_,incident_);
  }
  /** Return an iterator to 1 past the last edge in the graph
  *@return an iterator to 1 past the last node the input node is connected to via an edge 
  * if the graph has no edges .edge_begin()==.edge_end() 
  * Complexity O(1) 
  */
  edge_iterator edge_end() const{
     if(size()==0) return EdgeIterator(this,size(),0);
     return EdgeIterator(this,size(),0);
  }

 private:
    struct NodeStruct{
      Point NodeLocation; //stores the nodes(point)
      node_value_type NodeValue; //stores the nodes(value)
      size_type index_;//index of node
      NodeStruct(Point pt, node_value_type val,size_type ind):
         NodeLocation(pt),NodeValue(val),index_(ind){}
      NodeStruct():NodeLocation(),NodeValue(),index_(-1){};
    };
    std::vector<NodeStruct> NodeInfo; //vector to store Node position and value
    std::vector<size_type> index2Uid; //vector to store index mapping 
    std::vector<std::vector<std::pair<size_type,edge_value_type>>> adjList; //Store the nodes in adjointed list format
    //aka for each node, we all know all nodes to which it is connected
};

#endif // CME212_GRAPH_HPP
