#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <utility>
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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)


  //predeclare struct for nodes and edges
  struct Proxy_Node;
  struct Proxy_Edge;

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

  //template the class
 // using node_value_type = V;
//  using edge_value_type = E;
    typedef V node_value_type;
    typedef E edge_value_type;

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
    Node() {

      // invalid node
      // HW0: YOUR CODE HERE

    }

    /** Return this node's position. */
    /** Position is not modifiab;e */
    const Point& position() const {
      assert(valid());
      //return Point in the uid position of Graphs vector
      return G->v[uid].p;



    }

    /**Position is modifiable */
    Point& position(){
     assert(valid());
      return G->v[uid].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE

      return uid;
    }



    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning

      if(uid == n.uid && G == n.G) return true;


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
      (void) n;           // Quiet compiler warning

      if(G == n.G){
        if(uid < n.uid) return true;

      }
      else{
        if(G < n.G) return true;

      }


      return false;
    }

    //HW1, part1,to use template
    /** Return the value of a Point
     *  @return the value of @a Point

    */
    node_value_type& value(){
        assert(valid());
        return G->v[uid].value;

    }
    /** Return the value of a constant Point
     *   @return the value of @a Point

         @pre @a Point calling this func must be const
     */
    const node_value_type& value() const{
      assert(valid());
      return G->v[uid].value;

    }


    //hw1 part 3
    /** Return the degree of a Node
     *  @return the number of edges from a @a Node

     *  @post calling object must remain unchanged
     */
    size_type degree() const{
      assert(valid());
      return G->v[uid].deg;
    }

    /** Return an IncidentIterator at the beginning of a nodes
          adjacent neighbors

     *  @return an IncdientIterator at the beginning of a nodes
          adjacent neighbors

     *  @post calling object must remain unchanged
     */
    IncidentIterator edge_begin() const{

      return IncidentIterator(G,uid,0);

    }


     /** Return an IncidentIterator at the end of a nodes
          adjacent neighbors

     *  @return an IncdientIterator at the end of a nodes
          adjacent neighbors

     *  @post calling object must remain unchanged
     */
    IncidentIterator edge_end() const{

      return IncidentIterator(G, uid,G->adj[uid].size());

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects


      //check for valid node
      bool valid() const{

        return uid>=0 && uid < G->v.size();


      }

      //points to the Graoh container
      Graph* G;
      //this tracks the idx for the point located in the idx position of
      //the vector in Graph
      size_type uid;
      //degree of this node

      //private constuctor, Graph can make a valid node
      Node(const Graph* g, size_type u)
        : G(const_cast<Graph*>(g)), uid(u){

      }


  };


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE

    return v.size();

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
/*i
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE


    size_type old = num_nodes();
    Proxy_Node new_node;
    new_node.p = position;
    v.emplace_back(new_node);
    size_type current = num_nodes();
    assert( current == old +1);
    assert(Node(this,current-1).index() == old);
    (void) position;       //Quiet compiler warning
    return Node(this, current-1);        // valid nlode



  }
*/
  //hw1, part 1, to use template
  Node add_node(const Point& position, const node_value_type& nv = node_value_type()){

      size_type old = num_nodes();
      Proxy_Node new_node;
      new_node.p = position;
      new_node.value = nv;
      v.emplace_back(new_node);
      adj.emplace_back(0);//allocate space in adj list for when adding edges
      size_type current = num_nodes();
      assert(current ==old+1);
      assert(Node(this,current-1).index()==old);
      v[current-1].pn_uid = current-1;
      return Node(this,current-1);




  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning

    //node is in the graph if its uid is less than the
    //number of nodes - 1
    size_type current = num_nodes();
    if (n.uid <= current - 1) return true;


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
    (void) i;             // Quiet compiler warning


    assert(i >=0 && i < num_nodes());
    assert(Node(this,i).index()== i );

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

      // invalid Edge

      // HW0: YOUR CODE HERE

    }

    /** Returns the edge length between two adjacent nodes

     */
    double length() const{

      return norm(node1().position()-node2().position());

    }


    /** Return edge value

     */
    edge_value_type& value(){

      Node n2 = node2();
      size_type found=0;
      for(unsigned i=0; i<e_G->adj[n1_uid].size();++i){
        size_type n= std::get<0>(e_G->adj[n1_uid][i]);
        if(n==n2.index()){
          found = i;
        }
      }

      return std::get<1>(e_G->adj[n1_uid][found]).value;
    }

    /**Return edge value

     @post calling object must remain unchanged
     */
    const edge_value_type& value() const{

      Node n2 = node2();
      size_type found =0;
      for(auto i=0; i<e_G->adj[n1_uid].size();++i){
        size_type n= std::get<0>(e_G->adj[n1_uid][i]);
        if(n==n2.index())
          found=i;

      }
      return std::get<1>(e_G->adj[n1_uid][found]).value;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE

      return Node(e_G,n1_uid);      // valid Node


   }
    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE

      return Node(e_G,N2_uid);      // valid Node



    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      //check the incoming edge nodes with its own uid
      if((e_G==e.e_G) &&
         (std::min(n1_uid,N2_uid)==std::min(e.n1_uid,e.N2_uid)) &&
         (std::max(n1_uid,N2_uid)==std::max(n1_uid,N2_uid))){
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
      if(e_G==e.e_G){
        if((std::min(n1_uid,N2_uid) < std::min(e.n1_uid,e.N2_uid)) &&
           (std::max(n1_uid,N2_uid) < std::max(e.N2_uid,e.n1_uid))){
           return true;
        }
      }
      else{
        if(e_G < e.e_G) return true;

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

    //point back to Graph
    Graph* e_G;
    //keep track of idx nodes this  edge vector has
    size_type n1_uid;
    size_type N2_uid;
    //Private constructor for Graph to make valid Edge
    Edge(const Graph* e_g, size_type n1_u, size_type N2_u)
      : e_G(const_cast<Graph*>(e_g)), n1_uid(n1_u), N2_uid(N2_u){

    }


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE

    return e.size();


  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

    assert(i >=0 && i < num_edges());
    std::pair<size_type,size_type> temp  = this->e[i];

    return Edge(this,std::get<0>(temp),std::get<1>(temp));        // Invalid Edge



  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    //loop through vector containing proxy_edges to get their nodes
    for(size_type i=0; i < e.size(); i++){
      std::pair<size_type,size_type> temp= e[i];
      size_type n1_temp = std::get<0>(temp);
      size_type n2_temp = std::get<1>(temp);
      if((std::min(n1_temp,n2_temp)==std::min(a.uid,b.uid)) &&
         (std::max(n1_temp,n2_temp)==std::max(a.uid,b.uid))){
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
  
    bool has_e = has_edge(a,b);
    if(has_e){
       return Edge(this,a.uid,b.uid);
      }

    else{
      std::pair<size_type,size_type> temp(a.uid,b.uid);
      this->e.push_back(temp);//add to unique e
      Proxy_Edge new_edge1;
      Proxy_Edge new_edge2;
      new_edge1.n2_uid=b.uid;
      new_edge2.n2_uid=a.uid;
      std::pair<size_type,Proxy_Edge> pair1(b.uid,new_edge1);
      std::pair<size_type,Proxy_Edge> pair2(a.uid,new_edge2);
      adj[a.uid].push_back(pair1);
      v[a.uid].deg +=1; //add to the degree
      adj[b.uid].push_back(pair2);
      v[b.uid].deg+=1;
      //track adj idx,fill up the adj map
      v[a.uid].adj_uid[b.uid] = adj[b.uid].size()-1;
      v[b.uid].adj_uid[a.uid] = adj[a.uid].size()-1;


      return Edge(this,a.uid,b.uid);
    }



  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

   v.clear();
   e.clear();

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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return a Node
     *  @return @a Node

     *  @post   @a Node must remain unchanged
     */
    Node operator*() const{
       return Node(this->g,curr);

    }

    /** Return a NodeIterator from incrementation
     *  @return @a NodeIterator from incrementation

     */
    NodeIterator& operator++(){
      this->curr = curr + 1;
      return *this; 

    }
    /** Return a boolean from equating NodeIterators
     *  @param[in] ni  NodeIterator
     *  @return @a bool from equating NodeIterators

     *  @post @a ni remains unchanged
     */
    bool operator==(const NodeIterator& ni) const{

      return curr == ni.curr;

    }




   private:
    friend class Graph;
    //point back to graph
    Graph* g;
    //keep track of current node we have
    size_type curr;
    // private constructor for NodeIterator
    NodeIterator(const Graph* ng, size_type c)
      : g(const_cast<Graph*>(ng)), curr(c){

    }


  };

  //NodeIterator funcs for Graph

  /** Return a NodeIterator at the head of a Graph
   *  @return @a NodeIterator point to head of Graph

   *  @post   @a Graph must not be changed
   */
  NodeIterator node_begin() const{
    return NodeIterator(this,0);

  }

  /** Return a NodeIterator at the end of a Graph
   *  @return @a NodeIterator point to end of  Graph

   *  @post   @a Graph must not be changed
   */
  NodeIterator node_end() const{
     return NodeIterator(this,size());

  }


  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

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

    /** Return an Edge
     *  @return an Edge

     #  @post calling object remains unchanged
     */
    Edge operator*() const{
      std::pair<size_type,Proxy_Edge>temp = g->adj[uid][curr];
      return Edge(this->g, uid, std::get<0>(temp));

    }

    /** Return an IncidentIterator
     *  @return an IncidentIterator from incrementation

     */
    IncidentIterator& operator++(){
      this->curr = curr +1;
      return *this;

    }

    /** Return a boolean from equating IncidentIterator
     *  @param[in] iit is an IncidentIterator
     *  @return a boolean from equating IncidentIterators

     *  @post @iit must remain unchanged
     */
    bool operator==(const IncidentIterator& iit) const{

      return curr == iit.curr;

    }


   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    //pointer back to graph
    Graph* g;
    //keep track of current node in adj list
    //keep track of node uid
    size_type uid;
    size_type curr;
    IncidentIterator(const Graph* eg,size_type u, size_type c)
      : g(const_cast<Graph*>(eg)), uid(u), curr(c){

    }


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

    /** Return an Edge
     *  @return an Edge

     *  @post calling object can not change
     */
     Edge operator*() const{
       std::pair<size_type,size_type> temp = g->e[uid];
       return Edge(g,std::get<0>(temp),std::get<1>(temp));
     }

     /** Return an EdgeIterator
      *  @return an EdgeIterator from incrementation

      */
     EdgeIterator& operator++(){
        uid+=1;
        return *this;

     }
      /** Return an boolean from equating EdgeIterators
       *  @return a boolean from equating EdgeIterators

      */
     bool operator==(const EdgeIterator& et) const{

        return this->uid==et.uid;
     }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* g;

    size_type uid;//keep track of current edge in graphs var e(pairs of edges)

    EdgeIterator(const Graph* gg, size_type u)
      : g(const_cast<Graph*>(gg)), uid(u){

    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   /** Return an EdgeIterator at the beginning of the Graph
    *  @return an EdgeIterator at the begining of the Graph

    *  @post calling object must remain unchanged
    */
   EdgeIterator edge_begin() const{

     return EdgeIterator(this,0);
   }

   /** Return an EdgeIterator at the end of a Graph
    * @return an Edgeiterator at the end of a Graph

    *  @post calling object must remain unchanged
    */
   EdgeIterator edge_end() const{

     return EdgeIterator(this,e.size());
   }


    //Node/Edge Removals 
    /**Return the new number of nodes after 1 node deletion
     *@param[in] n Node to be removed.
     *@return The new number of nodes after 1 node deletion

     *@pre @a n must be a valid node.
     *@post @a n must have been remain unchanged in the function.
     *@post All iterators,pointers, etc pointing at @a n before removal
            will be pointing at the end node that took the place of @a n.
     */
    size_type remove_node(const Node& n){
      assert(n.valid());


      //remove duplicate edges incident to node a-O(d)
      for(auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){

        auto e1 = *ei;
        auto n2 = e1.node2();
        v[n2.index()].adj_uid.erase(n.uid);//hash/unordered_map is O(1) to erase 1 element
      }


      //remove node from unique edges O(num_nodes)
      size_type i=0;

      for(;i < e.size();++i){
        std::pair<size_type,size_type> temp = e[i];
                  std::cout<<"hey"<<std::endl;
        if(std::get<0>(temp)==n.uid || std::get<1>(temp)==n.uid){
          //switch this pair with the end pair
          e[i] = e[e.size()-1];
          e.pop_back();

        }

        i = 0;
      }



       auto end = *--v.end(); //this is a proxy node
      //swap adj lists in the vector adj
      adj[n.index()].swap(adj[end.pn_uid]); //O(1)
      //remove back end of adj vector
      adj.pop_back(); //remove node n's adj edges
  
      size_type temp = v[n.index()].pn_uid; // switch the idxs/uids
      //switch places with end element and n element in node vector
      v[n.index()] = end;
      v[n.index()].pn_uid=temp;
      //remove back end element
      v.pop_back();
   

      return v.size();
      
    }

   /**Return a NodeIterator pointing to node following the current node removal.
     *@param[in] n_it NodeIterator pointing to node to be removed.
     *@return a NodeIterator pointing to node following the current node removal.

     *@pre @a n_it must be a valid NodeIterator, pointing to a valid node.
     *@post All other iterators,pointers, etc pointing at @a n before removal
            will be pointing at the end node that took the place of @a n.
     */
    NodeIterator remove_node(NodeIterator n_it){
      
      int new_size = remove_node(*n_it);
      return ++n_it; //return next element after the removed node

    }

   /**Return new number of edges after 1 edge removal.
     *@param[in] a First node on the edge.
     *@param[in] b Second node on the edge. 
     *@return new number of edges after 1 edge removal.

     *@pre @a a and @b must be valid nodes.
     *@post All other iterators,pointers, etc pointing at @a a and @a b before removal will be invalidated. 
     */
    size_type remove_edge(const Node& a, const Node& b){
      assert(a.valid());
      assert(b.valid());
      if(e.size()==0) return 0;

  
      std::pair<size_type,size_type> find1(a.uid,b.uid);
      std::pair<size_type,size_type> find2(b.uid,a.uid);
      //remove from unique edges-O(num_edges)
      for(size_type i =0; i < e.size(); ++i){
         if(e[i]==find1 || e[i]==find2){
           e[i] = e[e.size()-1];//switch with the end element
           e.pop_back();
           break;
         }
      }
      //remove duplicate edges from the node's adj list
      if(!v[a.uid].adj_uid.empty()) v[a.uid].adj_uid.erase(b.uid);    
      if(!v[b.uid].adj_uid.empty()) v[b.uid].adj_uid.erase(a.uid);    
    


      return e.size();//new unique # of edges

    }
   /**Return new number of edges after 1 edge removal.
     *@param[in] e Edge to be removed.
     *@return new number of edges after 1 edge removal.

     *@pre @a e must be a valid edge (it's nodes must be valid nodes).
     *@post All other iterators,pointers, etc pointing to @a e before removal will be invalidated. 
     */
    size_type remove_edge(const Edge& e){

      return remove_edge(e.node1(),e.node2());

    }

   /**Return an EdgeIterator pointing to edge followed by the edge removed.
     *@param[in] e_it EdgeIterator pointing to edge to be removed.
     *@return an EdgeIterator pointing to edge followed by the edge removed.

     *@pre @a e_it must point to a valid edge (it's nodes must be valid nodes).
     *@post All other iterators,pointers, etc pointing to @a e before removal will be invalidated. 
     */
     EdgeIterator remove_edge(EdgeIterator e_it){

         int new_edge_size = remove_edge(*e_it);
         return ++e_it;
     }
     

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  //Internal tpye for graph elements
  struct Proxy_Node {
    Point p; //node position
    node_value_type value;
    size_type deg=0;
    size_type pn_uid;
    std::unordered_map<size_type,size_type> adj_uid;  //idx to track in adj list of a node,key=n2 and val=idx in n2's adj list
  
  };
  struct Proxy_Edge{

   edge_value_type value;
   size_type n2_uid;


  };

  //vector to hold proxy element for Nodes,indexed by Node uid
  std::vector<Proxy_Node> v;
 
  //vector to hold Edges(unique) for edgeiterator
  std::vector<std::pair<size_type,size_type>> e;
  //adjacency list for incidentiterator
  //use node IDs to get their adjnodes, with those ids is first in pair
  std::vector<std::vector<std::pair<size_type,Proxy_Edge>>> adj;

};

#endif // CME212_GRAPH_HPP
