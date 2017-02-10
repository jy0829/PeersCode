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
class Graph : private totally_ordered<Graph<V>>{
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)


  //predeclare struct for nodes and edges
  struct Proxy_Node;


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
  using node_value_type = V;

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
    const Point& position() const {

      //return Point in the uid position of Graphs vector
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

        return G->v[uid].value;

    }
    /** Return the value of a constant Point
     *   @return the value of @a Point

         @pre @a Point calling this func must be const
     */
    const node_value_type& value() const{

      return G->v[uid].value;

    }


    //hw1 part 3
    /** Return the degree of a Node
     *  @return the number of edges from a @a Node

     *  @post calling object must remain unchanged
     */
    size_type degree() const{

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
/*
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
    return Node(this, current-1);        // valid node



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

      // invalid Edge

      // HW0: YOUR CODE HERE

    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE

      return Node(e_G,n1_uid);      // valid Node


   }
    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE

      return Node(e_G,n2_uid);      // valid Node



    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      //check the incoming edge nodes with its own uid
      if((e_G==e.e_G) &&
         (std::min(n1_uid,n2_uid)==std::min(e.n1_uid,e.n2_uid)) &&
         (std::max(n1_uid,n2_uid)==std::max(n1_uid,n2_uid))){
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
        if((std::min(n1_uid,n2_uid) < std::min(e.n1_uid,e.n2_uid)) &&
           (std::max(n1_uid,n2_uid) < std::max(e.n2_uid,e.n1_uid))){
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
    size_type n2_uid;
    //Private constructor for Graph to make valid Edge
    Edge(const Graph* e_g, size_type n1_u, size_type n2_u)
      : e_G(const_cast<Graph*>(e_g)), n1_uid(n1_u), n2_uid(n2_u){

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
      this->adj[a.uid].push_back(b.uid);
      this->v[a.uid].deg +=1; //add to the degree
      this->adj[b.uid].push_back(a.uid);
      this->v[b.uid].deg+=1;
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
      return *this; //THIS COULD BREAK

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
      return Edge(this->g, uid, g->adj[uid][curr]);

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

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  //Internal tpye for graph elements
  struct Proxy_Node {
    Point p; //node position

    node_value_type value;

    size_type deg=0;

  };

  //vector to hold proxy elements for Point
  std::vector<Proxy_Node> v;
  //vector to hold Edges(unique) for edgeiterator
  std::vector<std::pair<size_type,size_type>> e;
  //adjacency list for incidentiterator
  //use node IDs to get their vec. of pairs like (ni,edge id)
  std::vector<std::vector<size_type>> adj;

};

#endif // CME212_GRAPH_HPP
