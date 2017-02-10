#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
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
  
/**point_vect is a vector of nodes encoded by their position. We use an adjency matrix to encode 
the neighbors of node of index i at Adj_vect[i] */
std::vector<Point> Point_vect;
std::vector<std::pair<unsigned, unsigned>> Edge_vect;
std::vector<std::vector<unsigned>> Adj_vect;
/**HW1 a vector<V> is used to store the value of node[i] at the place i in the vector*/
std::vector<V> Value_vect;

  unsigned number_edges;


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

  using node_value_type = V;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */

  /**For designing graph, keep track of num_edges forthe siwe_edge function*/
  Graph(): number_edges(0),Point_vect(), Adj_vect(),Value_vect() {
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
  class Node : private totally_ordered<Node>  {
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
    Node():dx(0),pgraph(NULL) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      //return Point();
      return pgraph->Point_vect[dx];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //return size_type(-1);
      return dx;
    }

    /**I add a function degree which return the degree of a node for the sake of the
    function has_edge*/


    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /**node.value() will return the value of the object node.*/
    node_value_type& value(){
      return pgraph->Value_vect[dx];

    }
  /** const node_value_type& value() is the constant version of the node_value_type& value() function. it is
  still returning the value of the node */
    const node_value_type& value() const{
      return pgraph->Value_vect[dx];
    }
  /**degree() returns the degree of the node without modifying it*/
   size_type degree() const{
return pgraph->Adj_vect[dx].size();

   }
/**edge_begin() return an the begining iterator  that will go through all the edges connected to the node*/
     incident_iterator edge_begin() const{
      return IncidentIterator(pgraph,dx,0);
     }
     /**edge_end() return the ending iterator which has been through all edges connected to the node.
     dereferencing this pointer is a bad behaviour (as it means nothing because it is pointing to nothing 
     relevant)*/
    incident_iterator edge_end() const{
      return IncidentIterator(pgraph,dx,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //(void) n;          // Quiet compiler warning
      //return false;
      return (n.dx==dx) && (n.pgraph==pgraph);
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
      //(void) n;           // Quiet compiler warning
      //return false;
      bool a;
      if (pgraph == n.pgraph){a=(dx<n.dx);}
      else{a =(pgraph<n.pgraph);}
        return a;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
        Node(const Graph* _pgraph, size_type idx)
        : dx(idx), pgraph(const_cast<Graph*>(_pgraph)){}
    
    size_type dx; // the node index in Graph
    Graph* pgraph; //point to current Graph



  }; //end of class Node.

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    //return 0;
    return Point_vect.size(); /**here dont need this->Point_vect.size() because there is 
    no ambiguity */

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
  Node add_node(const Point& position, const node_value_type& myval = node_value_type()) {
    // HW0: YOUR CODE HERE
    //(void) position;      // Quiet compiler warning
    //return Node();        // Invalid node


    /**Use the Node interface of the friend class Graph*/
    Point_vect.push_back(position);
    Value_vect.push_back(myval);
    Adj_vect.push_back(std::vector<unsigned>());
    return Node(this, Point_vect.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //(void) n;            // Quiet compiler warning
    //return false;
    return (n.pgraph == this) && (n.dx < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    //return Node();        // Invalid node
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() : dx1(0),dx2(0),pgraph(NULL){
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      //return Node();      // Invalid Node
      return Node(pgraph,std::min(dx1,dx2));
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      //return Node();      // Invalid Node
      return Node(pgraph,std::max(dx1,dx2));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //return false;
      std::set<size_type> myset1;
      std::set<size_type> myset2;
      myset1.insert(dx1);
      myset1.insert(dx2);
      myset2.insert(e.dx1);
      myset2.insert(e.dx2);
      return (pgraph == e.pgraph) && (myset2==myset1);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //return false;
//in add_edge we will have to be careful about sorting indices

      bool a;
      if (pgraph==e.pgraph){a = (dx1<e.dx1);}
      else {a=(pgraph<e.pgraph);} 
      return a;   }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* _pgraph, size_type index_1, size_type index_2)
        : dx1(index_1), dx2(index_2), pgraph(_pgraph){} //pointer has not to be changed here: it is constant
    size_type dx1;
    size_type dx2; // the node indices in Graph
    const Graph* pgraph; 
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return number_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    //return Edge();        // Invalid Edge
    /**if it is a valid edge, return the basic edge with the null pointer*/
    assert(i >= 0 && i < num_edges());
    size_type dx1 =Edge_vect[i].first;
    size_type dx2 =Edge_vect[i].second;
    return Edge(this,dx1,dx2);
    

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //(void) a; (void) b;   // Quiet compiler warning
    //return false;

    //index of b has to be in the adjacency list of a before Adj_vect[a.index()].end()

    if (a.degree() < b.degree()) {
      return (std::find(Adj_vect[a.index()].begin(), Adj_vect[a.index()].end(),
                 b.index()) != Adj_vect[a.index()].end()); }
    else {
      return (std::find(Adj_vect[b.index()].begin(), Adj_vect[b.index()].end(),
                 a.index()) != Adj_vect[b.index()].end()); }
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
    //(void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
    if (!has_edge(a, b)) {
      Adj_vect[a.index()].push_back(b.index());
      Adj_vect[b.index()].push_back(a.index());
      ++number_edges;
  	Edge_vect.push_back(std::make_pair(a.index(),b.index()));}
    return Edge(this, a.index(), b.index());
  }
  

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    Point_vect.clear();
    Adj_vect.clear();
    number_edges=0;
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
    NodeIterator() : pgraph(nullptr),dx(0) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /**operator*() is dereferencing operator
    *@pre  legitimate iterator is not the end node.
    *@return node at the location of the iterator
    *complexity is O(1)
    */
    Node operator*() const{
      return Node(pgraph,dx);
    }
    /** operator++() is incrementing the iterator
    *@pre legitimate iterator is not equal to the end node
    *@return the iterator with the new position (old_position+1)
    *complexity is O(1)
    */
    NodeIterator& operator++(){
      ++dx; 
      return *this;
    }
    /** operator== tests for equlity between two iterators
    *@param[in] iter_ is the iterator to be compared to.
    *@pre the two iterator have to be valid
    *@return the relevant bool
    *complexity O(1)

    */
    bool operator==(const NodeIterator& iter_) const{
      return (pgraph == iter_.pgraph) && (dx==iter_.dx);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    NodeIterator(const Graph* _pgraph, size_type _dx):
    pgraph(const_cast<Graph*>(_pgraph)),dx(_dx){
    }
    const Graph* pgraph;
    size_type dx;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /**node_begin() returns the iterator which points at the beginning of the nodes of the current graph
  */
  node_iterator node_begin() const{
    return node_iterator(this,0);
  }
  /**node_end() returns the iterator "end" of the nodes of the current grpah
  *@pre it should not be dereferenced as it would mean nothing
  */
   node_iterator node_end() const{
    return node_iterator(this, Point_vect.size());
   }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : _pgraph(nullptr),_dx(0),rank_in_Adj(0){
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /**operator*()  is dereferencing for elements of the class IncidentIterator
    *@pre element of IncidentIterator should not be equal to the end iterator
    *@post the result of operator*() say res, has first node given by res.node1()
    *@return element of Edge class incident to the node in question  and which 
    position is the one pointed by the iterator
    */
    Edge operator*() const{
      return Edge(_pgraph,_dx,_pgraph->Adj_vect[_dx][rank_in_Adj]);
    }
    /** operator++() is incrementing the iterator which points to the list of adjacency of the specific node
    *@pre iterator should not be equal to the end iterator
    *@return the iterator incremented by 1

    */
     IncidentIterator& operator++(){
      ++rank_in_Adj;
      return *this;
     }
     /**operator== is comparing two iterators of the class IncidentIterator.
     *@param[in] iter, the iterator we want to compare to our iterator (this)
     *@pre this and iter have to be constructed with the valid method
     *@ return the bool assessing wether or not the two iterators are equals
     */
    bool operator==(const IncidentIterator& iter) const{
      return (_pgraph == iter._pgraph) && (_dx==iter._dx)&& (rank_in_Adj==iter.rank_in_Adj);
    }

   private:
    friend class Graph;
    friend class Node; /**as we will call the iterators at a specific node, the class Node have to inherit 
    of these  properties*/
    // HW1 #3: YOUR CODE HERE
    const Graph* _pgraph;
    size_type _dx;
    size_type rank_in_Adj;
    IncidentIterator(const Graph* pgraph_,size_type dx_,size_type dx2_)
    :_pgraph(const_cast<Graph*>(pgraph_)), _dx(dx_),rank_in_Adj(dx2_) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : _pgraph(nullptr),_dx(0),rank_in_Adj(0){
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /**operator*()  is dereferencing for elements of the class EdgeIterator
    *@pre element of EdgeIterator should not be equal to the end iterator
    *@post the result of operator*() say res, has first node given by res.node1()
    *@return element of Edge class incident to the node in question  and which 
    position is the one pointed by the iterator
    */
     Edge operator*() const{
      return Edge(_pgraph,_dx,_pgraph->Adj_vect[_dx][rank_in_Adj]);
     }
    /** operator++() is incrementing the iterator which points to the list of adjacency of the specific node
    *@pre iterator should not be equal to the end iterator
    *@return the iterator incremented by 1

    */
    EdgeIterator& operator++(){
      ++rank_in_Adj;
      aux();
      return *this;
    }
    /**operator==() tests for equality for elements of
    the EdgIterator class
    *@param[in] iterator to be compared to
    *@param[out] boolean that states for the equality
    */
    bool operator==(const EdgeIterator& iter) const{
    	return (_pgraph==iter._pgraph) && (_dx==iter._dx) && (rank_in_Adj==iter.rank_in_Adj);
    }

   private:
    friend class Graph;
    friend class Edge;

    const Graph* _pgraph;
    size_type _dx;
    size_type rank_in_Adj;
    EdgeIterator(const Graph* pgraph_,size_type dx_,size_type dx2_):
    _pgraph(const_cast<Graph*>(pgraph_)),
    _dx(dx_),rank_in_Adj(dx2_){aux();}

    /**aux() makes sure that not both Edge(a,b) and Edge(b,a) are traversed by the iterator
      by just providing edges of type (a,b=adj_of_a) with value(a)<value(adj_of_a).
    */
      void aux(){
        while(_dx<_pgraph->Adj_vect.size()){
          while (rank_in_Adj<_pgraph->Adj_vect[_dx].size()){
            if(_dx<_pgraph->Adj_vect[_dx][rank_in_Adj]) 
              {return;}
            ++rank_in_Adj;
          }
          ++_dx;
          rank_in_Adj=0;
        }
        return;
      }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /**edge_begin()  return an iterator initiated at the beginning of the edges enumeration*/
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0,0);
  }
    /**edge_begin()  return an iterator initiated at the end of the edges enumeration*/

  edge_iterator edge_end() const{
    return EdgeIterator(this,Adj_vect.size(),0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
