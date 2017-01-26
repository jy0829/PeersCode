#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    nid = 0;
    eid=0;
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
  class Node {
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
      nid_ = -1;	//invalid value for unsigned
      e_=0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g_->node_element[nid_];
      //return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nid_;
      //return size_type(-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      if((this->g_ == n.g_) && (this->nid_==n.nid_))
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
      (void) n;           // Quiet compiler warning
      if(this->g_<n.g_)
        return true;
      else if((this->g_==n.g_) && (this->nid_<n.nid_))
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type nid_;
    Graph* g_;
    std::vector<node_type> nodes;
    size_type e_; //index for vector
    Node(const Graph* g, size_type nid): nid_(nid),g_(const_cast<Graph*>(g)){
    }

    //vector<node_type> edges;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_element.size();
    return 0;
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    node_element[nid] = position;
    nid++;
    //node_edge.push_back(0);
    return Node(this, nid-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    size_type i = n.index();
    if(i>num_nodes()-1)
        return false;
    return true;
    //return false;
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      eid_ = -1;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return n1_;
      //return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return n2_;
      //return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if((this->g_ == e.g_) && (this->eid_==e.eid_))
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if(this->g_<e.g_)
        return true;
      else if((this->g_==e.g_) && (this->eid_<e.eid_))
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type eid_;
    node_type n1_;
    node_type n2_;
    Graph* g_;
    Edge(const Graph* g, size_type eid, node_type n1, node_type n2): eid_(eid),n1_(n1),n2_(n2),g_(const_cast<Graph*>(g)){
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_element.size();
    //return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    if(edge_element.find(i)!=edge_element.end())
        return Edge(this,i,edge_element.at(i).first,edge_element.at(i).second);

    return Edge();        // Invalid Edge
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
    for(size_type i=0;i<num_edges();i++)
    {
	if((edge_element.at(i).first == a) && (edge_element.at(i).second == b))
	{
	    return true;
        }
        else if((edge_element.at(i).first == b) && (edge_element.at(i).second == a))
	{
	    return true;
        }
    }
    /*std::pair<std::multimap<node_type, std::pair<size_type, node_type>>::iterator, std::multimap<node_type, std::pair<size_type, node_type>>::iterator> range;
    range = node_edge.equal_range(a);
    for(std::multimap<node_type,std::pair<size_type,node_type>>::iterator it=range.first; it!=range.second;++it)
    {
        std::pair<size_type,node_type> p = it->second;
        if(p.second==b)
		return true;
    }
    range = node_edge.equal_range(b);
    for(std::multimap<node_type,std::pair<size_type,node_type>>::iterator it=range.first; it!=range.second;++it)
    {
        std::pair<size_type,node_type> p = it->second;
        if(p.second==a)
		return true;
    }*/
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
    //Search if edge exists
    std::pair<std::multimap<node_type, std::pair<size_type, node_type>>::iterator, std::multimap<node_type, std::pair<size_type, node_type>>::iterator> range;
    range = node_edge.equal_range(a);
    for(std::multimap<node_type,std::pair<size_type,node_type>>::iterator it=range.first; it!=range.second;++it)
    {
        std::pair<size_type,node_type> p = it->second;
        if(p.second==b)
		return Edge(this, p.first,a,b);
    }
    range = node_edge.equal_range(b);
    for(std::multimap<node_type,std::pair<size_type,node_type>>::iterator it=range.first; it!=range.second;++it)
    {
        std::pair<size_type,node_type> p = it->second;
        if(p.second==a)
		return Edge(this, p.first,a,b);
    }
    //if edge not found, add edge
    (void) a, (void) b;   // Quiet compiler warning
    std::pair<node_type,node_type> new_edge(a,b);
    edge_element[eid] = new_edge;
    eid++;

    std::pair<size_type,node_type> p(eid-1,b);
    node_edge.insert(std::pair<node_type,std::pair<size_type,node_type>>(a,p));

    return Edge(this, eid-1,a,b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_element.clear();
    edge_element.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::unordered_map<size_type, Point> node_element;
  std::unordered_map<size_type, std::pair<node_type,node_type>> edge_element;
  std::multimap<node_type,std::pair<size_type,node_type>> node_edge;
  size_type nid;
  size_type eid;
};

#endif // CME212_GRAPH_HPP
