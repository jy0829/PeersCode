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

 public:

  using node_value_type = V;

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

  using adj_type = std::vector<std::vector<size_type>>;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Declaring private attributes points, edge class and vector

 std::vector<std::pair<Point,node_value_type>> Points;


 struct edge_items
 {
     size_type NodeId1;
     size_type NodeId2;
 
     edge_items(size_type id1,size_type id2)

     {
          NodeId1=id1;
          NodeId2=id2;
     }
 };

 std::vector<std::vector<size_type>> AdjList;

 std::vector<edge_items> Edges;

 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //


  /** Construct an empty graph. */

  Graph() : Points(),Edges() {}
 

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
    }


    /** Returns the number of adjacent edges to this node
    * @return @a degree, which is the number of adjacent edges of a valid node
    *
    * Complexity: No more than O(1)
    */

      size_type degree() const
    {
       return GraphPointer->AdjList[NodeId].size();
    }

    /** Sets the incident iterator to the beginning of adjacent edges
    * @post @a *ii corresponding first edge of the current graph
    * @return @a ii, a valid incident iterator object of the current graph
    *
    * Complexity: No more than O(1)
    */

    incident_iterator edge_begin() const
    {
       IncidentIterator IncIterObject(GraphPointer,NodeId,0);
       return IncIterObject;
    }

    /** Sets the incident iterator to the end
    * @post @a iiend is one past the last incident iterator object, @a ii, in which @ *ii was valid operation
    * @return @a iiend, a valid incident iterator object of the current graph
    *
    * Complexity: No more than O(1)
    */

    incident_iterator edge_end() const
    {
        IncidentIterator IncIterObject(GraphPointer,NodeId,degree());
        return IncIterObject;
    }


    /** Return this node's accessible value. */
    node_value_type& value()
    {
        return GraphPointer->Points[NodeId].second;
    }

    /** Return this node's non-accesible value. */

    const node_value_type& value() const
    {
        return GraphPointer->Points[NodeId].second;  
    }
    


    /** Return this node's position. */
    const Point& position() const {
        return GraphPointer->Points[NodeId].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(GraphPointer->Points.size()>NodeId);
      return NodeId;
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
      return GraphPointer == n.GraphPointer && index() == n.index();
      }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const 
    {

        if (GraphPointer == n.GraphPointer) 
            return NodeId < n.NodeId;
        else 
            return GraphPointer < n.GraphPointer;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Constructor that initalizes Node object

    Node(const Graph* currentgraph, size_type index)
    {
        NodeId = index;
        GraphPointer = const_cast<Graph*>(currentgraph);
       // GraphPointer->Points[index].second = nodevalue;
    }

    // Declaring private variables

    size_type NodeId;
    Graph* GraphPointer;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const 
  {
      return Points.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const 
  {
      return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) 
  {

      // Pushing back position to points vector
      std::pair<Point,node_value_type> nodedata(position,node_value);
     
      Points.push_back(nodedata); // Add x,y,z and value to point vector
      std::vector<size_type> v;
      AdjList.push_back(v);
      Node NodeObject(this, size()-1);
      return NodeObject;

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */

  bool has_node(const Node& n) const 
  {
      return (n.GraphPointer==this && n.index() < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */

  Node node(size_type i) const 
  {
      assert(Points.size()>i); // Asserting numnodes > i
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
    }

    Edge(const Graph* currentgraph, size_type id1, size_type id2)
    {
        NodeId1=id1;
        NodeId2=id2;
        GraphPointer=const_cast<Graph*>(currentgraph);
    }

    /** Return a node of this Edge */
    Node node1() const 
    {
      return Node(GraphPointer,NodeId1);      
    }

    /** Return the other node of this Edge */
    Node node2() const 
    {
      return Node(GraphPointer,NodeId2);   
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const 
    {
   
    // Checking for equality
 
    return (GraphPointer==e.GraphPointer &&
            std::min(NodeId1,NodeId2)==std::min(e.NodeId1,e.NodeId2) &&
            std::max(NodeId1,NodeId2)==std::max(e.NodeId1,e.NodeId2));
     
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const 
    {

      return (GraphPointer==e.GraphPointer && 
              std::min(NodeId1,NodeId2)<=std::min(e.NodeId1,e.NodeId2) &&
              std::max(NodeId1,NodeId2)<=std::max(e.NodeId1,e.NodeId2));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* GraphPointer;
    size_type NodeId1;
    size_type NodeId2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  /** Fills out the adjacency list
  * @pre @a NodeId1 and @a NodeId2 correspond to valid nodes of the graph
  * @pre @a NodeId1 != @a NodeId2
  * @param[in] @a AdjList, which is private data variable of graph class
  * @param[in] @a NodeId1,NodeId2 which are nodes connected by one edge of graph
  * @post @a AdjList[NodeId1].size() = 1 + oldsize
  * @post @a AdjList[NodeId2].size() = 1 + oldsize

  */
  void Adjacency(adj_type& AdjList ,size_type NodeId1, size_type NodeId2)
  {
      assert(NodeId1 < num_nodes() && NodeId2 < num_nodes());
      assert(NodeId1 != NodeId2);

      AdjList[NodeId1].push_back(NodeId2);
      AdjList[NodeId2].push_back(NodeId1);

  }

  size_type num_edges() const 
  {
    return Edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const 
  {

    assert(Edges.size()>i); //Asseting that i < number of edges
    Edge EdgeObject(this,Edges[i].NodeId1,Edges[i].NodeId2);
    return EdgeObject;    

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const 
  {

    assert(a.GraphPointer== this and b.GraphPointer == this); // asserting nodes in graph
    assert(a.NodeId < size() && b.NodeId < size()); // asserting nodes are valid

    std::vector<size_type> connected_nodes = AdjList[a.NodeId];

    for (size_type ind = 0; ind < connected_nodes.size(); ind++)
    {
        if (connected_nodes[ind] == b.NodeId)
        {
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

  Edge add_edge(const Node& a, const Node& b) 
  {
    assert(this==a.GraphPointer && this == b.GraphPointer && // asserting preconditions
    a.NodeId < Points.size() && b.NodeId < Points.size() && a.NodeId!=b.NodeId);
    
    if (has_edge(a,b))
    {
        Edge EdgeObject(this,a.NodeId,b.NodeId);
        return EdgeObject;
    }
      assert(a.GraphPointer != nullptr && b.GraphPointer != nullptr);
      assert(a.NodeId != b.NodeId);



      edge_items edgeData(a.NodeId,b.NodeId);
    Edges.push_back(edgeData);
    Adjacency(AdjList,a.NodeId,b.NodeId); // Adding adjacency list
    Edge EdgeObject(this,a.NodeId,b.NodeId);
    return EdgeObject;
    

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Points.clear();
    Edges.clear();
    AdjList.clear();
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

    // Custom constructor
  
    NodeIterator(const Graph* currentgraph, size_type id)
    {
         GraphPointer = const_cast<Graph*>(currentgraph);
         NodeId = id;
    }
    
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Dereferences the node iterator
    * @pre @a ni is a valid iterator of this graph and != to node_end()
    * @return A Node iterator object @a ni with @a *ni == corresponding node
    * @post @a *ni == corresponding node of current graph
    *
    * Complexity: No more than O(1) complexity
    */

     Node operator*() const
     {
         assert(NodeId < GraphPointer->num_nodes());
         return Node(GraphPointer, NodeId);
     }

    /** Forwards the node iterator
    * @pre @a ni is a valid iterator of this graph and != to node_end().
    * @pre @a ni is dereferencible
    * @return A Node iterator object @a ninext
    * @post @a *ni == node of current graph @a ni != @a ninext
    *
    * Complexity: No more than O(1) complexity
    */

     node_iterator& operator++() 
     {
        assert(NodeId < GraphPointer->num_nodes());
        NodeId++;
	    return *this;
     }

      /** Testing whether node iterators are equal
     * @pre @a ni and @a nit are valid iterators of the current graph
     * @param[in] @a nit, which is a node iterator object
     * @return A boolean object which returns whether @a ni == @a nit
     *
     * Complexity: No more than O(1)
     */

     bool operator==(const node_iterator& nit) const 
     {
        assert(NodeId <= GraphPointer->num_nodes() &&
               nit.NodeId <= nit.GraphPointer->num_nodes());

        return (GraphPointer == nit.GraphPointer) && (NodeId == nit.NodeId);
     }


    private: 


    friend class Graph;

    size_type NodeId;
    const Graph* GraphPointer;
  };

    /** Sets the node iterator to the beginning
    * @post @a *ni corresponding node of the current graph
    * @return @a ni, a valid node iterator object of the current graph
    *
    * Complexity: No more than O(1)
    */

    node_iterator node_begin() const
     {
         assert(has_node(Node(this,0)));
         return NodeIterator(this, 0);
     }

    /** Sets the node iterator to the end
    * @post @a niend is one past the last node iterator object, @a ni, in which @ ni
    * is dereferencible
    * @return @a niend, a valid node iterator object of the current graph
    *
    * Complexity: No more than O(1)
    */

     node_iterator node_end() const
     {
         assert(has_node(Node(this,size()-1)));
         return NodeIterator(this,size());
     }

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

    // Custom Constructor

    IncidentIterator(const Graph* currentgraph, size_type nodeid, size_type id2pos)
    {
        GraphPointer = const_cast<Graph*>(currentgraph);
        NodeId1 = nodeid;
        NodeId2Idx = id2pos;
    }

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Dereferences the incident iterator
    * @pre @a ii is a valid incident iterator of this graph and != to edge_end()
    * @return An incident iterator object @a ii with @a *ii == corresponding adjacent edge
    * @post @a *ii == corresponding adjacent edge of current graph
    *
    * Complexity: No more than O(1) complexity
    */

     Edge operator*() const
     {
         assert(NodeId1 != GraphPointer->AdjList[NodeId1][NodeId2Idx]
         && NodeId1 < GraphPointer->num_nodes());
         assert(GraphPointer->AdjList[NodeId1][NodeId2Idx] <= GraphPointer->num_nodes());
         return Edge(GraphPointer,NodeId1,GraphPointer->AdjList[NodeId1][NodeId2Idx]);
     }

    /** Forwards the incident iterator
    * @pre @a ii is a valid iterator of this graph and != to edge_end()
    * @return An incident iterator object @a iinext
    * @post @a *ii == valid adjacent edge of the specific node of
    * current graph @a ii != @a iinext
    *
    * Complexity: No more than O(1) complexity
    */

     incident_iterator& operator++()
     {
         assert(GraphPointer->AdjList[NodeId1][NodeId2Idx]<=GraphPointer->num_nodes());
         NodeId2Idx++;
         return *this;
     }

     /** Testing whether incident iterators are equal
     * @pre @a ii and @a iit are valid incident iterators of the current graph
     * @param[in] @a iit, which is a valid incident iterator object
     * @return A boolean object which returns whether @a ii == @a iit
     *
     * Complexity: No more than O(1)
     */

     bool operator == (const incident_iterator& iit) const
     {
          assert(NodeId1 <= GraphPointer->num_nodes());
          assert(iit.NodeId1 <= iit.GraphPointer->num_nodes());

          return (GraphPointer==iit.GraphPointer && NodeId1==iit.NodeId1 && NodeId2Idx == iit.NodeId2Idx);
     }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
     size_type NodeId2Idx;
     size_type NodeId1;
     const Graph* GraphPointer;

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

    EdgeIterator(const Graph* currentgraph, size_type eid)
    {
        GraphPointer = const_cast<Graph*>(currentgraph);
        EdgeId = eid;
    }

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const


    /** Dereferences the edge iterator
    * @pre @a ei is a valid edge iterator of this graph and != to edge_end()
    * @return An edge iterator object @a ei with @a *ei == corresponding edge
    * @post @a *ei == corresponding edge of current graph
    *
    * Complexity: No more than O(1) complexity
    */

    Edge operator* () const
    {
        assert(EdgeId < GraphPointer->num_edges());
        return Edge(GraphPointer,GraphPointer->Edges[EdgeId].NodeId1,GraphPointer->Edges[EdgeId].NodeId2);
    }

    /** Forwards the edge iterator
    * @pre @a ei is a valid edge iterator of this graph and != to edge_end()
    * @return An edge iterator object @a einext
    * @post @a *ei == valid edge of current graph @a ei != @a einext
    *
    * Complexity: No more than O(1) complexity
    */

    edge_iterator& operator++()
    {
        assert(EdgeId < GraphPointer->num_edges());
        EdgeId++;
        return *this;
    }

    /** Testing whether edge iterators are equal
    * @pre @a ei and @a eit are valid edge iterators of the current graph
    * @param[in] @a eit, which is a valid edge iterator object
    * @return A boolean object which returns whether @a ei == @a eit
    *
    * Complexity: No more than O(1)
    */

    bool operator == (const edge_iterator& eit) const
    {
        assert(EdgeId <= GraphPointer->num_edges() &&
               eit.EdgeId <= eit.GraphPointer->num_edges());
        return (GraphPointer == eit.GraphPointer && EdgeId == eit.EdgeId);
    }
   
   private:
    friend class Graph;

    const Graph* GraphPointer;
    size_type EdgeId;

  };


  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

   /** Sets the edge iterator to the beginning
   * @post @a *ei corresponding edge of the current graph
   * @return @a ni, a valid edge iterator object of the current graph
   *
   * Complexity: No more than O(1)
   */

   edge_iterator edge_begin() const
   {
       return EdgeIterator(this, 0);
   }

   /** Sets the edge iterator to the end
   * @post @a eiend is one past the last edge iterator object, @a ei, in which @ *ei was valid
   * @return @a eiend, a valid edge iterator object of the current graph
   *
   * Complexity: No more than O(1)
   */

   edge_iterator edge_end() const
   { 
       return EdgeIterator(this,Edges.size());
   }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

 public:

    void test_function() {
        int count = 0;
        for (auto ei = edge_begin(); ei != edge_end(); ++ei) {

            std::cout << "This is index of node1 " << (*ei).node1().index() << " of edge " << count << std::endl;
            std::cout << "This is index of node2 " << (*ei).node2().index() << " of edge " << count << std::endl;
            count++;
        }

        for (auto ni = node_begin(); ni != node_end(); ++ni) {
            std::cout << "This is the index of node " << (*ni).index() << std::endl;
            std::cout << "The size of the vector within AdjList is " << AdjList[(*ni).index()].size() << std::endl;
            std::cout << "First element is" << AdjList[(*ni).index()][0] << std::endl;
            std::cout << "First element is" << AdjList[(*ni).index()][1] << std::endl;
            std::cout << "First element is" << AdjList[(*ni).index()][2] << std::endl;

        }

        std::cout << "The number of nodes are " << num_nodes() << std::endl;
        std::cout << "The number of edges are " << num_edges() << std::endl;


        std::cout << "This is the size of the AdjList " << Edges.size() << std::endl;

    }
};
#endif // CME212_GRAPH_HPP



