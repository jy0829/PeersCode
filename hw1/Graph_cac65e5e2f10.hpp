#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP
//#define NDEBUG
/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include <unordered_map>
#include <cassert>
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V>
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
//  const std::vector<Point>* PointsDataPointer;

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
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
    // HW0: YOUR CODE HERE
    //PointsDataPointer = &PointsData;
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
      nodeIndex = -1;
      nodePoints = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return nodePoints->PointsData[nodeIndex];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return nodeIndex;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
   /** Return a reference of the "value" of the node.
    */
    node_value_type& value() {
      return const_cast<Graph*>(nodePoints)->allValues[nodeIndex];
    }
   /** Return a reference of the 'const value' of the node.
    */
    const node_value_type& value() const {
      return nodePoints->allValues[nodeIndex];
    }
   /** Return the degree of a given node.*/
    size_type degree() const {
      return IncidentMap[nodeIndex].size();
    }
    incident_iterator edge_begin() const {
      return IncidentIterator(nodeIndex,0,this->nodePoints);
    }
    incident_iterator edge_end() const {
      return IncidentIterator(nodeIndex,const_cast<Graph*>(this->nodePoints)->IncidentMap[nodeIndex].size(), this->nodePoints);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((nodePoints == n.nodePoints)&&(nodeIndex == n.nodeIndex)) {
        return true;
      }
      (void) n;          // Quiet compiler warning
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
      if (nodeIndex < n.nodeIndex) {
        return true;
      }
      (void) n;           // Quiet compiler warning
      return false;
   // Comment: here I didn't check if the two nodes are in the same graph.
    }

   private:
    // Allow Graph to access Node's private member data and functions.
   // Point
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph::size_type nodeIndex;
    const Graph* nodePoints;
    /** A valid constructor that a graph can use.*/
    Node(const Graph::size_type _nodeIndex,const Graph* _nodePoints) {
      nodeIndex = _nodeIndex;
      nodePoints = _nodePoints;
//      std::cout << "new node index: "<<nodeIndex<<"; nodePoints.size(): "<<(nodePoints->PointsData.size())<<std::endl;
//      std::cout << "nodePoints:"<<nodePoints<<std::endl;
      assert(nodePoints!=nullptr);
      assert((nodeIndex>=0)&&(nodeIndex<(nodePoints->PointsData.size())));
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return PointsData.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return PointsData.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value_ = node_value_type()) {
    // HW0: YOUR CODE HERE
    PointsData.push_back(position);
    allValues.push_back(node_value_);
    Node newNode = Node(PointsData.size()-1,this);
 //   Node newNode = Node();
    // Put info into the new node.
 //   newNode.nodeIndex = allNodesVector.size();
 //   newNode.nodePoints = PointsDataPointer;
 //   newNode.node_value = node_value_;
  //  std::cout<<sizeof(Node)<<std::endl;
    (void) position;      // Quiet compiler warning
    return newNode;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if ((n.nodePoints == this) && (n.nodeIndex < PointsData.size())) {
      return true;
    }
    (void) n;            // Quiet compiler warning
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
    if (i < PointsData.size()) {
      return Node(i,this); 
    }
    (void) i;             // Quiet compiler warning
   // std::cout << "Index out of bound!" << std::endl;
    return Node();        // Invalid node
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
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      n1 = -1;
      n2 = -1; 
      graphpointer = nullptr;     
   //   std::cout<<sizeof(Edge)<<std::endl;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graphpointer->node(n1);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graphpointer->node(n2);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (((n1==e.node1().index())&&(n2==e.node2().index()))||((n1==e.node2().index())&&(n2==e.node1().index()))) {
        return true;
      }
      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type smaller_node = (n1 < n2)?n1:n2;
      size_type larger_node = (n1 < n2)?n2:n1;
      if ((smaller_node < e.node1().index())&&(smaller_node < e.node2().index())) {
        return true;
      } else if (smaller_node == e.node1().index()) {
        return (larger_node < e.node2().index());
      } else if (smaller_node == e.node2().index()) {
        return (larger_node < e.node1().index());
      } else {
        (void) e;           // Quiet compiler warning
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type n1;
    size_type n2;
    const Graph* graphpointer;
    /** Construct a valid undirected edge. A graph can use this.
    * @param _n1: the first node of the edge
    * @param _n2: the second node of the edge.
    * @param _graphpointer: a pointer to the graph to which the edge belongs.
    * @pre @a _graphpointer!=nullptr
    * @pre @a _n1 and @a _n2 are not less than 0, less than @a _graphpointer->PointsData.size(), and not equal.
    * @post @a _n1 == @a n1, @a _n2==@a n2, @a graphpointer == @a _graphpointer,and the preconditions.*/
    Edge(const size_type _n1, const size_type _n2, const Graph* _graphpointer) {
      n1 = _n1;
      n2 = _n2;
      graphpointer = _graphpointer;
      assert(graphpointer!=nullptr);
      assert(n1!=n2);
      assert((n1>=0)&&(n1 < graphpointer->PointsData.size())&&(n2>=0)&&(n2 < graphpointer->PointsData.size()));
    }
    /** Another valid constructor that returns the kth edge in the graph. */
    Edge(const size_type k, const Graph* _graphpointer) { 
      assert(k<(_graphpointer->allEdges1.size()));
      n1 = _graphpointer->allEdges1[k];
      n2 = _graphpointer->allEdges2[k];
      graphpointer = _graphpointer;
    //  std::cout<<"Edge made!"<<n1<<" , "<<n2<<std::endl;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return allEdges1.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
//    if (i < allEdges.size()) {
//      std::set<Edge>::iterator it = allEdges.begin();
//      std::advance(it,i);
//      return *it;
//    }
    if (i < allEdges1.size()) {
      Edge newedge = Edge(allEdges1[i],allEdges2[i],this);
      return newedge;
    }
    (void) i;             // Quiet compiler warning
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
  //  for (int i = 0; i < allNodes1.size(); i++) {
  //    if (((allNodes1[i] == a)&&(allNodes2[i] == b))||((allNodes1[i] == b)&&(allNodes2[i] == a))) {
  //      return true;
  //    }
  //  }
    size_type A = a.index();
    size_type B = b.index();
    for (size_type i=0; i<allEdges1.size(); i++) {
      if (((allEdges1[i]==A)&&(allEdges2[i]==B))||((allEdges1[i]==B)&&(allEdges2[i]==A))) {
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
//    Edge newEdge = Edge();
//    newEdge.n1 = a.index();
//    newEdge.n2 = b.index();
//    newEdge.graphpointer = this;
//    if (!has_edge(a,b)) {
//      allEdges.insert(newEdge);
 //     std::cout<<"Inserted an edge!"<< std::endl;
  //}
//    return newEdge;
 //   (void) a, (void) b;   // Quiet compiler warning
 //   return Edge();        // Invalid Edge
    assert(allEdges1.size()==allEdges2.size());
    size_type A = a.index();
    size_type B = b.index();
    for (size_type i=0; i<allEdges1.size(); i++) {
      if (((allEdges1[i]==A)&&(allEdges2[i]==B))||((allEdges1[i]==B)&&(allEdges2[i]==A))) {
        return Edge(allEdges1[i],allEdges2[i],this);
      }
    }
    allEdges1.push_back(A);
    allEdges2.push_back(B);
//    newEdge.n1 = A;
//    newEdge.n2 = B;
//    newEdge.graphpointer = this;
    if (IncidentMap.find(A)!=IncidentMap.end()) {
      IncidentMap[A].push_back(B);
    } else {
      IncidentMap[A] = std::vector<size_type>(1,B);
    }
    if (IncidentMap.find(B)!=IncidentMap.end()) {
      IncidentMap[B].push_back(A);
    } else {
      IncidentMap[B] = std::vector<size_type>(1,A);
    }

    return Edge(A,B,this);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    PointsData.clear();
  //  allNodesVector.clear();
  //  allEdges.clear();
    allEdges1.clear();
    allEdges2.clear();
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
    /** Return the node_index.
    */
    size_type get_index() const {
      return node_index;
    }
    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      node_index = -1;
      graph_pointer = nullptr;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference an iterator to the node.
    */
    Node operator*() const {
      return graph_pointer->node(node_index);
    }
    /** Increment the iterator to next Node.
    */
    NodeIterator& operator++() {
      ++(this->node_index);
      return *this;
    }
    /** Compare if two node iterators are equal
    */
    bool operator==(const NodeIterator& node_iter_to_compare) const {
      return (this->graph_pointer == node_iter_to_compare.graph_pointer)&&((this->node_index) == (node_iter_to_compare.get_index()));
    }
    bool operator<(const NodeIterator& node_iter_to_compare) const {
      return (this->graph_pointer == node_iter_to_compare.graph_pointer)&&((this->node_index) < (node_iter_to_compare.get_index()));
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type node_index;
    const Graph* graph_pointer;
    NodeIterator(size_type _node_index, const Graph* _graph_pointer) {
      node_index = _node_index;
      graph_pointer = _graph_pointer;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return an iterator that points to the first node in the graph.
  */
  node_iterator node_begin() const {
 //   node_iterator to_return = NodeIterator();
 //   to_return.graph_pointer = this;
 //   to_return.node_index = 0;
    return NodeIterator(0,this);
  }
  /** Return an iterator that points to the end (i.e. just after the
  * last) node of the graph.
  */ 
  node_iterator node_end() const {
//    node_iterator to_return = NodeIterator();
//    to_return.graph_pointer = this;
//    to_return.node_index = (this->allNodesVector).size();
    return NodeIterator((this->PointsData).size(),this);
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      spawn = -1;
      incident_index = -1;
      graphpointer = nullptr;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereference the IncidentIterator to get the edge.*/ 
    Edge operator*() const {
      assert((spawn>=0)&&(spawn<(graphpointer->PointsData.size())));
      assert((incident_index>=0)&&(incident_index<=(const_cast<Graph*>(graphpointer)->IncidentMap[spawn].size())));
      return Edge(spawn, const_cast<Graph*>(graphpointer)->IncidentMap[spawn][incident_index],graphpointer);
    }
    /** Increment the incident operator. */
    IncidentIterator& operator++() {
      ++(this->incident_index);
      return *this;
    }
    /** Compare if two incident operators are equal.*/
    bool operator==(const IncidentIterator& it) const {
      return (this->spawn == it.spawn)&&(this->incident_index == it.incident_index)&&(this->graphpointer == it.graphpointer);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type spawn;  //The node that spawns the iterator
    size_type incident_index;
    const Graph* graphpointer;
    IncidentIterator(size_type _spawn,size_type _incident_index, const Graph* _pointer) {
      spawn = _spawn;
      incident_index = _incident_index;
      graphpointer = _pointer;
      assert(graphpointer!=nullptr);
      assert((spawn>=0)&&(spawn<(graphpointer->PointsData.size())));
      assert((incident_index>=0)&&(incident_index<=(const_cast<Graph*>(graphpointer)->IncidentMap[spawn].size())));
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
      graphpointer_ = nullptr;
      edge_index_ = -1;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Return the edge an edge iterator points to.
    Edge operator*() const {
      assert(graphpointer_!=nullptr);
      assert((edge_index_>=0)&&(edge_index_<(graphpointer_->allEdges1.size())));
   //   std::cout<<"Dereferencing an edge iterator: "<<edge_index_<<" in "<<graphpointer_->allEdges1.size()<<std::endl;
      return Edge(edge_index_,graphpointer_);
    }
    // Increment the edge iterator.
    EdgeIterator& operator++() {
      edge_index_ = edge_index_+1;
   //   std::cout<<"Incremented!!"<<edge_index_<<std::endl;
      return *this;
    }
    // Compare two edge iterators.
    bool operator==(const EdgeIterator& other) const {
      return (edge_index_==other.edge_index_)&&(graphpointer_==other.graphpointer_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graphpointer_;
    size_type edge_index_;
    // A valid edge iterator.
    EdgeIterator(const Graph* graphpointer,size_type edge_index):graphpointer_(graphpointer),edge_index_(edge_index) {
      assert(graphpointer_!=nullptr);
      assert((edge_index_>=0)&&(edge_index_<=(graphpointer_->allEdges1.size())));
   // std::cout<<"A new edge iterator:"<<edge_index_<<" in "<<graphpointer_->allEdges1.size()<<std::endl;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // An edge iterator that points to the first edge.
  edge_iterator edge_begin() const {
    return EdgeIterator(this,0);
  }
  // An edge iterator that points just after the last edge.
  edge_iterator edge_end() const {
    return EdgeIterator(this,allEdges1.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> PointsData;
  std::vector<node_value_type> allValues;
//  std::vector<Node> allNodesVector;
//  std::set<Edge> allEdges;
  std::vector<size_type> allEdges1;
  std::vector<size_type> allEdges2;
  std::unordered_map<size_type,std::vector<size_type> > IncidentMap;
};

#endif // CME212_GRAPH_HPP
