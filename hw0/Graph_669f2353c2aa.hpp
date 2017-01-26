#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <set>

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
  // Store a Vector of Points, indexed by Node's index
  std::vector<Point> NodeLocations; 
  // Graph's storage of Edges is defined at the bottom
  // it is a map of sets. The key is the lower valued index of
  // one of the two nodes that comprise the edge. The set is the larger
  // valued indices that edge connects to

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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->NodeLocations[id]; 
      // HW0: YOUR CODE HERE
      //return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //return size_type(-1);
      return id;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //Still need to check in same graph
      if(id==n.index() and graph_==n.graph_){
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
      if(graph_==n.graph_ and id<n.index()){
         return true; 
      } 
      if(graph_<n.graph_){
         return true; 
      } 
      (void) n;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    //Index of a given node 
    size_type id; 
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Private Node constructor
    Node(const Graph* graph,size_type index) :graph_(const_cast<Graph*>(graph)),id(index){
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return NodeLocations.size();
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
    // HW0: YOUR CODE HERE
  Node add_node(const Point& position){ 
    Node NewNode;
    NewNode.id=NodeLocations.size();
    //Pushback is O(1) amortized 
    NodeLocations.push_back(position);       
    NewNode.graph_=this;
    //(void) position;      // Quiet compiler warning
    //return Node();        // Invalid node
    return NewNode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if(n.index()<this->size())
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
    Node SpecificNode(this,i);
//    (void) i;             // Quiet compiler warning
//    return Node();        // Invalid node
    return SpecificNode;
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return FirstNode;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return SecondNode;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(FirstNode==e.node1() and SecondNode==e.node2()){ 
         return true;
      }
      if(FirstNode==e.node2() and SecondNode==e.node1()){ 
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
      Node SmallThis; 
      Node LargeThis; 
      Node SmallE; 
      Node LargeE; 
      if(FirstNode<SecondNode){
         SmallThis=FirstNode;
         LargeThis=SecondNode;
      }
      else{ 
         SmallThis=SecondNode;
         LargeThis=FirstNode;
      }
      if(e.node1()<e.node2()){
         SmallE=e.node1();
         LargeE=e.node2();
      }
      else{ 
         SmallE=e.node2();
         LargeE=e.node1();
      }


      if(SmallE<SmallThis){ 
         return false; 
      }
      else if(SmallE==SmallThis and LargeE<LargeThis){ 
         return false; 
      }
      (void) e;           // Quiet compiler warning
      return true;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Edge(Node First,Node Second) :FirstNode(First),SecondNode(Second){
    }
    Node FirstNode;
    Node SecondNode;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type Counter=0;
    // O(num_nodes()) constant iterator over number of keys in map
    for(auto it=EdgeIndices.cbegin();it!=EdgeIndices.cend();++it){
       Counter=Counter+(it->second).size();
    }
    return Counter;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
   //O(num_nodes())
    size_type Counter=0; 
    size_type FirstIndex=0; 
    size_type SecondIndex=0;
    // Outer for loop iterates through keys of map O(num_nodes) 
    for(auto iterator=EdgeIndices.cbegin();iterator!=EdgeIndices.cend();++iterator){
       size_type SizeOfNextSet=(iterator->second).size(); 
       Counter=Counter+SizeOfNextSet;
       if(i<Counter){
          FirstIndex=iterator->first;
          Counter=Counter-SizeOfNextSet;
          auto SetIterator=(iterator->second).cbegin();
          //Inner for loop iterates through set of nodes O(num_nodes)
          for(size_type count=Counter;count<i;count++){
             ++SetIterator; 
          }
          SecondIndex=*SetIterator;
          break;
       }
    }
    //Iterate through map until we reach the key we want
    // then iterate through set until we have the appropriate index
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    Node First(this,FirstIndex);
    Node Second(this,SecondIndex); 
    Edge DummyReturn(First,Second); 
    return DummyReturn;        // Invalid Edge
    //return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
   //O(log num_nodes)
   if(b.index()>a.index()){
      //O(log num_nodes) searches keys in map
      auto iter=EdgeIndices.find(b.index());
      //O(log num_nodes) searches nodes in set   
      if((iter->second).find(a.index())!=(iter->second).cend()){
        return true;
     }
   }
   else if(a.index()==b.index()){
     //Self edge always true 
     return true;
   } 
   else{
      //O(log num_nodes) searches keys in map
      auto iter=EdgeIndices.find(a.index());   
      //O(log num_nodes) searches nodes in set   
      if((iter->second).find(b.index())!=(iter->second).cend()){
        return true;
     }
      
   }
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
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
    if(a.index()<b.index()){
       //[]operator is O(log num_nodes)
       // insert is O(log num_nodes) 
       EdgeIndices[a.index()].insert(b.index()); 
    }
    else if (a.index()>b.index()){
       //[]operator is O(log num_nodes)
       // insert is O(log num_nodes) 
       EdgeIndices[b.index()].insert(a.index());
    }
    //Since we are using sets insert only increases the size of set if entry is unique 
    (void) a, (void) b;   // Quiet compiler warning
    Edge DummyReturn(a,b); 
    return DummyReturn;        // Invalid Edge
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
   NodeLocations.clear();
   EdgeIndices.clear();  
   // HW0: YOUR CODE HERE
  }

 private:
  //Store edge as a map of sets. The lower indexed node is the key
  // then there is a set of the higher indexed node it is connected to
  std::map < size_type ,std::set<size_type> > EdgeIndices; 
};

#endif // CME212_GRAPH_HPP
