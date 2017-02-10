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

// Citations/Collaboration:
// Worked with Ian Shaw, Jonathon Roth, and Amery Martin (HW0) directly
// Attended office hours with Amy Shoemaker (TA) (HW0)
// Referenced (and followed the structure closely) of proxy_example.cpp
// Also referenced public github reposititories for previous students in 
// Harvard CS207: Lan, Tian; Piasevoli, Filip; Tran, Dustin; Zacarias, Lisa; 
// **Code was not copied direcly, only referenced for ideas when I was stuck

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
 
template <typename V>
class Graph {
 private:
  // Predeclaring the internal structures: internal_node and internal_edge
    struct internal_node;
 //   struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
    
  /** Type of value in Node */
  using node_value_type = V;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // Construct a Graph with zero size
      size_ = 0; 
      esize_ = 0; // Initialize edge size = 0
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
        return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    // Return the value of the given node
      node_value_type& value(){
          return fetch().value;
      }
      const node_value_type& value() const{
          return fetch().value;
      }
      
    // Return number of edges corresponding to the given node
      size_type degree() const{
          return fetch().neighbors.size();
      }
      
    // Return iterator that points to first edge
      incident_iterator edge_begin() const{
          return IncidentIterator(graph_, uid_, 0);
      }
    
    // Return iterator that points to last edge
      incident_iterator edge_end() const{
          return IncidentIterator(graph_, uid_, degree());
      }
     
      
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Check for same graph and uid
        if (n.graph_ == graph_ && n.uid_ == uid_)
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
      // Assign order by uid
      // Checking for graph equality as per hw0 feedback
        if (n.graph_ == graph_)
            if (uid_ < n.uid_)
                return true;
            return false;
        if (graph_ < n.graph_)
            return true;
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // This space declares private data members and methods for Node
    // that are not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // **Using proxy_example.cpp as guide**:
    // Pointer back to the graph_type container
      graph_type* graph_;
    // This element's unique identification number
      size_type uid_;
    // Private Constructor 
      Node(const graph_type* graph, size_type uid)
          : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
      }
    // Get internal node that is stored in graph
    /** Helper method to return the appropriate node.
     * This checks if the uid is less than graph size and returns the
     * proper uid
     */
      internal_node& fetch() const {
          // We do not wish to loop over all nodes
          // We have commented out the section from proxy_example.cpp
          // we eliminate the for loop here
          //for (size_type i = 0; i < set_->size(); ++i)
              //if (set_->nodes_[i].uid == uid_)
                  //return set_->nodes_[i];
          //assert(false);
          if (uid_ < graph_->size())
              return graph_->nodes[uid_];
          assert(false);
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  // Return Graph's size
    size_type size() const {
        return size_;
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
  Node add_node(const Point& position, const node_value_type& value_in = node_value_type()) { 
    // initialize new_node of type internal_node 
      internal_node new_node;
    // Each element of the new_node has three features: point, uid, and value
      new_node.point = position;
      new_node.uid = size_;
      new_node.value = value_in; 
    // Use push_back function to append the new node to vector nodes
      nodes.push_back(new_node);
    // Increment our size 
      ++size_;
    // Returns a node that points to the new node
      return Node(this, size_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // n.uid is unsigned, so we need not check non-negativity here
    //  if ((n.uid_ < nodes.size()) && (n.index() < size())) 
    //      return true;
    //  return false;
      return (n < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // Check that i is appropriate
    //  assert (i < size());
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
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      // No code needed at this time
    }

    /** Return a node of this Edge */
    Node node1() const {
      // Return appropriate node of the edge corresponding to uid
      //  return graph_->node(graph_->edges[uid_].node1); 
        return graph_-> node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Return appropriate node2 
      //  return graph_->node(graph_->edges[uid_].node2); 
        return graph_-> node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Follow structure for equality operator in Node
        if (graph_ == e.graph_){
            if (node1_ == e.node1_ && node2_ == e.node2_)
                return true;
            if (node1_ == e.node2_ && node2_ == e.node1)
                return true;
        } 
        else{
            return false;
        }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //  Follow structure for ordering as in class Node
      //  Order by uid
      // Check for equality in graph
        if (( graph_ == e.graph_ && node1_ < e.node1_ ) || (graph_ < e.graph_))
            return true;
        return false;
    }

   private:
      friend class Graph;
    // Pointer back to the Graph container
      graph_type* graph_;
    // Two edges' unique identification numbers
      size_type node1_;
      size_type node2_;
    // Private Constructor
      Edge(const graph_type* graph, size_type node1, size_type node2)
          : graph_(const_cast<graph_type*>(graph)), node1_(node1), node2_(node2) {
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // Return the total number of edges
      return esize_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Check if i is appropriate
    //  assert(i < size()); // (from proxy_example)
    // Return new Edge
    //  return Edge(this, i);   
      
    // Here, count will be the number of edges that have been iterated
    // Initialize count at 0 and then increment
      unsigned int count = 0;
      ++i;
    // Here, j is the current node and k is the uid of neighbors
    // Initialize both j and k.
      unsigned int j = 0;
      unsigned int k = 0;
      while (count < i){
        // if k < size of the neighbors of the current node, and
        // if the current node > neighbors
          if (k < nodes[j].neighbors.size() && nodes[j].neighbors[k] > j){
            // Increase number of of edges by 1
              ++count;
          }
        // Move to next uid
          ++k;
        // Check if uid exceeds size of neighbor
          if (k >= nodes[j].neighbors.size()){
            // Restart uid 
              k = 0;
            // Move to next node
              ++j;
          }
      }
    // If we have the last neighbor of a node
      if (k == 0){
          --j;
          return Edge(this, j, nodes[j].neighbors.back());
      }
      else {
          --k;
          return Edge(this, j, nodes[j].neighbors[k]);
      }
          
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Define temporary internal_edge vector to see if a and b form an edge
    // and call this vector is_edge
    //  internal_edge is_edge;
    // Assign index of a and b to node1 and node2 respectively in the is_edge vector
    //  is_edge.node1 = a.index();
    //  is_edge.node2 = b.index();
    // Search the vector edges to see if is_edge exists
   //  if (std::find(edges.begin(), edges.end(), is_edge) == edges.end())
     //     return false;
     // return true;
      
    // Initialize node1_uid and node2_uid
      size_type node1_uid = a.index();
      size_type node2_uid = b.index();
    // Define iterator that finds node b in the neighbors of node a
      auto it = std::find(nodes[node1_uid].neighbors.begin(), nodes[node1_uid].neighbors.end(), node2_uid);
      if (it == nodes[node1_uid].neighbors.end())
          return false;
      return true;
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
    // Returns an Edge that points to a new edge
    // Initialize new_edge of type internal_edge and assign node1 and node2 
    // the values of a and b respectively
    //  internal_edge new_edge;
      size_type node1_uid = a.index();
      size_type node2_uid = b.index();
    // If the edge exists,
      if (has_edge(a,b)) 
          ; // do nothing
      else{ 
    // Otherwise add the edge      
          nodes[node1_uid].neighbors.push_back(node2_uid);
          nodes[node2_uid].neighbors.push_back(node1_uid);
          ++esize_;
      }
    // Return the edge
      if (node1_uid < node2_uid)
          return Edge(this, node1_uid, node2_uid);
      else
          return Edge(this, node2_uid, node1_uid);
  }
      

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      nodes.clear();
      // edges.clear();
      size_ = 0;
      esize_ = 0; //<--------------------------------------------------------------
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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
   
    // Return the node indicated by the iterator
      Node operator*() const{
          return graph_->node(uid_);
      }
    // Find iterator that points to next node  
      NodeIterator& operator++(){
          if (uid_ < graph_->size())
              ++uid_;
          return *this;
      }
      
    // Check if two nodes are equal
      bool operator ==(const NodeIterator& n) const{
          if(graph_== n.graph_ && uid_ == n.uid_)
              return true;
          return false;
      }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
      graph_type* graph_;
      size_type uid_;
      
      NodeIterator(const graph_type* graph, size_type uid)
          : graph_(const_cast<graph_type*>(graph)), uid_(uid){
          }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
    
  // Return the node iterator that points to the first
  // node of the graph
    node_iterator node_begin() const{
        return NodeIterator(this, 0);
    }
    
  // Return the node iterator that points to the last
  // node of the graph
    node_iterator node_end() const{
        return NodeIterator(this, nodes.size());
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
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
      
    // Return the Edge given by incident_iterator
      Edge operator*() const{
          return Edge(graph_, uid_, graph_->nodes[uid_].neighbors[idx_]);
      }
      
    // Return the IncidentIterator that points to the next edge
      IncidentIterator& operator++(){
          if (idx_ == graph_-> nodes[uid_].neighbors.size())
              return *this;
          else{
              ++idx_;
              return *this;
          }
      }
      
    // Return a bool that checks if incident iterator is the same
    // as edge @ e (i.e they point to same thing)
      bool operator==(const IncidentIterator& e) const{
          return (graph_ == e.graph_ && uid_ == e.uid_ && idx_ == e.idx_);
      }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
      graph_type* graph_;
      size_type uid_;
      size_type idx_;
      
      IncidentIterator(const graph_type* graph, size_type uid, size_type idx)
          : graph_(const_cast<graph_type*>(graph)), uid_(uid), idx_(idx){
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
      
    // Return the Edge that is indicated by iterator
      Edge operator*() const{
          return Edge(graph_, node_, graph_->nodes[node_].neighbors[nuid_]);
      }
      
    // Return the EdgeIterator that points to the next edge
      EdgeIterator& operator++(){
          assert(node_ < graph_->size());
        // Go through nodes until next edge is found 
        // Use a do-while loop here to ensure that the loop
        // is executed at least once
          do{
            // If nuid within the neighbors of node, 
            // increment nuid
              if (nuid_ < (graph_->nodes[node_].neighbors.size()-1)){
                  ++nuid_;
              }
            // Otherwise reset nuid and look at next node
              else {
                  nuid_ = 0;
                  ++node_; 
                  while (node_ != graph_->nodes.size() && graph_->nodes[node_].neighbors.size() == 0)
                      ++node_;
              }
          } while (node_ != graph_->nodes.size() && graph_->nodes[node_].neighbors.size() == 0);
        // Return output
          return *this;
      }
      
    // Test if edge_iterator points to same edge object as @a e
      bool operator==(const EdgeIterator& e) const{
          return (graph_ == e.graph_ && node_ == e.node_ && nuid_ == e.nuid_);
      }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
      graph_type* graph_;
      size_type node_;
      size_type nuid_;
      
    // Private constructor
      EdgeIterator(const graph_type* graph, size_type node, size_type nuid)
          : graph_(const_cast<graph_type*>(graph)), node_(node), nuid_(nuid){
          }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
    
    edge_iterator edge_begin() const{
        unsigned int i = 0;
        while (i < nodes.size()){
            if (nodes[i].neighbors.size() != 0)
                return EdgeIterator(this, i, 0);
        }
        return EdgeIterator(this, size(), 0);
    }
    
    edge_iterator edge_end() const{
        return EdgeIterator(this, size(), 0);
    }
    
    

 private:
  // Internal type for set nodes
    struct internal_node {
        Point point;
        size_type uid;
        std::vector<size_type> neighbors;
        node_value_type value;
       
    };
    
  // The following is no longer needed, as we now have a single 
  // structure (internal_node) that can also take in a vector
  // of neighbors, which correspond to edges
    
  // Internal type for set edges
  //  struct internal_edge {
   //     public:
   //         size_type node1;
   //         size_type node2;
    //        size_type uid;
    //        bool operator == (const internal_edge& n) const {
     //           if ((node1 == n.node1 && node2 == n.node2) || 
      //              (node1 == n.node2 && node2 == n.node1))
     ///               return true;
      //          return false;
      //      }
  //  };
        

    std::vector<internal_node> nodes;
   // std::vector<internal_edge> edges;
    size_type size_;
    size_type esize_;

  // Disable copy and assignment of a Graph (from proxy_example.cpp)
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;
};

#endif // CME212_GRAPH_HPP
