#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <utility>
#include <vector>

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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** Interpretation of nodes contained in graph. */
  using node_value_type = V;

  /** Interpretation of edges contained in graph*/
  using edge_value_type = E;

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

  /** Struct that stores node info to store in @a nodes container */
  struct node_info {
    Point p;
    node_value_type value;
    int index;
    node_info(Point p_, node_value_type v, size_type i) : p(p_),
                                                          value(v),
                                                          index(i) {
    }
  };

  /** Struct that stores edge info to store in @a edge container */ 
  struct edge_info {
    size_type e_id;
    edge_value_type value;
    edge_info(size_type id, edge_value_type v) : e_id(id),
                                                 value(v) {
    }
  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : number_of_edges(0) {
    
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // PUBLIC FUNCTIONS
  //
  
  /** Remove node @a n and incident edges if it is contained in graph.
   * @pre The index of node @a n must satisfy 0 <= n.index < size(graph). 
   * @return 1 if node @a n and incident edges are removed, 0 otherwise.
   * @post new num_nodes() <= old num_nodes
   * @post new num_edges <= old num_edges 
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  size_type remove_node(const Node& n) {
    if(!has_node(n)) {
      return 0; // Node was not deleted
    }
  
    // Find the correct index of @a n within nodes container
    size_type correct_index = index_to_uid[n.index()]; 
    
    // Delete all connections.
    for(size_type i = 0; i < connections[correct_index].size(); i++) {
      size_type correct_adj = connections[correct_index][i].e_id;
      for(size_type j = 0; j < connections[correct_adj].size(); j++) {
        if(connections[correct_adj][j].e_id == correct_index) {
          // Remove edge(i,j)
          connections[correct_adj].erase(connections[correct_adj].begin() + j);
        }
      }
      --number_of_edges; // Decrement number of edges in graph 
    }
   
    connections[correct_index].clear();
    
    index_to_uid[n.index()] = *(index_to_uid.end() - 1);
    index_to_uid.pop_back();
    nodes[index_to_uid[n.index()]].index = n.index(); // index updated for node moved    
     
    //nodes[correct_index].index = -1; // index now implies this is deleted
     
    return 1; // Return 1 if successfully deleted
  }
  
  /** Remove node @a n and incident edges if it is contained in graph.
   * @pre @a n_it is a valid iterator in graph. 
   * @return n_it.
   * @post new num_nodes() <= old num_nodes
   * @post new num_edges <= old num_edges 
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it; // Iterator pointing to position just deleted from
  }

  /** Remove edge from graph, delete @i.index() from @a j's adjacency container
   * and delete @j.index() from i's adjacency container. 
   * @pre If graph has this edge then 0 <= i.index(), j.index() < num_nodes()
   * @return 1 if edge ij deleted from graph, 0 otherwise. 
   * @post new num_edges <= old num_edges. 
   * @post old num_nodes == new num_nodes
   *
   * Complexity: O(num_nodes() + num_edges()) 
   */
  size_type remove_edge(const Node& i, const Node& j) {
    if(!has_edge(i,j)) return 0;
    
    // Find correct uid's (index within nodes container)
    size_type index_i = index_to_uid[i.index()];
    size_type index_j = index_to_uid[j.index()];
    
    // Remove index_j from @a i's adj list
    for(size_type i = 0; i < connections[index_i].size(); i++) {
      size_type adj = connections[index_i][i].e_id;
      if(adj == index_j) {
        auto remove_iter = connections[index_i].begin() + i;
        connections[index_i].erase(remove_iter);
        break;
      }
    }   

    // Remove index_i from @a j's adj list
    for(size_type j = 0; j < connections[index_j].size(); j++) {
      size_type adj = connections[index_j][j].e_id;
      if(adj == index_i) {
        auto remove_iter = connections[index_j].begin() + j;
        connections[index_j].erase(remove_iter);
        break;
      }
    }
   
    --number_of_edges; // Decrement number of edges in graph
    return 1;
  }

  /** Remove edge from graph, delete @i.index() from @a j's adjacency container
   * and delete @j.index() from i's adjacency container. 
   * @pre If graph has this edge then 0 <= ej.node1().index(), ej.node2().index() < num_nodes()
   * @return value of @a remove_edge(Node i, Node j). 
   * @post new num_edges <= old num_edges. 
   * @post old num_nodes == new num_nodes
   *
   * Complexity: O(num_nodes() + num_edges()) 
   */
  size_type remove_edge(const Edge& ej) {
    return remove_edge(node1(ej), node2(ej));
  }

  /** Remove edge from graph, delete @i.index() from @a j's adjacency container
   * and delete @j.index() from i's adjacency container. 
   * @pre @a *e_it must satisfy preconditions of remove_edge(Edge e).
   * @return iterator @a e_it.  
   * @post new num_edges <= old num_edges. 
   * @post old num_nodes == new num_nodes
   *
   * Complexity: O(num_nodes() + num_edges()) 
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
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
    Point& position() {
      return graph_->nodes[graph_->index_to_uid[index()]].p;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes[graph_->index_to_uid[index()]].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes[uid_].index;
    }
    
    /** Return @a value of node at @a index() */
    node_value_type& value() {
    return graph_->nodes[graph_->index_to_uid[index()]].value; 
    }
    /** Return @a value of node at @a index() in const graph */
    const node_value_type& value() const {
      return graph_->nodes[graph_->index_to_uid[index()]].value; 
    }
    /** Return number of nodes connected to current node (the degree) */
    size_type degree() const {
      return graph_->connections[graph_->index_to_uid[index()]].size(); 
    }
    /** Return IncidentIterator that starts at the first adjacent node to current node */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, graph_->index_to_uid[index()], 0);
    }
    /** Return IncidentIterator that signifies one position past last valid edge */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, graph_->index_to_uid[index()], degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ( (this->uid_ == n.uid_) && (n.graph_ == this->graph_) ) {
        return true;
      }
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
      if ((this->uid_ < n.uid_)) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  
    /** Pointer back to the Graph container */
    graph_type *graph_;
    /** This node's index in Graph container */
    size_type uid_;

    Node(const graph_type* g, size_type nn) 
      : graph_(const_cast<graph_type*>(g)), uid_(nn) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return index_to_uid.size();
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
    node_info vert(position, val, num_nodes()); // About to add node, num_nodes() gives correct index
    nodes.push_back(vert);
    index_to_uid.push_back(nodes.size() - 1); // push back uid (index)
    
    connections.push_back({}); // push back empty vector
    
    return Node(this, nodes.size() - 1);    
  }

  /** Determine if a Node belongs to this Graph
   * @pre @a 0 <= n.uid_ < M where M is the size of the nodes container.
   * @pre index of node @a n satisfies 0<= index < N, where N is the number of nodes in graph.
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /** Using the node invariant property node(i) => node.uid_ = i */
    if (index_to_uid[nodes[n.uid_].index] == n.uid_) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, index_to_uid[i]);
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
    Edge() { 

    }

    edge_value_type& value() {
      return graph_->connections[uid_a][uid_b_index].value;
    }

    const edge_value_type& value() const {
      return graph_->connections[uid_a][uid_b_index].value;
    }
  
    double length() const {
      return norm( node1().position() - node2().position() );
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->graph_, uid_a);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(this->graph_, get_node2_info().e_id);   
    }

    /** Helper function to retreive the id of the 'other' node. */
    edge_info get_node2_info() const {
      return graph_->connections[uid_a][uid_b_index];
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_) {
        return false;
      }
      Node a = node1();
      Node b = node2();

      Node e1 = e.node1();
      Node e2 = e.node2();
      if ( (a == e1 && b == e2) || (a == e2 && b == e1) ) {
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
      auto uid_b = get_node2_info().e_id;
      auto e_uid_b = e.get_node2_info().e_id;
      if ( (&graph_ < &e.graph_) || ((uid_a < e.uid_a) && (uid_b < e_uid_b) )) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type *graph_;
    size_type uid_a;
    size_type uid_b_index;

    Edge(const graph_type* g, size_type na, size_type nb) : graph_(const_cast<graph_type*>(g)), 
                                                            uid_a(na),
                                                            uid_b_index(nb) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return number_of_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *(std::next(edge_begin(), i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // iterate through edges and check if there exists an edge
    // that connects @a a and @a b.
    /** Use Connections map for complexity: O(num_adjacent_nodes) */
    size_type index_a = index_to_uid[a.index()];
    size_type index_b = index_to_uid[b.index()];
    if(connections.size() == 0) return false;
    for(size_type i = 0; i < connections[index_a].size(); i++) {
      if(connections[index_a][i].e_id == index_b) return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) { 
    size_type index_a = index_to_uid[a.index()];
    size_type index_b = index_to_uid[b.index()];
    if ( has_edge(a,b) ) {
      auto adj_vec = connections[index_a];
      auto it = adj_vec.begin();
      while((it != adj_vec.end()) && (it->e_id != index_b)) {
        ++it;
      }
      assert(it != adj_vec.end());
      size_type other_index = std::distance(adj_vec.begin(), it);
      
      return Edge(this, index_a, other_index);
    }
   
    edge_info a_info(index_b, val);
    edge_info b_info(index_a, val); 
    
    connections[index_a].push_back(a_info);
    connections[index_b].push_back(b_info);
    number_of_edges++;

    size_type b_index = connections[index_a].size() - 1; // location within a's adj vector
    return Edge(this, index_a, b_index);   
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    index_to_uid.clear();
    connections.clear();
    number_of_edges = 0;
   }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator(const graph_type *g, size_type x) : graph_(const_cast<graph_type*>(g)), 
                                                     index_(x) {

    }

    /** Return Node that belongs to graph @a graph_ and has index @a index_.
     *
     * Must be a valid node.
    */
    Node operator*() const {
      return Node(graph_, graph_->index_to_uid[index_]);
    }

    /** Increment node index and return pointer to iterator.
     *
     * Nodes must be totally ordered.
    */
    NodeIterator& operator++() {
      this->index_++;
      return *this;
    }

    /** Return true if two node iterators are equal, false otherwise.
     *
     * Nodes belonging to same graph must be totally ordered.
    */
    bool operator==(const NodeIterator& iter) const {
      if(this->graph_ == iter.graph_ && index_ == iter.index_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    /** Pointer back to the Graph container */
    graph_type* graph_;
    /** This node's index in Graph container */
    size_type index_;
  };

  /** Return NodeIterator pointing to the start of node container.
    *
    * Node container must be nonempty.
  */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return NodeIterator that points to one position past last valid Node
    *
    * This node iterator must never be dereferenced.
  */
  node_iterator node_end() const {    
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(const graph_type* g, size_type i, size_type j)
      : graph_(const_cast<graph_type*>(g)), 
        central_node(i),
        adj_node(j) {
    }

    /** Return valid Edge from graph.
     *
     *  Returned edge connects node located at index @a central_node in
     *  adjacency matrix to the node with index stored at @a adj_node in 
     *  adjacency matrix.
     */
    Edge operator*() const {
      return Edge(graph_, central_node, adj_node);
    }

    /** Increment the index of the adjacent node to node at @a central_node
     *  and return IncidentIterator.
     *
     *  adj_node must always be valid when incremented. So it cannot equal
     *  the edge_end() that returns an IncidentIterator.
     *
     */
    IncidentIterator& operator++() {      
      adj_node++;
      return *this;
    }

    /** IncidentIterator equality comparator function. 
     *  
     *  IncidentIterators are only equal if their graphs are equal,
     *  their @a central_nodes are equal, and their @a adj_nodes are equal.
     *
     *  Returns boolean value.
     */
    bool operator==(const IncidentIterator& iter) const {
      if(this->graph_ == iter.graph_ 
        && central_node == iter.central_node
        && adj_node == iter.adj_node) {
        return true;        
      }
      return false;
    }

   private:
    friend class Graph;
    /** Pointer back to the Graph container */
    graph_type *graph_;
    /**  central and adjacent nodes */
    size_type central_node;
    size_type adj_node;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy


    /** Construct an invalid EdgeIterator. */
    EdgeIterator(const graph_type *g, size_type r, size_type c) : 
                                                      graph_(const_cast<graph_type*>(g)), 
                                                      row_(r),
                                                      col_(c) {

    }

    /** Dereference EdgeIterator and return valid Edge. 
     *  
     *  Returns edge at index @a index_. There is no guarantee
     *  on the ordering of the nodes for a given edge.
     *
     *  EX: if Edge(a,b) is returned then node1() and node2()
     *  will each return one or the other. 
     */
    Edge operator*() const {
      return Edge(graph_, graph_->index_to_uid[row_], col_);
    }

    /** Increment and return EdgeIterator to point at the next valid edge. 
     *  
     *  The next valid edge must be valid to be dereferenced. 
     */
    EdgeIterator& operator++() {
      
      ++col_;
      // Fix position of @a col and @a row to valid place in adjacency container 
      while(true) {
        if(row_ >= graph_->num_nodes()) {
          row_ = graph_->index_to_uid[graph_->num_nodes()-1]; 
          col_ = graph_->connections[row_].size(); 
          break;
        }
        else if (col_ >= graph_->connections[graph_->index_to_uid[row_]].size()) {          
          row_++;
          while(row_ < graph_->num_nodes() && graph_->nodes[graph_->index_to_uid[row_]].index == -1) {
            row_++;
          }
          col_ = 0;
        }
        else if (graph_->index_to_uid[row_] > graph_->connections[graph_->index_to_uid[row_]][col_].e_id) {
          col_++;
        }
        else {
          break;
        }
      }
      
      return *this;
    }

    /** Equality comparator for EdgeIterators. 
     *  
     *  EdgeIteratos are equal if their @a graphs_ are the same and 
     *  their @a index_s are equal. 
     */
    bool operator==(const EdgeIterator& ej) const {
      if (graph_ == ej.graph_ && graph_->index_to_uid[row_] == graph_->index_to_uid[ej.row_] && col_ == ej.col_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    graph_type *graph_;
    size_type row_;
    size_type col_;
  };

 /** Returns EdgeIterator that points to the start of a container of edges. 
  *  
  *  Edge at this index must be valid.
  */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, index_to_uid[0], 0);
  }
    
 /** Returns EdgeIterator that points to one index past the last valid edge.
  *  
  *  Returned edge must not be dereferenced. This edge is not valid
  */
  edge_iterator edge_end() const {
    return EdgeIterator(this, index_to_uid[num_nodes()-1], connections[index_to_uid[num_nodes()-1]].size() ); 
  }

 private:


  /** Vector of node_info that represent each node in graph */
  std::vector<node_info> nodes;
  
  /** Vector containing vector of edge_info structs */
  std::vector<std::vector<edge_info>> connections;
  
  /** Number of edges in graph. */
  size_type number_of_edges;
  
  /** Vector that maps node index to its uid */
  std::vector<size_type> index_to_uid;
};

#endif // CME212_GRAPH_HPP
