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
 //Add Hw 1 Part 1
template <typename V>
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

    //Hw 1 Part 1
  using node_value_type = V;



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
      // HW0: YOUR CODE HERE

    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE

      return graph_->pts[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE

      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value. */
    /** Remark: I store the nodes' position and value in the the tuple position, value. */
    node_value_type& value(){
      return graph_->alpha_nodes[uid_].second;
     }

    /** Return this node's value as a const, aka unmodifiable. */
    /** Remark : extact part 2 of the tuple alpha_ndoes*/
    const node_value_type& value() const {
      return graph_->alpha_nodes[uid_].second;
     }

    /** Return the number of nodes which are directly connected to that node */
    size_type degree() const {
      return graph_->adjList[uid_].size();
    }

    /** Returns an IncidentIterator type iterator which points to the first
     * edge which is adjacent to this node (edge of: this_node, first node adjacent to that edge)*/
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,uid_,0);
    }

    /** Returns an IncidentIterator type iterator which does NOT points to
     * any edge adjacent to this node because degree = index for last edge + 1
     * aka returns a null_pointer*/
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,uid_,degree());
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE

      //COMPARE IF NOT THE SAME GRAPH - ADD IT
      //There are 2 cases - if on different graphs and if on same graph
      if (n.graph_ != graph_) return (n.graph_ < graph_);
      return (n.uid_ < uid_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // A node is defined by its pointer to the graph, and its id
    Graph* graph_;
      // This element's unique identification number
      size_type uid_;
      /** Private Constructor */
      Node(const Graph* graph, size_type uid)
              : graph_(const_cast<Graph*>(graph)), uid_(uid) {
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return pts.size();
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

  Node add_node(const Point& position, const node_value_type& my_val = node_value_type()) {

        //Hw 1 Ex 1 - added value to a vector of nodes
        // Add a vector of pairs (position, value)
        // create separetely the pairs
        alpha_nodes.push_back(std::make_pair(position, my_val));

        // HW0: YOUR CODE HERE
        //when we add a node, we need to add to the - vector of nodes, but also
        // to the edges via the adjoined List
        pts.push_back(position);
        adjList.push_back(std::vector<size_type>());
        return Node(this, (size()-1));
    }

   //OLD IMPLEMENTATION from hw 0 that was modified from 1
   /**Node add_node(const Point& position) {
     HW0: YOUR CODE HERE

    // when we add a node, we need to add to the - vector of nodes, but also
    // to the edges via the adjoined List
    pts.push_back(position);
    adjList.push_back(std::vector<size_type>());
    return Node(this, (size()-1));
  } */

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    return (n.graph_ == this && (n.uid_ >=0 && n.uid_<=size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

    assert(i < size());
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->node(uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->node(uid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //written code though no indication
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
      //no need to write code here - by the hw's indications; wrote it anyways

      //IDEA - look, in this order, at both graphs and nodes
      //If different graphs, look directly at the graphs id's
      // If same graphs, check first by minumum, and if same minimum, by maximum

      //If different graph, simiply look at their graphs ids
      if (graph_ != e.graph_) {return graph_ < e.graph_;}

      // if on same graph, check the min and max between 2 nodes
      if (graph_ == e.graph_)
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
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // An Edge is defined by its pointer to a gprah, and its 2 nodes(id's of nodes)
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
    // HW0: YOUR CODE HERE

    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    auto id1 = edges[i].first;
    auto id2 = edges[i].second;
    return Edge(this, id1, id2);

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //chevk via the adjointed list, where the adjointed list is like a vector
    // of vector; to remember later in terms of x[i][j], in terms of the search
    for (size_type i=0; i<adjList[a.uid_].size(); i++) {
        if (adjList[a.uid_][i]==b.uid_) return true;
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
    //Case 1: if the add we want to add already exists
    if (has_edge(a, b)) {return Edge(this, a.uid_, b.uid_);}
    // Case 2: if the edge does not exist
    else {
          // Add to the adjointed list
          adjList[a.uid_].push_back(b.uid_);
          adjList[b.uid_].push_back(a.uid_);
          // Add also to the pairs, in a sorted way
          if (a.uid_ < b.uid_) edges.push_back(std::make_pair(a.uid_, b.uid_));
          else edges.push_back(std::make_pair(b.uid_, a.uid_));

          return Edge(this, a.uid_, b.uid_);
          }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    edges.clear();
    pts.clear();
    adjList.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {



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


    /** Dereferencing the node
     * Returns node which points @graph and has index @index_points *
     * pre condition: 0<=index_points<number of nodes*/
    Node operator*() const {
      return graph_->node(index_points);
    }

    /** Increment the NodeIterator to point to the next Node  */
    /** Returns a NodeIterator which thus points to the next Node
     * or to the null pointer if the index_point was originally the last index*/
    NodeIterator& operator++() {
      index_points++;
      return *this;
    }

    /** Returns a boolean which checks if 2 node iterators are the same , looking
     * to see if they  :
     * 1) point to the same graph 2) have the same index*/
    bool operator==(const NodeIterator& iter) const {
      return (iter.graph_ == graph_) && (iter.index_points == index_points);
    }


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // The node iterator needs 2 key features - graph to which it is pointing, and index
    const Graph *graph_;
    size_type index_points;

    //Returns a NodeIterator which is pointing out to graph g with index i
    // where i is between 0<= and <number of nodes
    NodeIterator(const Graph* g, size_type i): graph_(g), index_points(i) {}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

    //Node Iterator itself is a pointer
    //use the constructor
    /** Returns a Node Iterator that points to the first Node in the graph     */
    node_iterator node_begin() const {
      NodeIterator it(this,0);
      return it;
    }

    /** Returns a Node Iterator that points to the nullpointer     */
    node_iterator node_end() const
    {
      NodeIterator it(this,alpha_nodes.size());
      return it;
    }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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

     /**An incident iterator will go through each node, and for each node will be able to tell me
      * what are are the edges adjacent to that node, think in terms of (node_index; edge_index)
      *  IncidedentIterator belongs to the graph, is at the node index and points to the edge index
      */

      /** Dereference an Incident Iterator to return an edge to which
       * this incident_iterator points to
       * */
      Edge operator*() const
      {
       return Edge(graph_, node_index_,graph_->adjList[node_index_][edge_index_]);
      }

      /** Increment operator */
      /** Increment the NodeIterator to point to the next adjacent Node or nullpointer  */
      /** Returns a NodeIterator which thus points to the next  adjacent NOde
       * or to the null pointer if the edge_index was originally the last index*/
      IncidentIterator& operator++()
      {
        edge_index_++;
        return *this;
      }

      /** Returns a boolean which checks if 2 incident_iterators are the same , looking
        * to see if they:
        * 1) point to the same graph
        * 2) point to the same node
        * 3) and that node points exactly to the same edge */
      bool operator==(const IncidentIterator& Inc_Iter) const
      {
        return (Inc_Iter.node_index_ == node_index_) &&
               (Inc_Iter.edge_index_ == edge_index_) &&
               (Inc_Iter.graph_ == graph_);
      }

  private:
    friend class Graph;

      size_type node_index_;
      size_type edge_index_;
      Graph *graph_;

      /** Constructor of type Incident Iterator which stores
       * 1. pointer to the graph
       * 2. node index
       * 3. edge index
       * (this is useful because an edge will be (node index, edge index) ("approximatelly" - see deference
       *  code for the exact version)*/
      IncidentIterator(const Graph* g, size_type node_index, size_type edge_index)
              : graph_(const_cast<Graph*>(g)), node_index_(node_index), edge_index_(edge_index) {}
      // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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


    /** Dereferencing the edge
    * Returns edge which points @graph and has index @index_edge *
    * pre condition: 0<=index_edge<number of edges*/
    Edge operator*() const {
       return graph_->edge(edge_index_);
      }

    /** Increment the EdgeIterator to point to the next Edge  */
    /** Returns an EdgeIterator which thus points to the next Edge */
    EdgeIterator& operator++() {
       edge_index_++;
       return *this;
      }

    /** Returns a boolean which checks if 2 edge iterators are the same , looking
     * to see if they:
     * 1) point to the same graph 2) same edge index*/
    bool operator==(const EdgeIterator& EdgeIter) const {
          return (EdgeIter.graph_ == graph_) && (EdgeIter.edge_index_ == edge_index_);
      }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

      size_type edge_index_;
      Graph *graph_;

      //Returns an EdgeIterator which is poiting out to graph g with index i
      // where i is between 0<= and <number of edges
      //Remark for myself - analogy to the node iterator
      EdgeIterator(const Graph* g, size_type edge_index)
              : graph_(const_cast<Graph*>(g)), edge_index_(edge_index) {}


  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

    /** Returns an Edge Iterator that points to the first Edge in the graph     */
    edge_iterator edge_begin() const
    {
      EdgeIterator it(this,0);
      return it;
    }

    /** Returns an Edge Iterator that points to the nullpointer     */
    edge_iterator edge_end() const
    {
      EdgeIterator it(this,edges.size());
      return it;
    }

 private:

    std::vector<Point> pts; //vector to store the nodes(points)
    std::vector<std::vector<size_type>> adjList; //Store the nodes in adjointed list format
    //aka for each node, we all know all nodes to which it is connected
    std::vector<std::pair<size_type, size_type>> edges; //store simple pairs of nodes, aka
    //"simple edges" - (i,j)

    //Hw 1 Ex 1
    std:: vector<std::pair<Point, V>> alpha_nodes; // store for the nodes both the position and value

};

#endif // CME212_GRAPH_HPP
