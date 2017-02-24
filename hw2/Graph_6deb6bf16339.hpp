#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>

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

  // HW1
  using node_value_type = V;

  // HW2
  using edge_value_type = E;

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
      // HW0: YOUR CODE HERE
      parentGraph_ = NULL;
      uid_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      return parentGraph_->nodes.at(uid_).point;
    }

    /** Return this node's position. */
    // HW2
    Point& position() {
      return parentGraph_->nodes.at(uid_).point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return parentGraph_->nodes[uid_].index;
    }

    // HW1: YOUR CODE HERE
    /** Return the value of this node */
    node_value_type& value() {
        return parentGraph_->nodes[uid_].n_value;
    }

    /** Return the value of this node as const */
    const node_value_type& value() const {
        return parentGraph_->nodes[uid_].n_value;
    }

    /** Return the degree of this node, the number of edges touching this node */
    size_type degree() const {
        return parentGraph_->nodes[uid_].adjacency_list.size();
    }

    /** Return an iterator to the beginning of the container of incident edges from this node
     * Complexity O(1)
     */
    incident_iterator edge_begin() const {
        incident_iterator it(parentGraph_, uid_, 0);
        return it;
    }

    /** Return an iterator to the end of the container of incident edges from this node
     * Complexity O(1)
     */
    incident_iterator edge_end() const {
        incident_iterator it(parentGraph_, uid_, parentGraph_->nodes[uid_].adjacency_list.size());
        return it;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // nodes are equal if they have the same index
      return parentGraph_!=NULL && parentGraph_==n.parentGraph_ && this->index()==n.index();
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
      // first compare index() if graphs are the same
      if (this->parentGraph_ == n.parentGraph_) {
        return this->index() < n.index();
      } else {
        return this->parentGraph_ < n.parentGraph_;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* parentGraph_;
    size_type uid_;
    Node(Graph* graph, size_type uid): parentGraph_(graph), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return n_index2uid.size();
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
  /* Is this supposed to be removed? new add_node method from HW1 seems to supersede it
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    Node new_node(this, size());;
    nodes.push_back(new_node);
    points.push_back(position);
    return new_node;
  }
  */

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value Value of the node, if not given will be set to type default
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value= node_value_type()) {
      size_type uid = nodes.size();
      Node new_node(this, uid);
      NodeInfo ni;
      ni.node = new_node;
      ni.point = position;
      ni.n_value = value;
      ni.index = n_index2uid.size();
      nodes.push_back(ni);
      n_index2uid.push_back(uid);

      return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      return n.parentGraph_!=nullptr and
             this==n.parentGraph_ and
             n.uid_==n_index2uid[nodes[n.uid_].index];
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return nodes[n_index2uid[i]].node;
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
  class Edge : totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      parentGraph_ = NULL;
      node1_ = 0;
      node2_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return parentGraph_->nodes[node1_].node;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return parentGraph_->nodes[node2_].node;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool direct   = (this->node1() == e.node1()) && (this->node2() == e.node2());
      bool indirect = (this->node1() == e.node2()) && (this->node2() == e.node1());
      return direct || indirect;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        // use node indices to order edges
        size_type this1 = this->node1().index();
        size_type this2 = this->node2().index();
        size_type this_small = std::min(this1, this2);
        size_type this_big = std::max(this1, this2);
        size_type other1 = e.node1().index();
        size_type other2 = e.node2().index();
        size_type other_small = std::min(other1, other2);
        size_type other_big = std::max(other1, other2);
        if (this->parentGraph_ == e.parentGraph_) {
            if (this_small == other_small) {
                return this_big < other_big;
            } else {
                return this_small < other_small;;
            }
        } else {
            return this->parentGraph_ < e.parentGraph_;
        }
    }

    // HW2
    edge_value_type& value() {
        return parentGraph_->edges.at(uid_).e_value;;
    }

    const edge_value_type& value() const {
        return parentGraph_->edges.at(uid_).e_value;;
    }

    double length() const {
        return norm(node1().position() - node2().position());
    }

    size_type index() const {
        return edges[uid_].index;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* parentGraph_;
    size_type node1_; // node uids
    size_type node2_;
    size_type uid_; // edge uid

    Edge(const Graph* parentGraph, const Node node1, const Node node2, const size_type uid) :
      parentGraph_(), node1_(node1.uid_), node2_(node2.uid_), uid_(uid) {
          parentGraph_ = const_cast<Graph*>(parentGraph);
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return e_index2uid.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return edges[e_index2uid[i]].edge;
  }

  /* abstract away almost identical logic in has_edge and add_edge */
  bool has_edge_helper(const Node& a, const Node& b, Edge& e) const {
    if (a.uid_ >= nodes.size()) {
        // not safe to check adjacency_list, so false
        return false;
    }
    auto adj_list = nodes[a.uid_].adjacency_list;
    for (auto it = adj_list.begin(); it != adj_list.end(); ++it) {
        if ( (*it).node2()==b ){
            e = *it;
            return true;
        }
    }
    return false;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      Edge e;
      return has_edge_helper(a, b, e);
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value=edge_value_type()) {
    // HW0: YOUR CODE HERE
    // check if edge already exists
    Edge e;
    if (has_edge_helper(a,b,e)) {
        return e;
    }
    // edge doesn't exist yet
    size_type uid = edges.size();
    Edge new_edge(this, a, b, uid);
    EdgeInfo ei;
    ei.edge = new_edge;
    ei.e_value = value;
    ei.index = e_index2uid.size();
    edges.push_back(ei);
    nodes[a.uid_].adjacency_list.push_back(Edge(this, a, b, uid));
    nodes[b.uid_].adjacency_list.push_back(Edge(this, b, a, uid));
    e_index2uid.push_back(uid);
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear();
    edges.clear();
    n_index2uid.clear();
    e_index2uid.clear();
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
    /** Returns the node at the current position of the iterator
     * @pre @a _pos refers to a valid node index
     * @post no iterator members are changed
     * Complexity O(1)
     */
    Node operator*() const {
        return parentGraph_->nodes[parentGraph_->n_index2uid[index_]].node;
    }

    /** Increments @a _pos and returns itself
     * Complexity O(1)
     */
    NodeIterator& operator++() {
        ++index_;
        return *this;
    }

    /** Returns true if this iterator points to the same index as @a otherIt
     * @post no iterator members are changed
     * Complexity O(1)
     */
    bool operator==(const NodeIterator& otherIt) const {
        return index_ == otherIt.index_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* parentGraph_;
    size_type index_;
    NodeIterator(const Graph* parentGraph) {
        parentGraph_ = const_cast<Graph*>(parentGraph);
    }
  };

  // HW1 #2: YOUR CODE HERE
  /* Returns iterator to the beginning of Container of nodes */
  node_iterator node_begin() const {
      NodeIterator it(this);
      it.index_ = 0;
      return it;
  }

  /* Returns iterator to the end of Container of nodes */
  node_iterator node_end() const {
      NodeIterator it(this);
      it.index_ = n_index2uid.size();
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
    /** Returns the edge at the current position of the iterator
     * @pre @a _startingNodeInd refers to a valid node index
     * @pre @a _pos refers to a valid index into adjacency_list[_startingNodeInd]
     * @post no iterator members are changed
     * @post edge.node1().index() == _startingNodeInd
     * Complexity O(1)
     */
    Edge operator*() const {
        return parentGraph_->nodes[startingNodeUID_].adjacency_list[index_];
    }

    /** Increments @a _pos and returns itself
     * Complexity O(1)
     * */
    IncidentIterator& operator++() {
        ++index_;
        return *this;
    }

    /** Returns true if this iterator points to the same index as @a otherIt
     * @post no iterator members are changed
     * Complexity O(1)
     */
    bool operator==(const IncidentIterator& otherIt) const {
        return index_ == otherIt.index_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* parentGraph_;
    size_type startingNodeUID_;
    size_type index_; // edge index
    IncidentIterator(const Graph* parentGraph, size_type startingNodeUID, size_type index) {
        parentGraph_ = const_cast<Graph*>(parentGraph);
        startingNodeUID_ = startingNodeUID;
        index_ = index;
    }
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
    /** Returns the edge at the current position of the iterator
     * @pre @a _pos refers to a valid edge index
     * @post no iterator members are changed
     * Complexity O(1)
     */
    Edge operator*() const {
        return parentGraph_->edges[parentGraph_->e_index2uid[index_]].edge;
    }

    /** Increments @a _pos and returns itself
     * Complexity O(1)
     * */
    EdgeIterator& operator++() {
        ++index_;
        return *this;
    }

    /** Returns true if this iterator points to the same index as @a otherIt
     * @post no iterator members are changed
     * Complexity O(1)
     */
    bool operator==(const EdgeIterator& otherIt) const {
        return index_ == otherIt.index_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* parentGraph_;
    size_type index_;
    EdgeIterator(const Graph* parentGraph) {
        parentGraph_ = const_cast<Graph*>(parentGraph);
    }
  };

  // HW1 #5: YOUR CODE HERE
  /** Returns iterator to the beginning of Container of edges */
  edge_iterator edge_begin() const {
      edge_iterator it(this);
      it.index_ = 0;
      return it;
  }

  /** Returns iterator to the end of Container of edges */
  edge_iterator edge_end() const {
      edge_iterator it(this);
      it.index_ = e_index2uid.size();
      return it;
  }

  /** removes node
   * @post this->has_node(n)==false
   * @post this->size() == old_this->size() - 1
   * @post nodes[n_index2uid[i]].index == i for all i in n_index2uid
   * invalidates all node iterators, edge iterators, and incident iterators
   * O(num_nodes + d*num_edges)
   */
  size_type remove_node(const Node& n) {
      if (has_node(n)) {
          size_type index = n.index();
          // decrement all indices higher than index
          for (auto it=index+1; it<n_index2uid.size(); ++it) {
              nodes[n_index2uid[it]].index -= 1;
          }
          // remove uid from n_index2uid
          n_index2uid.erase(n_index2uid.begin() + index);
          // remove incident edges
          for (auto it=n.edge_begin(); it!=n.edge_end();) {
              if (!remove_edge(*it)) {
                  ++it;
              }
          }
          return 1;
      }
      return 0;
  }
  node_iterator remove_node(node_iterator n_it) {
      remove_node(*n_it);
      return n_it;
  }
  /** remove edge between nodes n1 and n2
   * @post num_edges() = old_num_edges() - 1
   * @post edges[e_index2uid[i]].index == i for all i in e_index2uid
   * invalidates all edge iterators and incident iterators
   * O(num_edges)
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
      if (has_edge(n1,n2)) {
          Edge e;
          // delete the edge from adjacency list from n1 to n2
          for (auto it=nodes[n1.uid_].adjacency_list.begin(); it!=nodes[n1.uid_].adjacency_list.end(); ++it) {
              if ((*it).node2() == n2) {
                  e = *it;
                  nodes[n1.uid_].adjacency_list.erase(it);
                  break;
              }
          }
          // delete the edge from adjacency list from n2 to n1
          for (auto it=nodes[n2.uid_].adjacency_list.begin(); it!=nodes[n2.uid_].adjacency_list.end(); ++it) {
              if ((*it).node2() == n1) {
                  nodes[n2.uid_].adjacency_list.erase(it);
                  break;
              }
          }
          // decrement all indices higher than index
          size_type index = edges[e.uid_].index;
          for (auto it=index+1; it<e_index2uid.size(); ++it) {
              edges[e_index2uid[it]].index -= 1;
          }
          // remove uid from e_index2uid
          e_index2uid.erase(e_index2uid.begin() + index);
          return 1;
      }
      return 0;
  }
  size_type remove_edge(const Edge& e) {
      Node n1 = e.node1();
      Node n2 = e.node2();
      return remove_edge(n1, n2);
  }
  edge_iterator remove_edge(edge_iterator e_it) {
      Node n1 = (*e_it).node1();
      Node n2 = (*e_it).node2();
      remove_edge(n1, n2);
      return e_it;
  }

 private:
  struct NodeInfo {
      Node node;
      std::vector<Edge> adjacency_list;
      Point point;
      node_value_type n_value;
      size_type index;
  };

  struct EdgeInfo {
      Edge edge;
      edge_value_type e_value;
      size_type index;
  };

  // HW0: YOUR CODE HERE
  // Edges are store reduntantly, once in edges, twice in adjacency_list
  // Much more storage cost, but much more efficient to iteratore over incident edges and all edges
  // adjacency_list[node_ind] is a vector of edges from node "node_ind" to other nodes

  // indexed by node uid
  std::vector<NodeInfo> nodes;
  // indexed by edge uid
  std::vector<EdgeInfo> edges;
  // indexed by node index
  std::vector<size_type> n_index2uid;
  // indexed by edge index
  std::vector<size_type> e_index2uid;


};

#endif // CME212_GRAPH_HPP
