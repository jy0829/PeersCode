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

class Graph {
 private:

 public:

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node & Edge types and synonyms. */
  class Node;
  using node_type = Node;

  class Edge;
  using edge_type = Edge;

  /** Type of indexes and sizes. */
  using size_type = unsigned;

  /** Construct an empty graph. */
  Graph() : nodes_(), nodes_size_(0), next_node_uid_(0),
            edges_(), edges_size_(0), next_edge_uid_(0),
            adjacency_list_() {
  }

  /** Default destructor */
  ~Graph() {
  }

  //
  // NODES
  //

  class Node {
   public:

    /** Construct an invalid node. */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Test whether this node and @a n are equal.*/
    bool operator==(const Node& n) {
      return (graph_ == n.graph_ and uid_ == n.uid_);
    }

    /** Test whether this node is less than @a n in a global order.*/
    bool operator < (const Node& n) {
      if (graph_ == n.graph_ and uid_ < n.uid_) {return true;}
      if (graph_ < n.graph_) {return true;}
      return false;
    }

   private:

    // Pointer back to the Graph container, ID number, constructor
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  // Initiate containers and counters
  std::vector<Point> nodes_; // node vector of points
  size_type nodes_size_; // node vector size counter
  size_type next_node_uid_; // next node counter
  std::vector< std::pair<size_type, size_type> >
       edges_; // edge vector of pair of node indeces
  size_type edges_size_; // edge vector size counter
  size_type next_edge_uid_; // next edge counter
  std::vector< std::vector< std::pair<size_type, size_type> > >
       adjacency_list_; // adjacency list vector
  std::pair<size_type, size_type>
       nodepair_; // container for a pair of nodes
  std::vector< std::pair<size_type, size_type> >
       node_adj_vector_; // container for a node-edge pair

  /** Return the number of nodes in the graph. */
  size_type nodes_size() const {
    return nodes_size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return nodes_size();
  }

  /** Add a node to the graph, returning the added node. */
  Node add_node(const Point& position) {
    nodes_.push_back(position);
    adjacency_list_.push_back(node_adj_vector_);
    ++nodes_size_;
    ++next_node_uid_;
    return Node(this, next_node_uid_-1);
  }

  /** Determine if a Node belongs to this Graph*/
  bool has_node(const Node& n) const {
    return (n.uid_ < num_nodes());
  }

  /** Return the node with index @a i.*/
  Node node(size_type i) const {
    assert(i < nodes_size());
    return Node(this, i);
  }

  //
  // EDGES
  //

  class Edge {

   public:

    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, graph_->edges_[edge_uid_].first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, graph_->edges_[edge_uid_].second);
    }

    /** Return this edge's index, a number in the range [0, graph_size). */
    size_type index() const {
      return edge_uid_;
    }

    /** Test whether this edge and @a e are equal.*/
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_ and edge_uid_ == e.edge_uid_);
    }

    /** Test whether this edge is less than @a e in a global order.*/
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_ and edge_uid_ < e.edge_uid_) {return true;}
      if (graph_ < e.graph_) {return true;}
      return false;
    }

   private:

    // Pointer back to the Graph container, ID number, constructor
    Graph* graph_;
    size_type edge_uid_;
    Edge(const Graph* graph, size_type edge_uid)
        : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid) {
    }

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
  };

  /** Return the total number of edges in the graph.*/
  size_type num_edges() const {
    return edges_size_;
  }

  /** Return the edge with index @a i.*/
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(this, i);
  }

  //Test whether edge is in graph already
  bool has_edge(const Node& a, const Node& b) const {
    for (size_type i = 0; i < adjacency_list_[a.uid_].size(); ++i) {
      if (b.uid_ == adjacency_list_[a.uid_][i].first) {return true;}
    }
    return false;
  }

  //Add edge if it's not already in the graph
  Edge add_edge(const Node& a, const Node& b) {
    if (has_edge(a, b) == true) {
      return Edge(this, next_edge_uid_);
    }
    else {
      // Push back nodes to edge vector
      nodepair_.first = a.uid_;
      nodepair_.second = b.uid_;
      edges_.push_back(nodepair_);
      // Push back node-edge pairs to adjacency list
      nodepair_.first = b.uid_;
      nodepair_.second = next_edge_uid_;
      adjacency_list_[a.uid_].push_back(nodepair_);
      nodepair_.first = a.uid_;
      adjacency_list_[b.uid_].push_back(nodepair_);
      // Increment
      ++edges_size_;
      ++next_edge_uid_;
      return Edge(this, next_edge_uid_-1);
    }
  }

  /** Remove all nodes and edges from this graph.*/
  void clear() {
    nodes_.clear();
    edges_.clear();
    adjacency_list_.clear();
    nodes_size_ = 0;
    next_node_uid_ = 0;
    edges_size_ = 0;
    next_edge_uid_ = 0;
  }

private:

};

#endif // CME212_GRAPH_HPP
