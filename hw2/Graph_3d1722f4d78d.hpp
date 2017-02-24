#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

template <typename V, typename E>
class Graph {

public:
  /** Predeclaration of types, synonyms, sizes. */
  using graph_type = Graph<V, E>;
  class Node;             using node_type = Node;
  class Edge;             using edge_type = Edge;
  class NodeIterator;     using node_iterator = NodeIterator;
  class EdgeIterator;     using edge_iterator = EdgeIterator;
  class IncidentIterator; using incident_iterator = IncidentIterator;
  using size_type = unsigned;
  using node_value_type = V;
  using edge_value_type = E;

  /** Constructor & destructor */
  Graph() {}
  ~Graph() {}

  // NODE CLASS

  class Node : private totally_ordered<Node>{

   public:
    Node() {}

    /** Return this node's position, index and value. */
    const Point& position() const {
      return graph_->nodes_[uid_].first;
    }
    Point& position() {
      return graph_->nodes_[uid_].first;
    }
    size_type index() const {
      return uid_;
    }
    node_value_type& value() {
      return graph_->nodes_[uid_].second;
    }
    const node_value_type& value() const {
      return graph_->nodes_[uid_].second;
    }
    // Items for incident iterator
    size_type degree() const{
      return graph_->adjacency_list_[uid_].size();
    }
    IncidentIterator edge_begin() const {
      return IncidentIterator(graph_, 0, uid_);
    }
    IncidentIterator edge_end() const {
      return IncidentIterator(graph_, graph_->adjacency_list_[uid_].size(), uid_);
    }
    Node neighbor_node(size_type adj_index) const {
      return Node(graph_, graph_->adjacency_list_[uid_][adj_index].first);
    }
    double neighbor_node_init_length(size_type adj_index) const {
      return graph_->edges_lengths_[graph_->adjacency_list_[uid_][adj_index].second];
    }

    /** Test whether this node and @a n are equal, or < in global order*/
    bool operator==(const Node& n) {
      return (graph_ == n.graph_ and uid_ == n.uid_);
    }
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ and uid_ == n.uid_);
    }
    bool operator < (const Node& n) {
      if (graph_ == n.graph_ and uid_ < n.uid_) {return true;}
      if (graph_ < n.graph_) {return true;}
      return false;
    }
    bool operator < (const Node& n) const {
      if (graph_ == n.graph_ and uid_ < n.uid_) {return true;}
      if (graph_ < n.graph_) {return true;}
      return false;
    }

   private:
    Graph* graph_;
    size_type uid_;
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    friend class Graph;
  };

  /** Return the number of nodes in the graph. */
  size_type nodes_size() const { return nodes_.size(); }
  size_type num_nodes() const { return nodes_size(); }
  size_type size() const { return nodes_size(); }

  /** Add a node to the graph, returning the added node. */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    pos_val_pair_.first = position;
    pos_val_pair_.second = node_value;
    nodes_.push_back(pos_val_pair_);
    adjacency_list_.push_back(node_adj_vector_);
    return Node(this, nodes_.size()-1);
  }

  /** Determine if a Node belongs to this Graph*/
  bool has_node(const Node& n) const {
    return (this == n.graph_ and n.uid_ < num_nodes());
  }

  /** Return the node with index @a i.*/
  Node node(size_type i) const {
    assert(i < nodes_.size());
    return Node(this, i);
  }

  /** Find distance between nodes*/
  double norm(const Node& node1, const Node& node2) {
    double x = node1.position().x - node2.position().x;
    double y = node1.position().y - node2.position().y;
    double z = node1.position().z - node2.position().z;
    double dist = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
    return dist;
  }

  // EDGES
  class Edge : private totally_ordered<Edge>{

   public:
    Edge() {}

    /** Return nodes of this Edge, index*/
    Node node1() const {
      return Node(graph_, graph_->edges_[edge_uid_].first);
    }
    Node node2() const {
      return Node(graph_, graph_->edges_[edge_uid_].second);
    }
    size_type index() const {
      return edge_uid_;
    }
    double length() const{
      double x = node1().position().x - node2().position().x;
      double y = node1().position().y - node2().position().y;
      double z = node1().position().z - node2().position().z;
      double dist = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
      return dist;
    }

    /** Test whether this edge and @a e are equal or < in global order.*/
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_ and edge_uid_ == e.edge_uid_);
    }
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_ and edge_uid_ < e.edge_uid_) {return true;}
      if (graph_ < e.graph_) {return true;}
      return false;
    }

   private:
    Graph* graph_;
    size_type edge_uid_;
    Edge(const Graph* graph, size_type edge_uid)
        : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid) {
    }
    friend class Graph;
  };

  /** Return the total number of edges in the graph.*/
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.*/
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(this, i);
  }

  //Test whether edge is in graph already
  bool has_edge(const Node& a, const Node& b) const {
    // for (size_type i = 0; i < adjacency_list_[a.uid_].size(); ++i) {
    for (size_type i = 0; i < a.degree(); ++i) {
      if (b.uid_ == adjacency_list_[a.uid_][i].first) {return true;}
    }
    return false;
  }

  //Add edge if it's not already in the graph
  Edge add_edge(const Node& a, const Node& b) {
    if (has_edge(a, b) == true) {
      // size_type adj_size = adjacency_list_[a.uid_].size();
      for (size_type i = 0; i < a.degree()-1; ++i) {
        if (b.uid_ == adjacency_list_[a.uid_][i].first) {
          return Edge(this, adjacency_list_[a.uid_][i].second);
        }
      }
      return Edge(this, adjacency_list_[a.uid_][a.degree()].second);
    }
    else {
      // Push back nodes to edge vector
      nodepair_.first = a.uid_;
      nodepair_.second = b.uid_;
      edges_.push_back(nodepair_);
      // Push back node-edge pairs to adjacency list
      nodepair_.first = b.uid_;
      nodepair_.second = edges_.size()-1;
      adjacency_list_[a.uid_].push_back(nodepair_);
      nodepair_.first = a.uid_;
      adjacency_list_[b.uid_].push_back(nodepair_);
      // Push back initial distance between nodes to edge length vector
      edges_lengths_.push_back(norm(a,b));
      return Edge(this, edges_.size()-1);
    }
  }

  /* Remove_node function
   * The next 2 functions take in a node object or iterator, removing the
   * relevant node from the graph by removing connected edges, removing the
   * node from containers, and renumbering nodes and edges with higher indeces
   * than the nodes and edges that were removed. This works with complexity
   * O(num_nodes+num_edges) in the node/edge vectors and O(num_nodes*degree/2)
   * in the adjacency list. The remove_edge function nested inside here works
   * in O(num_edges) in the edge vectors and O(num_nodes*degree) in the
   * adjacency list.
   *
   * Preconditions include indicating a valid node in the graph for removal,
   * either passing it directly or as a node iterator object.
   *
   * Postconditions include: all existing nodes remain valid, the node index
   * should remain less than g.num_nodes, and num_edges() returns the number of
   * unique undirected edges. Existing edges not connected to the node remain
   * valid.
   */
  size_type remove_node(const Node& n) {
    // remove edges connected to this node in reverse order in adj_list
    size_type adj_size = adjacency_list_[n.uid_].size();
    for (size_type i = 0; i < adj_size; ++i) {
      size_type idx = adj_size - 1 - i;
      remove_edge(n, n.neighbor_node(idx));
    }
    // remove node from adjacency list and nodes vector
    adjacency_list_.erase(adjacency_list_.begin() + n.index());
    nodes_.erase(nodes_.begin() + n.index());
    // re-index nodes with larger uid_ values in edge vector and adj_list
    for (size_type i = 0; i < edges_.size(); ++i) {
      if (edges_[i].first > n.index()) {
        --edges_[i].first;
      }
      if (edges_[i].second > n.index()) {
        --edges_[i].second;
      }
    }
    for (size_type i = 0; i < adjacency_list_.size(); ++i) {
      for (size_type j = 0; j < adjacency_list_[i].size(); ++j) {
        if (adjacency_list_[i][j].first > n.index()) {
          --adjacency_list_[i][j].first;
        }
      }
    }
    return 0;
  }

  NodeIterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    remove_node(n);
    return 0;
  }

  /* Remove edge function
   * The next 3 functions take in an edge, two nodes, or an edge iterator,
   * remove the relevant edge from the containers and adjacency list, and
   * renumbering edges with higher edge_uid_ numbers. This works with complexity
   * O(num_edges) in the edge vectors and O(num_nodes*degree) in the adjacency
   * list.
   *
   * Preconditions include indicating a valid edge in the graph for removal.
   *
   * Postconditions include: all existing edges remain valid, the edge index
   * should remain less than g.num_edges, and num_edges() returns the number of
   * unique undirected edges.
   */

  size_type remove_edge(const Node& a, const Node& b) {
    // test if edge is in graph
    if (has_edge(a,b) == false) {
      return 0;
    }
    // loop through Node A & B's neighbors to remove edge from adjacency list
    size_type edge_idx = 0;
    for (size_type i = 0; i < adjacency_list_[a.index()].size(); ++i) {
      if (adjacency_list_[a.index()][i].first == b.index()) {
        edge_idx = adjacency_list_[a.index()][i].second;
        adjacency_list_[a.index()].erase(adjacency_list_[a.index()].begin() + i);
      }
    }
    for (size_type i = 0; i < adjacency_list_[b.index()].size(); ++i) {
      if (adjacency_list_[b.index()][i].first == a.index()) {
        adjacency_list_[b.index()].erase(adjacency_list_[b.index()].begin() + i);
      }
    }
    // erase edge from containers and renumber edges with higher edge_uid_ vals
    if (edges_.size() > 0) {
      edges_.erase(edges_.begin() + edge_idx);
    }
    if (edges_lengths_.size() > 0) {
      edges_lengths_.erase(edges_lengths_.begin() + edge_idx);
    }
    for (size_type i = 0; i < adjacency_list_.size(); ++i) {
      for (size_type j = 0; j < adjacency_list_[i].size(); ++j) {
        if (adjacency_list_[i][j].second > edge_idx) {
          --adjacency_list_[i][j].second;
        }
      }
    }
    return 0;
  }

  size_type remove_edge(const Edge& e) {
    Node a = e.node1();
    Node b = e.node2();
    remove_edge(a, b);
    return 0;
  }

  EdgeIterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    remove_edge(e);
    return 0;
  }

  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    NodeIterator() {}

     //Operators and index retrieval
     Node operator*() const {
       return Node(graph_, idx_);
     }
     NodeIterator& operator++() {
       idx_ = idx_ + 1;
       return *this;
     }
     bool operator==(const NodeIterator& n) const {
       return {idx_ == n.idx_ and graph_ == n.graph_};
     }
     size_type index() const {
       return idx_;
     }

   private:
     Graph* graph_;
     size_type idx_;
     NodeIterator(const Graph* graph, size_type idx)
         : graph_(const_cast<Graph*>(graph)), idx_(idx) {
         }

   friend class Graph;
  };

  //Methods for beginning and end of node iterator
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }
  NodeIterator node_end() const {
    return NodeIterator(this, this->nodes_.size());
  }

  class IncidentIterator : private totally_ordered<IncidentIterator>{
   public:
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    IncidentIterator() {}

    //Operators and index retrieval
    Edge operator*() const {
      return Edge(graph_, graph_->adjacency_list_[node_uid_][idx_].second);
    }
    IncidentIterator& operator++() {
      idx_ = idx_ + 1;
      return *this;
    }
    bool operator==(const IncidentIterator& iit) const {
      return {idx_ == iit.idx_ and graph_ == iit.graph_};
    }
    size_type index() const {
      return idx_;
    }

   private:
     Graph* graph_;
     size_type idx_;
     size_type node_uid_;
     IncidentIterator(const Graph* graph, size_type idx, size_type node_uid)
         : graph_(const_cast<Graph*>(graph)), idx_(idx), node_uid_(node_uid) {
         }
     friend class Graph;
  };

  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    EdgeIterator() {}

    //Operators and index retrieval
     Edge operator*() const {
       return Edge(graph_, idx_);
     }
     EdgeIterator& operator++() {
       idx_ = idx_ + 1;
       return *this;
     }
     bool operator==(const EdgeIterator& n) const {
       return {idx_ == n.idx_ and graph_ == n.graph_};
     }
     size_type index() const {
       return idx_;
     }

   private:
     Graph* graph_;
     size_type idx_;
     EdgeIterator(const Graph* graph, size_type idx)
         : graph_(const_cast<Graph*>(graph)), idx_(idx) {
         }
    friend class Graph;
  };

  //Methods for beginning and end of edge iterator
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  EdgeIterator edge_end() const {
    return EdgeIterator(this, this->edges_.size());
  }

/** Remove all nodes and edges from this graph.*/
 void clear() {
   nodes_.clear();
   edges_.clear();
   adjacency_list_.clear();
 }

 private:

 // Initiate containers and counters
 std::vector< std::pair<Point, node_value_type> >
      nodes_; // node vector of positions and values
 std::vector< std::pair<size_type, size_type>>
      edges_; // edge vector of pair of node indeces
 std::vector<edge_value_type>
      edges_lengths_; // initial lengths of edges
 std::vector< std::vector< std::pair<size_type, size_type> > >
      adjacency_list_; // adjacency list
 std::pair<size_type, size_type>
      nodepair_; // container for a pair of nodes
 std::vector< std::pair<size_type, size_type> >
      node_adj_vector_; // container for a node-edge pair
 std::pair<Point, node_value_type>
      pos_val_pair_; // container for a position-value pair
};

#endif // CME212_GRAPH_HPP
