#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

template <typename V>
class Graph {

public:
  /** Predeclaration of types, synonyms, sizes. */
  using graph_type = Graph<V>;
  class Node;             using node_type = Node;
  class Edge;             using edge_type = Edge;
  class NodeIterator;     using node_iterator = NodeIterator;
  class EdgeIterator;     using edge_iterator = EdgeIterator;
  class IncidentIterator; using incident_iterator = IncidentIterator;
  using size_type = unsigned;
  using node_value_type = V;

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

    /** Test whether this node and @a n are equal, or < in global order*/
    bool operator==(const Node& n) {
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
    for (size_type i = 0; i < adjacency_list_[a.uid_].size(); ++i) {
      if (b.uid_ == adjacency_list_[a.uid_][i].first) {return true;}
    }
    return false;
  }

  //Add edge if it's not already in the graph
  Edge add_edge(const Node& a, const Node& b) {
    if (has_edge(a, b) == true) {
      size_type adj_size = adjacency_list_[a.uid_].size();
      for (size_type i = 0; i < adj_size-1; ++i) {
        if (b.uid_ == adjacency_list_[a.uid_][i].first) {
          return Edge(this, adjacency_list_[a.uid_][i].second);
        }
      }
      return Edge(this, adjacency_list_[a.uid_][adj_size].second);
    }
    else {
      // Push back nodes to edge vector
      nodepair_.first = a.uid_;
      nodepair_.second = b.uid_;
      edges_.push_back(nodepair_);
      // Push back node-edge pairs to adjacency list
      nodepair_.first = b.uid_;
      nodepair_.second = edges_.size();
      adjacency_list_[a.uid_].push_back(nodepair_);
      nodepair_.first = a.uid_;
      adjacency_list_[b.uid_].push_back(nodepair_);
      return Edge(this, edges_.size());
    }
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
 std::vector< std::pair<size_type, size_type> >
      edges_; // edge vector of pair of node indeces
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
