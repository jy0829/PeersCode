#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph{
  public:
  /** Type of this graph. */
  using graph_type = Graph<V>;
  using node_value_type = V;
  using size_type = unsigned;
  private:
    struct internal_nodes {
      Point point;
      node_value_type val;
    };

    std::vector<internal_nodes> nodes;
    std::vector<std::vector<size_type>> adjacency;
    size_type numedges= 0;
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //


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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  // ASSIGNMENT CODE
  Graph()
    : nodes(), adjacency() {}
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return nodes_->nodes[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }
    /**value(): @Return @a val the node value of value_type V.
    * @pre @a uid_ must be less than the number of nodes in the 
    * graph uid_ [0,graph_size).
    * @pre @a must be called on a valid node.
    */
     node_value_type& value() {
       return nodes_->nodes[uid_].val;
     }
    /**const value(): @Return const val the node value of value_type V.
    * @pre @a uid_ must be less than the number of nodes in the 
    * graph uid_ [0,graph_size).
    * @pre @a must be called on a valid node.
    */
     const node_value_type& value() const {
      return nodes_->nodes[uid_].val;
    }
    /**degree(): Returns the number of neighbors of a given node. 
    * @Return of size_type.
    * @pre must be called on a valid node.
    */
    size_type degree() const {
      return nodes_->adjacency[uid_].size();
    }
    /** edge_begin() constructs @a EB, an incident iterator for a 
    *given node pointing the neighbor.
    * @Return of type incident_iterator.
    * @pre node must be valid.
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(nodes_,uid_,0);
    }
    /** edge_edge() constructs @a EB, an incident iterator for a 
    *given node pointing the neighbor.
    * @Return of type incident_iterator.
    * @pre node must be valid.
    */
    incident_iterator edge_end() const {
      return IncidentIterator(nodes_,uid_,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ((uid_ == n.uid_) && (nodes_ == n.nodes_)) {
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
      if ((uid_ < n.uid_)&&(nodes_ == n.nodes_)) {
        return true;
      }
      else if (nodes_ < n.nodes_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* nodes_;
    size_type uid_;
    Node(const graph_type* node_input, size_type uid)
        : nodes_(const_cast<graph_type*>(node_input)), uid_(uid) {}
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
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
  Node add_node(const Point& position, const node_value_type& val_= node_value_type()) {
    internal_nodes addition;
    addition.point = position;
    addition.val = val_;
    nodes.emplace_back(addition);
    adjacency.emplace_back(std::vector<size_type>());
    return Node(this, nodes.size()-1);
  }
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.uid_ < nodes.size() && this == n.nodes_) {return true;}
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return edges_->node(node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return edges_->node(node2_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((e.edges_ == this->edges_)&&
         (((e.node1_id == this->node1_id)&&(e.node2_id==this->node2_id))||
         ((e.node1_id == this->node2_id)&&(e.node2_id==this->node1_id)))) {
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
      if (std::tie(this->edges_,std::min(node1_id,node2_id),std::max(node1_id,node2_id))<
         std::tie(e.edges_,std::min(e.node1_id,e.node2_id),std::max(e.node1_id,e.node2_id))) {
        return true;
      }
      else if (edges_ < e.edges_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* edges_;
    size_type node1_id;
    size_type node2_id;
    Edge(const graph_type* edge_input, size_type node1, size_type node2)
        : edges_(const_cast<graph_type*>(edge_input)),node1_id(node1), node2_id(node2) {}
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return numedges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(),i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // checks for edge in same node order.
    for (unsigned int i = 0; i < adjacency[a.uid_].size(); i++) {
      if (b.uid_ == adjacency[a.uid_][i]) {
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
      if (has_edge(a,b)==true) {
         return Edge(this,a.uid_,b.uid_);
      }
    numedges+=1;
    adjacency[a.uid_].emplace_back(b.uid_);
    adjacency[b.uid_].emplace_back(a.uid_);
    return Edge(this,a.uid_,b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    adjacency.clear();
  }


  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator{
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

    /** NodeIterator, operator*: Returns the node of that the node_iterator is 
    * pointing to. 
    * @Return type node.
    * @pre iterator is != to node_end().
    */
     Node operator*() const {
       return graphs->node(current);
     }
    /** NodeIterator, operator++: Returns NodeIterator after incrimenting to the 
    * next iterator position.
    */
     NodeIterator& operator++() {
       this->current = current+1;
       return *this;
     }
    /** Tests if NodeIterator is in the current graph, and if the NodeIterator
    * comparator and current NodeIterator are pointing to the same node
    * position. 
    * @peram[in]
    * @Return boolean value.
    */
     bool operator==(const NodeIterator& comparator) const {
       if(comparator.current==this->current&&this->graphs == comparator.graphs){
         return true;
       }
       return false;
     }
    /** If operator!= is true, then !operator==.
    * @Return boolean value.  
    */
     bool operator!=(const NodeIterator& n) const {
      return !(*this ==n);
     }
   private:
    friend class Graph;
    Graph* graphs;
    size_type current;
    NodeIterator(const Graph *graphs,size_type p):graphs(const_cast<Graph *>(graphs)),current(p){}
  };
    /**node_begin(): Returns a NodeIterator of the current graph at node 
    * index 0.
    * @Return NodeIterator.
    */
   node_iterator node_begin() const {
     return NodeIterator(this,0);
   }
    /**node_end(): Returns a NodeIterator of the current graph one position  
    * after the last node index.
    * @Return NodeIterator.
    */
   node_iterator node_end() const {
     return NodeIterator(this,size());
   }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
    /** IncidentIterator, operator*: Returns the edge incident to the   
    *  indicated node.
    * @Return type edge.
    * @pre iterator is != to node.edge_end().
    */
     Edge operator*() const {
       return Edge(graphs,node_itr,
              graphs->adjacency[node_itr][curr_itr]);
     }
    /** IncidentIterator, operator++: Returns IncidentIterator after   
    * incrimenting to the next iterator position.
    */
     IncidentIterator& operator++() {
       this->curr_itr = curr_itr+1;
       return *this;
     }
    /** Tests if IncidentIterator is in the current graph, and if the Iterator
    * comparator and current Iterator are pointing to the same edge. 
    * @peram[in]
    * @Return boolean value.
    */
     bool operator==(const IncidentIterator& n) const {
       if(n.curr_itr==this->curr_itr && this->graphs== 
       n.graphs && this->node_itr == n.node_itr){
         return true;
       }
       return false;
     }
    /** If operator!= is true, then !operator==.
    * @Return boolean value.  
    */
     bool operator!=(const IncidentIterator& n) const {
      return !(*this ==n);
     }
   private:
    friend class Graph;
    Graph* graphs;
    size_type node_itr,curr_itr;
    IncidentIterator(const Graph *graphs,size_type p,size_type a):
    graphs(const_cast<Graph *>(graphs)),node_itr(p),curr_itr(a){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator{
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
    /** EdgeIterator, operator*: Returns the edge incident to the   
    *  indicated node.
    * @Return type edge.
    * @pre iterator is != to node.edge_end().
    */
    Edge operator*() const {
       return Edge(graphs,(*N_).index(),((*I_).node2()).index());
    };
    /** EdgeIterator, operator++: Returns EdgeIterator after incrimenting  
    * to the next iterator position.
    */
    EdgeIterator& operator++() {
      ++I_;
      fix();
      return *this;
    };
    /** Tests if EdgeIterator is in the current graph, and if the Iterator
    * comparator and current Iterator are pointing to the same edge. 
    * @peram[in]
    * @Return boolean value.
    */
    bool operator==(const EdgeIterator& e) const {
      if((e.graphs == graphs)&&(e.N_ == N_)
        &&(e.I_ == I_)){
        return true;
      }
      return false;
    };
    /** If operator!= is true, then !operator==.
    * @Return boolean value.  
    */
    bool operator!=(const EdgeIterator& n) const {
      return !(*this ==n);
    };
   private:
    friend class Graph;
    Graph* graphs;
    node_iterator N_ = (*graphs).node_begin();
    incident_iterator I_ = (*N_).edge_begin();
    void fix() {
      while(N_ != (*graphs).node_end()) {
        while(I_ != (*N_).edge_end()) {
          if(*N_ < (*I_).node2()) {
            ++I_;
          }
          else {
            break;
          }
        }
        if(I_ == (*N_).edge_end()) {
          ++N_;
          I_ = (*N_).edge_begin();
        }
        else {
          break;
        }
      }
      if(N_ == (*graphs).node_end()) {
        N_ = (*graphs).node_begin();
        I_ = (*N_).edge_begin();
      }
    };
      EdgeIterator(const Graph *graphs,NodeIterator N,IncidentIterator I)
      :graphs(const_cast<Graph *>(graphs)),N_(N),I_(I){
        fix();
    };
  };
  /**edge_begin(): Returns a EdgeIterator of the current node at edge 
  * index 0.
  * @Return EdgeIterator.
  */
  edge_iterator edge_begin() const{
    return EdgeIterator(this,this->node_begin(),(*node_begin()).edge_begin());
  };
  /**edge_end(): Returns a EdgeIterator of the current node one position  
  * after the last edge index.
  * @Return EdgeIterator.
  */
  edge_iterator edge_end() const{
    return EdgeIterator(this,this->node_end(),(*node_end()).edge_end());
  };
};

#endif // CME212_GRAPH_HPP

