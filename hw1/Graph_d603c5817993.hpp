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

template<typename T>
struct incident_element{
  T edge_id;
  incident_element<T> *next; 
  incident_element(T _edge_id):edge_id(_edge_id),next(nullptr){}
};

template<typename T>
struct incident_list{
  incident_element<T>* head;
  incident_element<T>* tail;
  int size_;
  incident_list():head(nullptr),tail(nullptr),size_(0){}
  int size(){
    return size_;
  };
  void push_back(T& node_id){
    incident_element<T>* current = new incident_element<T>(node_id);
    if(head == nullptr){
      head = current;
      tail = current;
    }else{
      tail->next = current;
      tail = current;
    }
    size_++;
  };    
  T operator[](int i){
    incident_element<T>* x = head;
    for(; i > 0; --i) {x = x->next;}
      return x->edge_id;
  }         
};

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

  /** Type of this graph. */
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
  class Node: private totally_ordered<Node> {
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
      return graph->pos[node_index];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
        return graph->node_value_list[node_index];	
    };
    const node_value_type& value() const{
    	return graph->node_value_list[node_index];
    };
	
    size_type degree() const{
    	return (graph->incident_edge_list[node_index]).size();  
    }
   
    incident_iterator edge_begin() const{
	incident_iterator bedge_iter(graph,node_index,0);   
	return bedge_iter;
    };
    incident_iterator edge_end() const{
   	incident_iterator eedge_iter(graph,node_index,degree());
        return eedge_iter;
    };
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE 
      return n.graph == graph && node_index == n.index() ? true:false ;
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
      return node_index < n.index() ? true:false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph;
    Node(const Graph* graph_ ,const size_type node_index_) 
      :graph(const_cast<Graph*>(graph_)),node_index(node_index_){
    }
    size_type node_index;
   };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return num_node;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    // HW0: YOUR CODE HERE
    return num_node;
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    pos.push_back(position);
    node_value_list.push_back(0);
    num_node++;
    incident_list<size_type> new_incident_list;
    incident_edge_list.push_back(new_incident_list);
    return Node(this,num_node-1);        // Invalid node
  }

  
  Node add_node(const Point& position, const node_value_type& value){
     Node new_node(this,num_node);
     pos.push_back(position);
     node_value_list.push_back(value);
     incident_list<size_type> new_incident_list;
     incident_edge_list.push_back(new_incident_list);
     num_node++;
     return Node(this,num_node-1);
  } 
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.graph == this ? true:false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0 && i < num_nodes());
    return Node(this,i);        
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
  class Edge: private totally_ordered<Edge> {
    public:
      /** Construct an invalid Edge. */
      Edge() {
      }

      /** Return a node of this Edge */
      Node node1() const {
        // HW0: YOUR CODE HERE   
        return Node(graph,graph->edges[edge_index].first);
      }

      /** Return the other node of this Edge */
      Node node2() const {
        // HW0: YOUR CODE HERE
        return Node(graph,graph->edges[edge_index].second);   
      }

      /** Test whether this edge and @a e are equal.
      *
      * Equal edges represent the same undirected edge between two nodes.
      */
      bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE
        return e.graph == graph && e.edge_index == edge_index ? true : false;
      }

      /** Test whether this edge is less than @a e in a global order.
       *
       * This ordering function is useful for STL containers such as
       * std::map<>. It need not have any interpretive meaning.
       */
      bool operator<(const Edge& e) const {
        // HW0: YOUR CODE HERE
        return edge_index < e.edge_index;
      }


      private:
      // Allow Graph to access Edge's private member data and functions.
      friend class Graph;
      // HW0: YOUR CODE HERE
      // Use this space to declare private data members and methods for Edge
      // that will not be visible to users, but may be useful within Graph.
      // i.e. Graph needs a way to construct valid Edge objects
      Graph* graph; 
      size_type edge_index;
      Edge(const Graph* graph_ ,const size_type edge_index_) 
      :graph(const_cast<Graph*>(graph_)),edge_index(edge_index_){
      }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
      // HW0: YOUR CODE HERE
      return num_edge;
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
      // HW0: YOUR CODE HERE 
      assert(i >= 0 && i < num_edges()); 
      return Edge(this,i);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
     bool has_edge(const Node& a, const Node& b) const {
       // HW0: YOUR CODE HERE 
       for(size_type i = 0; i < num_edge; ++i){
         if((edges[i].first == a.index() && edges[i].second == b.index()) || (edges[i].first ==b.index() && edges[i].second == a.index()))
         {
           return true;
         } 
       } 
       return false;
     }

     Edge edge_first_node(size_type edge_id, size_type node_id){
         size_type first = edges[edge_id].first;
         if(first != node_id){
            assert(edges[edge_id].second == node_id);
            edges[edge_id].first = edges[edge_id].second; 
            edges[edge_id].second = first; 
         }
         return Edge(this,edge_id);
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
       if(has_edge(a,b)){
       	 for(size_type i = 0; i < num_edge; ++i){
            if((edges[i].first == a.index() && edges[i].second == b.index()))
            {
               return Edge(this,i);
            } 
            else if((edges[i].first == b.index() && edges[i].second == a.index())){
               size_type first = edges[i].first;
               edges[i].first = edges[i].second; 
               edges[i].second = first;
               return Edge(this,i);
            }
         } 
       }
       std::pair<size_type,size_type> thisPair (a.index(),b.index());
       edges.push_back(thisPair);
       Edge e = Edge(this,num_edge); 
       incident_edge_list[a.index()].push_back(num_edge);
       incident_edge_list[b.index()].push_back(num_edge); 
       num_edge++;
       return e;
     }
 

     /** Remove all nodes and edges from this graph.
      * @post num_nodes() == 0 && num_edges() == 0
      *
      * Invalidates all outstanding Node and Edge objects.
      */
     void clear() {
      // HW0: YOUR CODE HERE
       num_node = 0;
       num_edge = 0;
       pos.empty();
       node1.empty();
       node2.empty();
     }
      
      //
      // Node Iterator
      //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private totally_ordered<NodeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const{
	return graph->node(node_iter); 
    };
    NodeIterator& operator++(){
        node_iter = node_iter + 1;
	return *this;
    };
    bool operator==(const NodeIterator& n) const{
    	return (graph == n.graph) && (node_iter == n.node_iter)? true:false;
    };

    size_type node_id(){
       return node_iter;
     }
   private:
     friend class Graph;
      // HW1 #2: YOUR CODE HERE
     Graph* graph; 
     NodeIterator(const Graph* graph_ ,size_type n_iter) 
      :graph(const_cast<Graph*>(graph_)),node_iter(n_iter){
      }
    size_type node_iter;
 };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   NodeIterator node_begin() const{
      NodeIterator biter(this,0);
      return biter;
   };
   NodeIterator node_end() const{ 
      NodeIterator eiter(this,num_node);
      return eiter;
   }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator>  {
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
    Edge operator*() const{
       auto edge_id_ = (graph->incident_edge_list[node_id][incident_edge_id]);     
       auto _edge = graph->edge_first_node(edge_id_,node_id);
      return _edge;
    };
    
    IncidentIterator& operator++(){
	incident_edge_id++;
        return *this;
    };
    
    bool operator==(const IncidentIterator& n) const{
	return (graph == n.graph) && (node_id == n.node_id) && (incident_edge_id == n.incident_edge_id)? true:false;
    };

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph; 
    IncidentIterator(const Graph* graph_,size_type n_id, size_type i_id) 
      :graph(const_cast<Graph*>(graph_)),node_id(n_id),incident_edge_id(i_id){
    }
    size_type node_id;
    size_type incident_edge_id;
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
    Edge operator*() const{
        return graph->edge(edge_id);
    };

    EdgeIterator& operator++(){
        edge_id++;
        return *this;
    };
    bool operator==(const EdgeIterator& n) const{
       return (graph == n.graph) && (edge_id == n.edge_id);
    };

   private:
     friend class Graph;
      // W1 #5: YOUR CODE HERE
    Graph* graph; 
    EdgeIterator(const Graph* graph_,size_type e_id) 
      :graph(const_cast<Graph*>(graph_)),edge_id(e_id){
    }
    size_type edge_id;
   };

    // HW1 #5: YOUR CODE HERE
    //Supply definitions AND SPECIFICATIONS for:
   edge_iterator edge_begin() const{
     EdgeIterator biter(this,0);
     return biter;
    }
   edge_iterator edge_end() const{
     EdgeIterator eiter(this,num_edge);
      return eiter;
   }
  
  private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    size_type num_node = 0;
    size_type num_edge = 0;
    std::vector<Point> pos;
    std::vector<Edge> edge_list;
    std::vector<std::pair<size_type,size_type> > edges;
    std::vector<Node> node1;
    std::vector<Node> node2;
    std::vector<node_value_type> node_value_list;
    std::vector<incident_list<size_type> > incident_edge_list;
};
#endif // CME212_GRAPH_HPP
