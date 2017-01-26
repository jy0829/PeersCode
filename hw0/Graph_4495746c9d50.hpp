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
class Graph {
private:
    // Added private variables to end of Graph
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
     * Return type of Graph::Node::index(), Graph::num_nodes(),
     * Graph::num_edges(), and argument type of Graph::node(size_type)
     */
    using size_type = unsigned;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {
        edge_id = 0;
        node_id = 0;
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
        // Construct an invalid node.
        Node() {
        }

        /** Construct a valid node with given graph and id */
        Node(const Graph* graph, size_type id) {
            this->graph = graph;
            uid = id;
        }

        /** Return this node's position. */
        const Point& position() const {
            return graph->nodes[uid];
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return uid;
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            return n.graph == graph && n.uid == uid;
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node& n) const {
            return uid < n.uid;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        // Pointer to a main graph class holding pointers
        const Graph* graph;

        /** Positive integer representing unique id of node */
        size_type uid; 
    };

    /** Return the number of nodes in the graph.
     * Complexity: O(1).
     */
    size_type size() const {
        return node_id;
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
    Node add_node(const Point& position) {
        nodes.push_back(position);

        // Add node to adjacency with initial empty vector
        adj_list.push_back(std::vector<std::pair<size_type, size_type>>());
        node_id++;  
        return Node(this, node_id-1);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        return this == n.graph && n.uid <= node_id;
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
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
    class Edge {
    public:
        /** Construct an invalid Edge. */
        Edge() {
        }

        /** Construct valid edge with a given graph and id
         * reverse_order specifies whether the edge has the reverse
         * node order than the stored edge
         */
        Edge(const Graph* graph, size_type id, bool reverse_order) {
            this->graph = graph;
            uid = id;
            reverse = reverse_order;
        }

        /** Return a node of this Edge */
        Node node1() const {
            if(reverse){
                return Node(graph, graph->edges[uid].first);
            } else{
                return Node(graph, graph->edges[uid].second);
            }
        }

        /** Return the other node of this Edge */
        Node node2() const {
            if(reverse){
                return Node(graph, graph->edges[uid].second);
            } else{
                return Node(graph, graph->edges[uid].first);
            }
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            return e.graph == graph && e.uid == uid; 
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            return uid < e.uid;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        // Pointer to graph class containing edge data
        const Graph* graph;

        // Unique id of edge object to gather data
        size_type uid;

        // Whether or not edge order is reversed from stored order
        bool reverse;
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return edge_id;
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        return Edge(this, i, false);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        for(std::pair<size_type, size_type>  p : adj_list[a.uid]){
            if(p.first == b.uid){ 
                return true;
            } 
        }
        return false;
    }

    /** Find the index of an edge. If edge does not exist, return -1
     * @pre @a a and @a b are valid nodes of this graph
     * @return index location if for some @a i, edge(@a i) connects @a a and @a b.
     * @return -1 if edge does not exist
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    int edge_index(const Node& a, const Node& b) const {
        for(std::pair<size_type, size_type> p : adj_list[a.uid]){
            if(p.first == b.uid){
                return p.second;
            }
        }
        return -1;
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
        int index = edge_index(a,b); 

        // Edge already exists in graph
        if(index >= 0) {
            if(this->edge(index).node1() == a && this->edge(index).node2() == b){
                return Edge(this, index, false);
            } else{
                return Edge(this, index, true);
            }
        }
        std::pair<size_type, size_type> p(a.uid, b.uid);
        edges.push_back(p);

        // Add edge to both sides of adjacency list
        adj_list[a.uid].push_back(std::pair<size_type, size_type>(b.uid, edge_id));
        adj_list[b.uid].push_back(std::pair<size_type, size_type>(a.uid, edge_id));

        edge_id++;
        return Edge(this, edge_id-1, false);
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes.clear();
        edges.clear();
        adj_list.clear();
        node_id = 0;
        edge_id = 0;
    }

private:
    // Carries position data of nodes
    std::vector<Point> nodes;

    // Carries edge information as a pair of two node indicies
    std::vector<std::pair<size_type, size_type>> edges;

    // Adjacency list as an alternative representation of edges
    std::vector<std::vector<std::pair<size_type, size_type>>> adj_list;

    // Unique node id
    size_type node_id;

    // Unique edge id
    size_type edge_id;
};

#endif // CME212_GRAPH_HPP
