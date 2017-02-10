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

    /** Type of this graph. */
    using graph_type = Graph;

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
        const Point &position() const {
            return graph_->nodes_[uid_].first;
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

        /** Return this node's value. */
        node_value_type& value() {
            return graph_->nodes_[uid_].second;
        }

        /** Return this node's value. */
        const node_value_type& value() const {
            return graph_->nodes_[uid_].second;
        }

        /** Return this node's degree. */
        size_type degree() const {
            return graph_->adjacency_[uid_].size();
        }

        /** Return the beginning of this node's incident iterator. */
        incident_iterator edge_begin() const {
            return IncidentIterator(graph_, uid_, 0);
        }

        /** Return the ending of this node's incident iterator. */
        incident_iterator edge_end() const {
            return IncidentIterator(graph_, uid_, graph_->adjacency_[uid_].size());
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            return (graph_ == n.graph_) and (uid_ == n.index());
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node &n) const {
            if (graph_ == n.graph_) {
                return uid_ < n.index();
            }
            return (graph_ < n.graph_);
        }



    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        Graph* graph_;
        size_type uid_;

        /** Valid Node constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a uid_ < graph_->num_nodes()
         */
        Node(const Graph* graph, size_type uid)
                : graph_(const_cast<Graph*>(graph)), uid_(uid) {
            assert(graph_ != NULL);
            assert(uid_ < graph_->num_nodes());
            assert(uid_ >= 0);
        }
    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return nodes_.size();
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
    Node add_node(const Point &position, const node_value_type &node_value = node_value_type()) {
        std::pair <Point, node_value_type> node_info (position, node_value);
        nodes_.push_back(node_info);
        std::vector<size_type> connected_nodes;
        adjacency_.push_back(connected_nodes);
        return Node(this, nodes_.size() - 1);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node &n) const {
        return (this == n.graph_) and (n.index() < size());
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        assert(i < num_nodes());
        assert(i >= 0);
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
    class Edge : private totally_ordered<Edge>{
    public:
        /** Construct an invalid Edge. */
        Edge() {
            // HW0: YOUR CODE HERE
        }

        /** Return a node of this Edge */
        Node node1() const {
            return graph_->node(node1_uid_);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return graph_->node(node2_uid_);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge &e) const {
            bool same_graph = (graph_ == e.graph_);
            bool direction1 = ((node1_uid_ == e.node1_uid_) and (node2_uid_ == e.node2_uid_));
            bool direction2 = ((node1_uid_ == e.node2_uid_) and (node2_uid_ == e.node1_uid_));
            return same_graph and (direction1 or direction2);
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge &e) const {
            if (graph_ == e.graph_) {
                if (node1_uid_ == e.node1_uid_) {
                    return (node2_uid_ < e.node2_uid_);
                }
                return (node1_uid_ < e.node1_uid_);
            }
            return (graph_ < e.graph_);
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        Graph* graph_;
        size_type node1_uid_;
        size_type node2_uid_;

        /** Valid Edge constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a node1_uid_, @a node2_uid_ < graph_->num_nodes()
         * @pre @a node1_uid_ != @a node2_uid_
         */
        Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
                : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {
            assert(graph_ != NULL);
            assert(node1_uid_ != node2_uid_);
            assert(node1_uid_ < graph_->num_nodes());
            assert(node1_uid_ >= 0);
            assert(node2_uid_ < graph_->num_nodes());
            assert(node2_uid_ >= 0);
        }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return std::distance(edge_begin(), edge_end());
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        assert(i < num_edges());
        assert(i >= 0);
        return *std::next(edge_begin(), i);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node &a, const Node &b) const {
        size_type node1_uid = a.index();
        size_type node2_uid = b.index();
        std::vector<size_type> connected_nodes = adjacency_[node1_uid];
        for (size_type j = 0; j < connected_nodes.size(); j++) {
            if (connected_nodes[j] == node2_uid) {
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
    Edge add_edge(const Node &a, const Node &b) {
        size_type node1_uid = a.index();
        size_type node2_uid = b.index();
        if (!has_edge(a, b)) {
            adjacency_[node1_uid].push_back(node2_uid);
            adjacency_[node2_uid].push_back(node1_uid);
        }
        return Edge(this, node1_uid, node2_uid);
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes_.clear();
        adjacency_.clear();
    }

    //
    // Node Iterator
    //

    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator : private totally_ordered<NodeIterator>{
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

        /** Return the node at the node iterator.
         * @pre node_uid_ != graph_->num_nodes()
         */
        Node operator*() const {
            assert(node_uid_ != graph_->num_nodes());
            return graph_->node(node_uid_);
        }

        /** Return the node iterator to the next element. */
        NodeIterator& operator++() {
            ++node_uid_;
            return *this;
        }

        /** Test whether this node iterator and @a ni are equal.
         *
         * Equal node iterators have the same graph and the same node id.
         */
        bool operator==(const NodeIterator& ni) const {
            return (graph_ == ni.graph_) and (node_uid_ == ni.node_uid_);
        }

    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
        Graph* graph_;
        size_type node_uid_;

        /** Valid NodeIterator constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a node_uid_ < graph_->num_nodes()
         */
        NodeIterator(const Graph* graph, size_type node_uid)
                : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid) {
            assert(graph_ != NULL);
            assert(node_uid_ >= 0);
            assert(node_uid_ <= graph_->num_nodes());
        }

    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_iterator node_begin() const
    // node_iterator node_end() const

    /** Return a node iterator pointing to the first element of node sequence.
     * @return NodeIterator pointing to the first element
     * @post If node sequence is empty, dereference is invalid for the return value.
     *
     * Complexity: O(1)
     */
    node_iterator node_begin() const {
        return NodeIterator(this, 0);
    }

    /** Return a node iterator pointing to the past-the-end element in node sequence.
     * @return NodeIterator pointing to the past-the-end element
     * @post If node sequence is empty, dereference is invalid for the return value.
     *
     * Complexity: O(1)
     */
    node_iterator node_end() const {
        return NodeIterator(this, num_nodes());
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

        /** Return the edge at the incident iterator.
         * @pre adj_row_ != graph_->num_nodes()
         */
        Edge operator*() {
            assert(adj_row_ != graph_->num_nodes());
            size_type node1_uid = adj_row_;
            size_type node2_uid = graph_->adjacency_[adj_row_][adj_col_];
            return Edge(graph_, node1_uid, node2_uid);
        }

        /** Return the incident iterator to the next element. */
        IncidentIterator& operator++() {
            ++adj_col_;
            return *this;
        }

        /** Test whether this incident iterator and @a iit are equal.
         *
         * Equal incident iterators have the same graph and the same row
         * and column on the adjacency list.
         */
        bool operator==(const IncidentIterator& iit) const {
            return (graph_ == iit.graph_) and (adj_row_ == iit.adj_row_) and (adj_col_ == iit.adj_col_);
        }

    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE
        Graph* graph_;
        size_type adj_row_;
        size_type adj_col_;

        /** Valid IncidentIterator constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a adj_row_ <= graph_->num_nodes()
         * @pre 0 <= @a adj_col_ <= graph_->adjacency_[adj_row_].size()
         */
        IncidentIterator(const Graph* graph, size_type adj_row, size_type adj_col)
                : graph_(const_cast<Graph*>(graph)), adj_row_(adj_row), adj_col_(adj_col) {
            assert(graph_ != NULL);
            assert(adj_row_ >= 0);
            assert(adj_row_ <= graph_->num_nodes());
            assert(adj_col_ >= 0);
            assert(adj_col_ <= graph_->adjacency_[adj_row_].size());
        }
    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator : private totally_ordered<EdgeIterator>{
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

        /** Return the edge at the edge iterator.
         * @pre adj_row_ != graph_->num_nodes()
         */
        Edge operator*() const {
            assert(adj_row_ != graph_->num_nodes());
            size_type node1_uid = adj_row_;
            size_type node2_uid = graph_->adjacency_[adj_row_][adj_col_];
            return Edge(graph_, node1_uid, node2_uid);
        }

        /** Return the edge iterator to the next element. */
        EdgeIterator& operator++() {
            ++adj_col_;
            fix();
            return *this;
        }

        /** Test whether this edge iterator and @a eit are equal.
         *
         * Equal edge iterators have the same graph and the same row
         * and column on the adjacency list.
         */
        bool operator==(const EdgeIterator& eit) const {
            return (graph_ == eit.graph_) and (adj_row_ == eit.adj_row_) and (adj_col_ == eit.adj_col_);
        }


    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
        Graph* graph_;
        size_type adj_row_;
        size_type adj_col_;

        /** Valid EdgeIterator constructor
         * @pre @a graph_ != NULL
         * @pre 0 <= @a adj_row_ <= graph_->num_nodes()
         * @pre 0 <= @a adj_col_ <= graph_->adjacency_[adj_row_].size()
         * @post The iterator will start at the first valid Edge
         */
        EdgeIterator(const Graph* graph, size_type adj_row)
                : graph_(const_cast<Graph*>(graph)), adj_row_(adj_row), adj_col_(0) {
            assert(graph_ != NULL);
            assert(adj_row_ >= 0);
            assert(adj_row_ <= graph_->num_nodes());
            assert(adj_col_ >= 0);
            assert(adj_col_ <= graph_->adjacency_[adj_row_].size());
            fix();
        }

        /** Move the edge iterator to the next valid element.
         *  (Amy helped me with this. Thanks Amy!)
         */
        void fix() {
            /* Move the row index of adjacency list */
            while (adj_row_ < graph_->adjacency_.size()) {
                /* Move the column index of adjacency list */
                while (adj_col_ < graph_->adjacency_[adj_row_].size()) {
                    /* Edge is valid if node1_uid < node2_uid */
                    size_type node1_uid = adj_row_;
                    size_type node2_uid = graph_->adjacency_[adj_row_][adj_col_];
                    if (node1_uid < node2_uid) {
                        return;
                    }
                    ++adj_col_;
                }
                ++adj_row_;
                adj_col_ = 0;
            }
            return;
        }

    };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // edge_iterator edge_begin() const
    // edge_iterator edge_end() const

    /** Return an edge iterator pointing to the first element of edge sequence.
     * @return EdgeIterator pointing to the first element
     * @post If edge sequence is empty, dereference is invalid for the return value.
     *
     * Complexity: O(1)
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }

    /** Return an edge iterator pointing to the past-the-end element in edge sequence.
     * @return EdgeIterator pointing to the past-the-end element
     * @post If edge sequence is empty, dereference is invalid for the return value.
     *
     * Complexity: O(1)
     */
    edge_iterator edge_end() const {
        return EdgeIterator(this, adjacency_.size());
    }

private:

    std::vector<std::pair<Point, node_value_type>> nodes_;
    std::vector<std::vector<size_type>> adjacency_;

};

#endif // CME212_GRAPH_HPP
