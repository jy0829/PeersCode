/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
    Point vel;       //< Node velocity
    double mass;     //< Node mass
    NodeData() : vel(0), mass(0) {}
};

struct EdgeData {
    double k;     //< Spring constant
    double l;     //< Length
    EdgeData() : k(0), l(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template<typename G, typename F>
double symp_euler_step(G &g, double t, double dt, F force) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }

    return t + dt;
}

//
// Force
//

/** Functor that return the force of gravity */
struct GravityForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        Point f_grav(0, 0, -grav);
        f_grav = n.value().mass * f_grav;
        (void) t;   // silence compiler warnings
        return f_grav;
    }
};

/** Functor that return the spring force */
struct MassSpringForce {
    template<typename NODE>
    Point operator()(NODE n, double t) {
        Point node1_pos = n.position();
        Point f_spring(0, 0, 0);
        for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
            Edge adjacent_edge = *iit;
            Node adjacent_node = adjacent_edge.node2();

            double k = adjacent_edge.value().k;
            double l = adjacent_edge.value().l;

            Point node2_pos = adjacent_node.position();
            Point diff_pos = node1_pos - node2_pos;
            double distance = norm(diff_pos);
            f_spring += -k * (diff_pos / distance) * (distance - l);
        }
        (void) t;    // silence compiler warnings
        return f_spring;
    }
};

/** Functor that return the damping force */
struct DampingForce {

    double damping_const_;

    DampingForce(double damping_const)
            : damping_const_(damping_const) {
    }

    template<typename NODE>
    Point operator()(NODE n, double t) {
        (void) t; // silence compiler warnings
        return -damping_const_ * n.value().vel;
    }

};

/** Functor that return the combination (sum)
 *  of the two forces above.
 */
template<typename F1, typename F2>
struct CombinedForce {

    F1 f1_;
    F2 f2_;

    CombinedForce(F1 f1, F2 f2)
            : f1_(f1), f2_(f2) {
    }

    template<typename NODE>
    Point operator()(NODE n, double t) {
        return f1_(n, t) + f2_(n, t);
    }

};

/** Create a combination of two forces
 * @param[in] f1  Function object defining the force per node
 * @param[in] f2  Function object defining the force per node
 * @return the functor that define force per node equals to f1 + f2
 *
 */
template<typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
    return CombinedForce<F1, F2>(f1, f2);
};

/** Create a combination of three forces
 * @param[in] f1  Function object defining the force per node
 * @param[in] f2  Function object defining the force per node
 * @param[in] f3  Function object defining the force per node
 * @return the functor that define force per node equals to f1 + f2 + f3
 *
 */
template<typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
    return make_combined_force(CombinedForce<F1, F2>(f1, f2), f3);
};

/** Force function object for HW2 #1. */
struct Problem1Force {
    /** Return the force applying to @a n at time @a t.
     *
     * For HW2 #1, this is a combination of mass-spring force and gravity,
     * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
     * model that by returning a zero-valued force. */
    template<typename NODE>
    Point operator()(NODE n, double t) {
        Point node1_pos = n.position();
        if (node1_pos == Point(0, 0, 0) || node1_pos == Point(1, 0, 0)) {
            return Point(0, 0, 0);
        }
        Point f_spring(0, 0, 0);
        Point f_grav(0, 0, -grav);
        f_grav = n.value().mass * f_grav;
        for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
            Edge adjacent_edge = *iit;
            Node adjacent_node = adjacent_edge.node2();

            double k = adjacent_edge.value().k;
            double l = adjacent_edge.value().l;

            Point node2_pos = adjacent_node.position();
            Point diff_pos = node1_pos - node2_pos;
            double distance = norm(diff_pos);

            f_spring += (-1.0) * k * (diff_pos / distance) * (distance - l);
        }
        (void) t;    // silence compiler warnings
        return f_spring + f_grav;
    }
};

//
// Constraint
//

/** Functor that represent the fixed point constraint
 *
 *  Fix a node that violate the constraint by setting the
 *  node velocity to 0.
 *
 *  NOTE: I feel that this is easier for the user if she
 *        wants to define many fixed point constraints
 *        (says 100 points) than having separate constraint
 *        for each point and have to combine them later.
 */
struct FixedConstraint {
    std::vector <Point> fixed_;

    FixedConstraint(const std::vector <Point> fixed)
            : fixed_(fixed) {
    }

    void operator()(GraphType &g, double t) {
        for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
            Node n = *nit;
            auto it = find(fixed_.begin(), fixed_.end(), n.position());
            if (it != fixed_.end()) {
                n.value().vel = Point(0, 0, 0);
            }
        }
        (void) t;
    }
};

/** Functor that represent the plane constraint
 *
 *  Fix a node that violate the constraint by setting the
 *  node z-position to @a plane_ and velocity to 0.
 */
struct PlaneConstraint {

    double plane_;

    PlaneConstraint(double plane)
            : plane_(plane) {
    }

    void operator()(GraphType &g, double t) {
        for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
            Node n = *nit;
            if (dot(n.position(), Point(0, 0, 1)) < plane_) {
                n.position().z = plane_;
                n.value().vel.z = 0.0;
            }
        }
        (void) t;
    }

};

/** Functor that represent the sphere constraint
 *
 *  Fix a node that violate the constraint by setting the
 *  node position to the nearest point on the sphere and
 *  velocity to that of normal to the sphere's surface.
 */
struct SphereConstraint {

    Point center_;
    double radius_;

    SphereConstraint(const Point &center, double radius)
            : center_(center), radius_(radius) {
    }

    void operator()(GraphType &g, double t) {
        for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
            Node n = *nit;
            Point R = (n.position() - center_) / norm(n.position() - center_);
            if (norm(n.position() - center_) < radius_) {
                n.position() = radius_ * R + center_;
                n.value().vel -= dot(n.value().vel, R) * R;
            }
        }
        (void) t;
    }

};

/** Functor that represent the sphere constraint
 *
 *  Fix a node that violate the constraint by remove
 *  the node and edges that hit the sphere.
 */
struct SphereConstraint2 {

    Point center_;
    double radius_;

    SphereConstraint2(const Point &center, double radius)
            : center_(center), radius_(radius) {
    }

    void operator()(GraphType &g, double t) {
        for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
            Node n = *nit;
            if (norm(n.position() - center_) < radius_) {
                g.remove_node(n);
            }
        }
        (void) t;
    }

};

/** Functor that return the combination
 *  of the two constraints above.
 */
template<typename C1, typename C2>
struct CombinedConstraint {

    C1 c1_;
    C2 c2_;

    CombinedConstraint(C1 c1, C2 c2)
            : c1_(c1), c2_(c2) {
    }

    void operator()(GraphType &g, double t) {
        c1_(g, t);
        c2_(g, t);
    }

};

/** Create a combination of two constraints
 * @param[in] c1  Functor defining constraint on graph
 * @param[in] c2  Functor defining constraint on graph
 * @return the functor that represent c1 and c2 constraints
 *
 */
template<typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
    return CombinedConstraint<C1, C2>(c1, c2);
};

/** Create a combination of three constraints
 * @param[in] c1  Functor defining constraint on graph
 * @param[in] c2  Functor defining constraint on graph
 * @param[in] c3  Functor defining constraint on graph
 * @return the functor that represent c1, c2, and c3 constraints
 *
 */
template<typename C1, typename C2, typename C3>
CombinedConstraint<CombinedConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
    return make_combined_constraint(CombinedConstraint<C1, C2>(c1, c2), c3);
};

int main(int argc, char **argv) {
    // Check arguments
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
        exit(1);
    }

    // Construct an empty graph
    GraphType graph;

    // Create a nodes_file from the first input argument
    std::ifstream nodes_file(argv[1]);
    // Interpret each line of the nodes_file as a 3D Point and add to the Graph
    Point p;
    std::vector<typename GraphType::node_type> nodes;
    while (CME212::getline_parsed(nodes_file, p))
        nodes.push_back(graph.add_node(p));

    // Create a tets_file from the second input argument
    std::ifstream tets_file(argv[2]);
    // Interpret each line of the tets_file as four ints which refer to nodes
    std::array<int, 4> t;
    while (CME212::getline_parsed(tets_file, t)) {
        graph.add_edge(nodes[t[0]], nodes[t[1]]);
        graph.add_edge(nodes[t[0]], nodes[t[2]]);
//#if 0
        // Diagonal edges: include as of HW2 #2
        graph.add_edge(nodes[t[0]], nodes[t[3]]);
        graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
        graph.add_edge(nodes[t[1]], nodes[t[3]]);
        graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }

    // HW2 #1 YOUR CODE HERE
    // Set initial conditions for nodes.
    for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit) {
        Node node = *nit;
        node.value().mass = 1.0 / graph.num_nodes();
        node.value().vel = Point(0, 0, 0);
    }

    for (auto iit = graph.edge_begin(); iit != graph.edge_end(); ++iit) {
        Edge edge = *iit;
        edge.value().k = 100.0;
        edge.value().l = norm(edge.node1().position() - edge.node2().position());
    }

    // Set the constraints
    std::vector <Point> fixed;
    fixed.push_back(Point(0, 0, 0));
    fixed.push_back(Point(1, 0, 0));
    FixedConstraint constraint1 = FixedConstraint(fixed);
    PlaneConstraint constraint2 = PlaneConstraint(-0.75);
    SphereConstraint2 constraint3 = SphereConstraint2(Point(0.5, 0.5, -0.5), 0.15);
    auto constraint = make_combined_constraint(constraint1, constraint2, constraint3);

    // Print out the stats
    std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

    // Launch the Viewer
    CME212::SFML_Viewer viewer;
    auto node_map = viewer.empty_node_map(graph);

    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

    viewer.center_view();

    // We want viewer interaction and the simulation at the same time
    // Viewer is thread-safe, so launch the simulation in a child thread
    bool interrupt_sim_thread = false;
    auto sim_thread = std::thread([&]() {

        // Begin the mass-spring simulation
        double dt = 0.001;
        double t_start = 0;
        double t_end = 5;

        for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
            //std::cout << "t = " << t << std::endl;

            // Set up forces
            auto force = make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0 / graph.num_nodes()));

            symp_euler_step(graph, t, dt, force);
            constraint(graph, t);

            viewer.clear();
            node_map.clear();

            // Update viewer with nodes' new positions
            viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
            viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
            viewer.set_label(t);

            // These lines slow down the animation for small graphs, like grid0_*.
            // Feel free to remove them or tweak the constants.
            if (graph.size() < 100)
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

    });  // simulation thread

    viewer.event_loop();

    // If we return from the event loop, we've killed the window.
    interrupt_sim_thread = true;
    sim_thread.join();

    return 0;
}
