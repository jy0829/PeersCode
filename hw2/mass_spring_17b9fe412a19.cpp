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
  double c;        //< Damping constant
  NodeData() : vel(0), mass(1), c(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double length;    //< Edge rest length L_{ij} (initial length)
  double spring;    //< Edge spring constant K_{ij}
  EdgeData() : length(0.0), spring(100.0) {}
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
 * @tparam G::node_value_type contains values (vel, mass)
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a constraint changes the position and velocity of n
 *           if n violates the constraint.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  // Apply the constraints on every node
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    // Check the constraint on n, which updates its position and velocity
    constraint(n, t, dt);
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Object to combine two forces.
 */
template <typename F1, typename F2>
struct CombinedForce2 {
  CombinedForce2(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}

  /** Return the combined gravity force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n, t);
  }

  F1 f1_;
  F2 f2_;
};


/** Object to combine three forces.
 */
template <typename F1, typename F2, typename F3>
struct CombinedForce3 {
  CombinedForce3(F1 f1, F2 f2, F3 f3) : f1_(f1), f2_(f2), f3_(f3) {}

  /** Return the combined gravity force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n, t) + f3_(n, t);
  }

  F1 f1_;
  F2 f2_;
  F3 f3_;
};


/** Combine two forces into one.
 * @param[in] f1 first force to add in
 * @param[in] f2 second force to add in
 *
 * @return the sum of these two forces
 */
template <typename F1, typename F2>
CombinedForce2<F1, F2> make_combined_force(F1 f1, F2 f2) {
  /** Return the combined force of f1 and f2 applied to @a n at time @a t..
   */
  return CombinedForce2<F1, F2>(f1, f2);
}


/** Combine three forces into one.
 * @param[in] f1 first force to add in
 * @param[in] f2 second force to add in
 * @param[in] f3 third force to add in
 *
 * @return the sum of these three forces
 */
template <typename F1, typename F2, typename F3>
CombinedForce3<F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  /** Return the combined force of f1, f2 and f3 applied to @a n at time @a t..
   */
  return CombinedForce3<F1, F2, F3>(f1, f2, f3);
}


/** Gravity Force function object for HW2 #3. */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;  // silent warning
    return Point(0, 0, - grav * n.value().mass);
  }
};


/** Spring Force function object for HW2 #3. */
struct SpringForce {
  /** Return the mass-spring force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;  // silent warning
    Point spring_force = Point(0);
    NODE n2;
    double dist;
    double K;
    double length;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      n2 = (*it).node2();
      K = (*it).value().spring;
      length = (*it).value().length;
      dist = norm(n.position() - n2.position());
      //assert(dist > 0.00001);
      spring_force -= (K / dist) * (dist - length) * (n.position() - n2.position());
    }
    return spring_force;
  }
};


/** Damping Force function object for HW2 #3. */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;  // silent warning
    return - n.value().c * n.value().vel;
  }
};


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;  // silent warning
    if ((n.position() == Point(0)) || (n.position() == Point(1, 0, 0))) {
      return Point(0);
    }
    else {
      Point spring_force = Point(0);
      NODE n2;
      double dist;
      double K;
      double length;
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        n2 = (*it).node2();
        K = (*it).value().spring;
        length = (*it).value().length;
        dist = norm(n.position() - n2.position());
        spring_force -= (K / dist) * (dist - length) * (n.position() - n2.position());
      }
      return Point(0, 0, - grav * n.value().mass) + spring_force;
    }
  }
};


/** Plane Constraint object for HW2#4 */
struct PlaneConstraint {
  /** Update the node's position and velocity if it violates the constraint.
   *
   * For a plane constraint, the violation is if n.position().z < -0.75
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  void operator()(NODE n, double t, double dt) {
    // HW2 #1: YOUR CODE HERE
    (void) t; // silent warning
    if ((n.position() != Point(0)) && (n.position() != Point(1, 0, 0))) {
      n.position() += n.value().vel * dt;
      if (n.position().z < -0.75) {
        n.position().z = -0.75;
        n.value().vel.z = 0.0;
      }
    }
  }
};


/** Sphere Constraint object for HW2#4 */
struct SphereConstraint {
  /** Update the node's position and velocity if it violates the constraint.
   *
   * For a sphere constraint, the violation is if |n.position() - c| < radius
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  void operator()(NODE n, double t, double dt) {
    // HW2 #1: YOUR CODE HERE
    (void) t; // silent warning
    if ((n.position() != Point(0)) && (n.position() != Point(1, 0, 0))) {
      n.position() += n.value().vel * dt;
      Point c = Point(0.5, 0.5, -0.5);
      double radius = 0.15;
      Point cn = n.position() - c;
      if (norm(cn) < radius) {
        Point r = cn / norm(cn);
        n.position() = c + radius * r;
        n.value().vel = n.value().vel - dot(n.value().vel, r) * r;
      }
    }
  }
};


/** Sphere Constraint with node deletion object for HW2#5 */
struct SphereDeleteConstraint {
  /** Delete the node if it violates the constraint.
   *
   * For a sphere constraint, the violation is if |n.position() - c| < radius
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  void operator()(NODE n, double t, double dt) {
    // HW2 #1: YOUR CODE HERE
    (void) t; // silent warning
    if ((n.position() != Point(0)) && (n.position() != Point(1, 0, 0))) {
      n.position() += n.value().vel * dt;
      Point c = Point(0.5, 0.5, -0.5);
      double radius = 0.15;
      Point cn = n.position() - c;
      if (norm(cn) < radius) {
        n.remove_node();
      }
    }
  }
};


int main(int argc, char** argv)
{
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
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().vel = Point(0.0);
    (*it).value().mass = 1.0/graph.num_nodes();
    (*it).value().c = 1.0/graph.num_nodes();
  }

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    for (auto it2 = (*it).edge_begin(); it2 != (*it).edge_end(); ++it2) {
      (*it2).value().spring = 100.0;
      (*it2).value().length = (*it2).length();
    }
  }

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
  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;

        symp_euler_step(graph, t, dt,
                        make_combined_force(GravityForce(), SpringForce(), DampingForce()),
                        SphereDeleteConstraint());

        // Clear everything for SphereDelete
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);

        // Add edges for SphereDelete
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
