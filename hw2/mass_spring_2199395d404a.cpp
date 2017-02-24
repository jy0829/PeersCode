/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <chrono>
#include <fstream>
#include <thread>

#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"
static double N;
// Gravity in meters/sec^2
static constexpr double grav = 9.81;
/** Custom structure of data to store with Nodes */

struct NodeData {
  Point vel;   //< Node velocity
  double mass; //< Node mass
  NodeData() : vel(0), mass(1) {}
};
struct EdgeData {
  double L0;
  double K0;
  EdgeData() : L0(1), K0(100) {}
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
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F>
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
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Force function object for HW2 #1. */
template <typename F1, typename F2> struct ReturnForce2 {
  F1 force1;
  F2 force2;
  template <typename NODE> Point operator()(NODE n, double t) {
    return force1(n, t) + force2(n, t);
  }
};

struct GravityForce {
  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;

    Point grav_force = (n.value().mass) * Point(0, 0, -grav);
    return grav_force;
  }
};
struct MassSpringForce {

  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;

    Point total_force = Point(0);
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      auto n2 = (*i).node2();
      Point L1 = n.position() - n2.position();
      double K1 = (*i).value().K0;
      double L = (*i).value().L0;
      total_force += -(K1 * L1 * (norm(L1) - L)) / norm(L1);
    }
    return total_force;
  }
};

struct DampingForce {
  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    double c = 1.0 / N;
    Point total_force = -c * n.value().vel;
    return total_force;
  }
};

struct Problem1Force {
  /* Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }

    Point grav_force = (n.value().mass) * Point(0, 0, -grav);
    Point total_force = Point(0);
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      auto n2 = (*i).node2();
      Point L1 = n.position() - n2.position();
      double K1 = (*i).value().K0;
      double L = (*i).value().L0;
      total_force += -(K1 * L1 * (norm(L1) - L)) / norm(L1);
    }
    return (total_force + grav_force);
  }
};

template <typename F1, typename F2>
ReturnForce2<F1, F2> make_combined_force(F1 fa1, F2 fa2) {
  return ReturnForce2<F1, F2>{.force1 = fa1, .force2 = fa2};
}

template <typename F1, typename F2, typename F3>
ReturnForce2<ReturnForce2<F1, F2>, F3> make_combined_force(F1 fa1, F2 fa2,
                                                           F3 fa3) {
  return {{fa1, fa2}, fa3};
}

template <typename C1, typename C2> struct ReturnConstr {
  C1 c1;
  C2 c2;
  template <typename Graph> void operator()(Graph &g, double t) {
    c1(g, t);
    c2(g, t);
  }
};

template <typename C1, typename C2>
ReturnConstr<C1, C2> combine_constr(C1 ca1, C2 ca2) {
  return {ca1, ca2};
}

template <typename C1, typename C2, typename C3>
ReturnConstr<ReturnConstr<C1, C2>, C3> combine_constr(C1 ca1, C2 ca2, C3 ca3) {
  return {{ca1, ca2}, ca3};
}

template <typename C1, typename C2, typename C3, typename C4>
ReturnConstr<ReturnConstr<ReturnConstr<C1, C2>, C3>, C4>
combine_constr(C1 ca1, C2 ca2, C3 ca3, C4 ca4) {
  return {{{ca1, ca2}, ca3}, ca4};
}

template <typename C1, typename C2, typename C3, typename C4, typename C5>
ReturnConstr<ReturnConstr<ReturnConstr<ReturnConstr<C1, C2>, C3>, C4>, C5>
combine_constr(C1 ca1, C2 ca2, C3 ca3, C4 ca4, C5 ca5) {
  return {{{{ca1, ca2}, ca3}, ca4}, ca5};
}

template <typename NODE> struct BaseConstraint {
  NODE n;
  Point p;
  template <typename Graph> void operator()(Graph &G, double t) {
    (void)t;
    G.node(n.index()).position() = p;
    G.node(n.index()).value().vel = Point(0, 0, 0);
  }
};

struct PlaneConstraint {
  template <typename Graph> void operator()(Graph &G, double t) {
    (void)t;
    for (auto it = G.node_begin(); it != G.node_end(); ++it) {
      auto n = *it;
      float z = n.position().z;
      if (z < -0.75) {
        n.position().z = -0.75;
        n.value().vel.z = 0;
      }
    }
  }
};

struct SphereConstraint1 {
  template <typename Graph> void operator()(Graph &G, double t) {
    (void)t;
    Point center = Point(0.5, 0.5, -0.5);
    float r = 0.15;
    for (auto it = G.node_begin(); it != G.node_end(); ++it) {
      auto n = *it;
      if (norm(n.position() - center) < r) {
        Point R = (n.position() - center) / norm(n.position() - center);
        n.position() = R * r + center;
        n.value().vel = n.value().vel - dot(n.value().vel, R) * R;
      }
    }
  }
};

struct SphereConstraint2 {
  template <typename Graph> void operator()(Graph &G, double t) {
    (void)t;
    Point center = Point(0.5, 0.5, -0.5);
    float r = 0.15;
    for (auto it = G.node_begin(); it != G.node_end(); ++it) {
      auto n = *it;

      if (norm(n.position() - center) < r) {
        G.remove_node(n);
        // G.remove_node(it);
      }
    }
  }
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

    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions for your nodes, if necessary.
  N = (double)graph.num_nodes();
  BaseConstraint<Node> f1, f2;
  for (auto i = graph.node_begin(); i != graph.node_end(); ++i) {
    (*i).value().mass = 1.0 / N;
    if ((*i).position() == Point(0, 0, 0)) {
      f1.n = *i;
      f1.p = Point(0, 0, 0);
    }

    if ((*i).position() == Point(1, 0, 0)) {
      f2.n = *i;
      f2.p = Point(1, 0, 0);
    }

    for (auto j = (*i).edge_begin(); j != (*i).edge_end(); ++j) {

      double len = (*j).length();
      (*j).value().L0 = len;
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
  auto sim_thread = std::thread([&]() {

    // Begin the mass-spring simulation
    double dt = 0.0005; // 0.001;
    double t_start = 0;
    double t_end = 5; // 5

    for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
      // std::cout << "t = " << t << std::endl;
      // symp_euler_step(graph, t, dt, Problem1Force());
      auto Force = make_combined_force(GravityForce(), MassSpringForce(),
                                       DampingForce());
      symp_euler_step(graph, t, dt, Force);
      combine_constr(f1, f2, PlaneConstraint(), SphereConstraint2(),
                     SphereConstraint1())(graph, t);
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

  }); // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
