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

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// Spring constant
static constexpr double K = 100.0;

// length
static double init_length;

// length
static double init_c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;   //< Node velocity
  double mass; //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double spring; //< Node velocity
  double len;    //< Node mass
  EdgeData() : spring(K), len(1.0) {}
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
    // if (n.position() == Point (0,0,0) || n.position() == Point(1,0,0))
    //  continue;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // if (n.position() == Point (0,0,0) || n.position() == Point(1,0,0))
    //  continue;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/**
 * FORCES AND COMBINERS
 */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t for HW #1
   */
  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);

    Point f = Point(0, 0, 0);
    for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
      Node n2 = (*iit).node2();
      f += -(*iit).value().spring * (n.position() - n2.position()) *
           ((*iit).length() - (*iit).value().len) / (*iit).length();
    }
    f += n.value().mass * Point(0, 0, -grav);
    return f;
  }
};

struct GravityForce {
  /** Return the gravitational force applying to @a n at time @a t.
   */
  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    return n.value().mass * Point(0, 0, -grav);
  }
};

struct MassSpringForce {
  /** Return the spring force applying to @a n at time @a t.
   */
  template <typename NODE> Point operator()(NODE n, double t) {
    Point f = Point(0, 0, 0);
    for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
      Node n2 = (*iit).node2();
      f += -(*iit).value().spring * (n.position() - n2.position()) *
           ((*iit).length() - (*iit).value().len) / (*iit).length();
    }
    (void)t;
    return f;
  }
};

struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.
   */
  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    return -n.value().vel * init_c;
  }
};

/** Struct to combine forces
 */
template <typename T, typename U> struct CombinedForce {
  T f1;
  U f2;
  template <typename NODE> Point operator()(NODE n, double t) {
    return f1(n, t) + f2(n, t);
  }
};

/** 
 * Functions to accept different number of forces
 */
template <typename T, typename U, typename V>
CombinedForce<CombinedForce<T, U>, V> make_combined_force(T f1, U f2, V f3) {
  return {{f1, f2}, f3};
};

template <typename T, typename U>
CombinedForce<T, U> make_combined_force(T f1, U f2) {
  return {f1, f2};
};

/**
 * CONSTRAINTS AND COMBINERS
 */

/** Struct to contraint to a plane
 */
struct PlaneC {
  template <typename Graph> void operator()(Graph &g, double t) {
    (void)t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i) {
      auto n = *i;
      if (inner_prod(n.position(), Point(0, 0, 1)) < -0.75) {
        n.value().vel.back() = 0.0;
        n.position().back() = -0.75;
      }
    }
  }
};

/** Struct to contraint to a sphere
 */
struct SphereC {
  template <typename Graph> void operator()(Graph &g, double t) {
    (void)t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i) {
      auto n = *i;
      if (norm(n.position() - Point(0.5, 0.5, -0.5)) < 0.15) {

        // calculate R
        Point R = (n.position() - Point(0.5, 0.5, -0.5)) /
                  norm(n.position() - Point(0.5, 0.5, -0.5));

        // update velocity
        n.value().vel -= inner_prod(n.value().vel, R) * R;

        // update position
        n.position() = 0.15 * (n.position() - Point(0.5, 0.5, -0.5)) /
                           norm(n.position() - Point(0.5, 0.5, -0.5)) +
                       Point(0.5, 0.5, -0.5);
      }
    }
  }
};

/** Struct to contraint and remove nodes
 */
struct SphereRemC {
  template <typename Graph> void operator()(Graph &g, double t) {
    (void)t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i) {
      if (norm((*i).position() - Point(0.5, 0.5, -0.5)) < 0.15) {
        g.remove_node(i);
      }
    }
  }
};

/** Struct to contraint to a point
 */
template <typename NODE> struct Fixed {

  NODE n;
  Point pos;

  template <typename Graph> void operator()(Graph &g, double t) {
    (void)t;
    g.node(n.index()).position() = pos;
    g.node(n.index()).value().vel = Point(0, 0, 0);
  }
};

/** Struct to combine two constraints
 */
template <typename T, typename U> struct CombinedCon {
  T C1;
  U C2;
  template <typename Graph> void operator()(Graph &g, double t) {
    C1(g, t);
    C2(g, t);
  }
};

/** Functions to combine different number of forces
 */
template <typename T, typename U, typename V>
CombinedCon<CombinedCon<T, U>, V> make_combined_const(T f1, U f2, V f3) {
  return {{f1, f2}, f3};
};

template <typename T, typename U>
CombinedCon<T, U> make_combined_const(T f1, U f2) {
  return {f1, f2};
};

template <typename T, typename U, typename V, typename W>
CombinedCon<CombinedCon<CombinedCon<T, U>, V>, W>
make_combined_const(T f1, U f2, V f3, W f4) {
  return {{{f1, f2}, f3}, f4};
};

template <typename T, typename U, typename V, typename W, typename X>
CombinedCon<CombinedCon<CombinedCon<CombinedCon<T, U>, V>, W>, X>
make_combined_const(T f1, U f2, V f3, W f4, X f5) {
  return {{{{f1, f2}, f3}, f4}, f5};
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
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    //#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions for your nodes, if necessary.
  Fixed<Node> p1;
  Fixed<Node> p2;
  init_length = norm(graph.node(0).position() - graph.node(1).position());
  init_c = 1.0 / graph.num_nodes();
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    (*ni).value().mass = 1.0 / graph.num_nodes();
    for (auto nii = (*ni).edge_begin(); nii != (*ni).edge_end(); ++nii) {
      (*nii).value().len = (double)(*nii).length();
    }
    if ((*ni).position() == Point(0.0, 0.0, 0.0)) {
      p1.n = (*ni);
      p1.pos = (*ni).position();
    }
    if ((*ni).position() == Point(1.0, 0.0, 0.0)) {
      p2.n = (*ni);
      p2.pos = (*ni).position();
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
    // double dt = 0.0005;
    double dt = 0.001;
    double t_start = 0;
    double t_end = 5.00;

    for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {

      symp_euler_step(graph, t, dt,
                      make_combined_force(GravityForce(), MassSpringForce(),
                                          DampingForce()));
                                          // Problem1Force());

      auto constraint =
          make_combined_const(p1, p2, SphereRemC(), SphereC(), PlaneC());
      constraint(graph, t);

      // Clear the viewer’s nodes and edges
      viewer.clear();
      node_map.clear();
      
      // Update viewer with nodes' new positions and new edges
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