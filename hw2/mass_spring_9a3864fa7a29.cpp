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
#include <typeinfo>

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
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = unsigned;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports a position-value pair
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraints(g, @a t), where g
 *           is a graph and @a t is the current time. @a constraints are void
 *           functions that affect the position, velocity or existence of a
 *           Node in the graph, along with its incident edges.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {

  constraints(g,t);

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraints(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

double norm(Node& node1, Node& node2) {
  double x = node1.position().x - node2.position().x;
  double y = node1.position().y - node2.position().y;
  double z = node1.position().z - node2.position().z;
  double dist = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
  return dist;
}

double norm(Node& node1) {
  double x = node1.position().x;
  double y = node1.position().y;
  double z = node1.position().z;
  double dist = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
  return dist;
}

Point f_spring(Node& node1, Node& node2, double& K, double& L) {
  Point position_diff = node1.position() - node2.position();
  Point calc = -K*(position_diff)/norm(position_diff)*
               (norm(position_diff) - L);
  return calc;
}

/** Force function object for HW2 Problem #1.
 * This force object sets the nodes at (0,0,0) and (1,0,0) and imposes spring
 * and gravity forces on the remaining nodes.
 */
struct Problem1Force {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }
    double K = 100; double L;
    double m = n.value().mass;
    double gm = -grav*m;
    Point f_spr = Point(0,0,0);
    Node node2;
    for (size_type i = 0; i < n.degree(); ++i) {
      node2 = n.neighbor_node(i);
      L = n.neighbor_node_init_length(i);
      f_spr += f_spring(n, node2, K, L);
    }
    Point f_grav = Point(0,0,gm);
    (void) t;
    return f_spr + f_grav;
  }
};

/** The following three force objects take in a node and time step as
 * arguments and return gravity, spring, and damping forces.
 */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    double m = n.value().mass;
    double gm = -grav*m;
    Point f_grav = Point(0,0,0);
    f_grav = f_grav + Point(0,0,gm);
    (void) t;
    return f_grav;
  }
};

struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    double K = 100;
    double L;
    Point f_spr = Point(0,0,0);
    Node node2;
    for (size_type i = 0; i < n.degree(); ++i) {
      node2 = n.neighbor_node(i);
      L = n.neighbor_node_init_length(i);
      f_spr += f_spring(n, node2, K, L);
    }
    (void) t;
    return f_spr;
  }
};

struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    double c = n.value().mass;
    Point f_dam = -c*n.value().vel;
    (void) t;
    return f_dam;
  }
};

/** This struct takes two forces as arguments and returns the combined force.
 */
template <typename f1, typename f2>
struct MakeCombinedForce {
  f1 force1_;
  f2 force2_;
  MakeCombinedForce<f1, f2>(f1 force1, f2 force2)
    : force1_(force1), force2_(force2) { }
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return force1_(n, t) + force2_(n, t);
  }
};

/** The following two functions take in forces as arguments and returns
 * the combined force. The first function is nested in the second function
 * to allow for three forces to be included.
 */
template <typename f1, typename f2>
MakeCombinedForce<f1, f2> make_combined_force(f1 force1, f2 force2) {
  return MakeCombinedForce<f1, f2>(force1, force2);
}

template <typename f1, typename f2, typename f3>
MakeCombinedForce<MakeCombinedForce<f1, f2>, f3>
    make_combined_force(f1 force1, f2 force2, f3 force3) {
  MakeCombinedForce<f1, f2> force_interim =
    MakeCombinedForce<f1, f2>(force1, force2);
  return MakeCombinedForce<MakeCombinedForce<f1, f2>, f3>(force_interim, force3);
}

/** Following a similar to pattern to above, the following are four constraints
 * that we can set on the position or existence of the nodes and edges.
 * 1) PlaneConstraint: Ensure all nodes/edges stay at or above z = -0.75
 * 2) SphereConstraint: Place a sphere of radius 0.15 at (0.5, 0.5, -0.5)
 *    around which the grid falls
 * 3) SphereConstraint2: Place an identical sphere as in #2, except nodes and
 *    edges disappear when they make contact with the sphere
 * 4) PinConstraint: Set nodes at positions (0,0,0) and (1,0,0) so they don't
 *    move due to gravity, spring or damping forces
 */
struct PlaneConstraint {
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto n = *ni;
      if (n.position().z < -0.75) {
        n.position().z = -0.75;
        n.value().vel.z = 0;
      }
    }
    (void) t;
  }
};

double inner_product(Point& a, Point& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

struct SphereConstraint {
  Point c_ = Point(0.5, 0.5, -0.5);
  double r_ = 0.15;
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto n = *ni;
      Point position_diff = n.position() - c_;
      if (norm(position_diff) < r_) {
        Point R_i = position_diff/norm(position_diff);
        n.position() = r_*(R_i) + c_;
        n.value().vel = n.value().vel -
           inner_product(n.value().vel, R_i) * R_i;
      }
    }
    (void) t;
  }
};

struct SphereConstraint2 {
  Point c_ = Point(0.5, 0.5, -0.5);
  double r_ = 0.15;
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto n = *ni;
      Point position_diff = n.position() - c_;
      if (norm(position_diff) < r_) {
        g.remove_node(n);
      }
    }
    (void) t;
  }
};

struct PinConstraint {
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto n = *ni;
      if (n.position() == Point(0,0,0) ||
          n.position() == Point(1,0,0)) {
        n.value().vel = Point(0,0,0);
      }
    }
    (void) t;
  }
};

/** Again similar to above, the following struct and function allow for
 * multiple constraints to be imposed using a single entry in the
 * symp_euler_step function.
 */
template <typename c1, typename c2>
struct MakeCombinedConstraints {
  c1 constraint1_;
  c2 constraint2_;
  MakeCombinedConstraints<c1, c2>(c1 constraint1, c2 constraint2)
    : constraint1_(constraint1), constraint2_(constraint2) { }
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    constraint1_(g, t);
    constraint2_(g, t);
  }
};

template <typename c1, typename c2>
MakeCombinedConstraints<c1, c2> make_combined_constraint(c1 constraint1, c2 constraint2) {
  return MakeCombinedConstraints<c1, c2>(constraint1, constraint2);
}


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
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions: node mass is 1/N, with N = number of nodes.
  double node_mass = (double)1/graph.num_nodes();
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().mass = node_mass;
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
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;

        // symp_euler_step(graph, t, dt, Problem1Force());

        symp_euler_step(graph, t, dt,
          make_combined_force(GravityForce(),
                              MassSpringForce(),
                              DampingForce()),
          make_combined_constraint(PinConstraint(),
                                   SphereConstraint2()));

        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 1000)
          std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
