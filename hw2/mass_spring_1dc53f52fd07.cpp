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

#include <functional>

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;                     //< Spring constant
  double L;                     //< Spring rest-length
  EdgeData() : K(100), L(0) {}
};


// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


//
// UPDATE STEPS
//

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and given constraint.
 * @param[in,out] g           Graph
 * @param[in]     t           The current time (useful for time-dependent forces)
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    
    // Ignore fixed points
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0))
      continue;
    
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Ignore fixed points
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0))
      continue;
    
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  
  return t + dt;
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and given constraint.
 * @param[in,out] g           Graph
 * @param[in]     t           The current time (useful for time-dependent forces)
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @param[in]     constraint  Constraint object applied on a node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint (n, @a t, @a g).
 *           where n is a node of the graph and @a t is the current time.
 *           @a constraint changes the node value without returning anything.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    
    // Ignore fixed points
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0))
      continue;
    
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Ignore fixed points
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0))
      continue;
    
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  
  // Apply constraint on nodes
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
    constraint(*it, t, g);

  return t + dt;
}


//
// FORCE OBJECTS
//

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator() (NODE n, double t) {
    // Return null force for static extremities
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);
    
    // Sum of gravity and spring forces
    Point force(0, 0, -grav * n.value().mass);
    for (GraphType::incident_iterator iit = n.edge_begin() ; iit != n.edge_end() ; ++iit) {
      Point difference = n.position() - (*iit).node2().position();
      force -= (*iit).value().K * difference / norm_2(difference) * (norm_2(difference) - (*iit).value().L);
    }
    return force;
  }
};

struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator() (NODE n, double t) {
    return Point(0, 0, -grav * n.value().mass);
  }
};

struct MassSpringForce {
  /** Return the spring forces applying to @a n at time @a t. */
  template <typename NODE>
  Point operator() (NODE n, double t) {
    Point mass_spring_force(0, 0, 0);
    for (GraphType::incident_iterator iit = n.edge_begin() ; iit != n.edge_end() ; ++iit) {
      Point difference = n.position() - (*iit).node2().position();
      mass_spring_force -= (*iit).value().K * difference / norm_2(difference) * (norm_2(difference) - (*iit).value().L);
    }
    return mass_spring_force;
  }
};

struct DampingForce {
  /** Initialize the damping constant. */
  double damp_;
  DampingForce(double damp) : damp_(damp) {}
  /** Return the damping force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator() (NODE n, double t) {
    return -damp_ * n.value().vel; 
  }
};


template <typename F1, typename F2>
struct combined_force {
  F1 force1_;
  F2 force2_;
  combined_force(F1 force1, F2 force2)
      : force1_(force1), force2_(force2) {}
  template <typename NODE>
  Point operator() (NODE n, double t) {
    return force1_(n, t) + force2_(n, t);
  }
};

template <typename F1, typename F2>
combined_force<F1, F2> make_combined_force(F1 force1, F2 force2) {
  return combined_force <F1, F2> (force1, force2);
}

template <typename F1, typename F2, typename F3>
combined_force<combined_force<F1, F2>, F3> make_combined_force(F1 force1, F2 force2, F3 force3) {
  return combined_force <combined_force<F1, F2>, F3> (make_combined_force(force1, force2), force3);
}


//
// CONSTRAINT OBJECTS
//

/** Constraint on a vertical plane. */
struct VerticalPlaneConstraint {
  double value_;
  VerticalPlaneConstraint(double value) : value_(value) {}
  template <typename NODE, typename GRAPH>
  void operator() (NODE n, double t, GRAPH& g) {
    if (n.position().z < value_) {
      n.position().z = value_;
      n.value().vel.z = 0;
    }
  }
};

/** Constraint on a sphere. */
struct SphereConstraint {
  Point center_;
  double radius_;
  SphereConstraint(Point center, double radius)
      : center_(center), radius_(radius) {}
  template <typename NODE, typename GRAPH>
  void operator() (NODE n, double t, GRAPH& g) {
    if (norm_2(n.position() - center_) < radius_) {
      Point R = (n.position() - center_) / norm_2(n.position() - center_);
      n.position() = center_ + radius_ * R;
      n.value().vel -= dot(n.value().vel, R) * R;
    }
  }
};

/** Constraint on a sphere and disappear on contact. */
struct DisappearSphereConstraint {
  Point center_;
  double radius_;
  DisappearSphereConstraint(Point center, double radius)
      : center_(center), radius_(radius) {}
  template <typename NODE, typename GRAPH>
  void operator() (NODE n, double t, GRAPH& g) {
    if (norm_2(n.position() - center_) < radius_) {
      g.remove_node(n);
    }
  }
};

template <typename C1, typename C2>
struct combined_constraint {
  C1 constraint1_;
  C2 constraint2_;
  combined_constraint(C1 constraint1, C2 constraint2)
      : constraint1_(constraint1), constraint2_(constraint2) {}
  template <typename NODE, typename GRAPH>
  void operator() (NODE n, double t, GRAPH& g) {
    constraint1_(n, t, g);
    constraint2_(n, t, g);
  }
};

template <typename C1, typename C2>
combined_constraint<C1, C2> make_combined_constraint(C1 constraint1, C2 constraint2) {
  return combined_constraint <C1, C2> (constraint1, constraint2);
}

template <typename C1, typename C2, typename C3>
combined_constraint<combined_constraint<C1, C2>, C3> make_combined_constraint(C1 constraint1, C2 constraint2, C3 constraint3) {
  return combined_constraint <combined_constraint<C1, C2>, C3> (make_combined_constraint(constraint1, constraint2), constraint3);
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
// #if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
// #endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }
  
  // Set initial conditions for your nodes and edges.
  for (auto nit = graph.node_begin() ; nit != graph.node_end() ; ++nit) {
    // Update the mass of the node
    auto n = *nit;
    n.value().mass = (double) 1 / graph.num_nodes();
    // Update the value of the edges
    for (auto iit = n.edge_begin() ; iit != n.edge_end() ; ++iit) {
      auto e = *iit;
      e.value().K = 100;
      e.value().L = norm_2(e.node1().position() - e.node2().position());
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
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        // std::cout << "t = " << t << std::endl;
        
        // Apply forces and constraints
        // symp_euler_step(graph, t, dt, Problem1Force());
        // symp_euler_step(graph, t, dt,
                        // make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1/graph.size())));
        // symp_euler_step(graph, t, dt,
                        // make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1/graph.size())),
                        // make_combined_constraint(VerticalPlaneConstraint(-0.75), SphereConstraint(Point(0.5, 0.5, -0.5), 0.15)));
        symp_euler_step(graph, t, dt,
                        make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1/graph.size())),
                        DisappearSphereConstraint(Point(0.5, 0.5, -0.5), 0.15));
        
        // Clear the viewer's nodes and edges
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

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
