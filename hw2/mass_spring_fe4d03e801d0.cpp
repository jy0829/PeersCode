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
  NodeData() : vel(0), mass(1) {}
};
/** Custom structure of data to store with Edges */
struct EdgeData {
  double K = 100;
  double length;
  EdgeData(double k, double l) : K(k), length(l) {}
  EdgeData() {}
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C c) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) continue;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    
  }
  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) continue;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    
  }

  c(g, t);

  return t + dt;
}


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
    // if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      // return Point(0, 0, 0);
    
    Point fGrav = Point(0, 0, (-1) * n.value().mass * grav);
    Point fSpring = Point(0, 0, 0);
    // int i = 0;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      // i++;
      double L = (*it).value().length;
      double K = (*it).value().K;
      Node n2;
      if (n == (*it).node1()) {
        n2 = (*it).node2();
      } else n2 = (*it).node1();
      double distance = norm(n.position() - n2.position());
      fSpring = fSpring + (-1) * K * (distance - L) * (n.position() - n2.position()) / distance;
      // std::cout << "t = " << t << ", i = " << i << ": (*it).node2() " << (*it).node2().position() << std::endl;
    }
    (void) t;   // silence compiler warnings

    return fGrav + fSpring;
    // return Point(0);
  }
};

/** Gravity Force function object. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point fGrav = Point(0, 0, (-1) * n.value().mass * grav);
    (void) t;
    return fGrav;    
  }
};
/** Mass Spring Force function object. */
struct MassSpringForce {
  MassSpringForce(){}
  template<typename NODE>
  Point operator()(NODE n, double t) {

    Point fSpring = Point(0, 0, 0);
    // int i = 0;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      // i++;
      double L = (*it).value().length;
      double K = (*it).value().K;
      Node n2;
      if (n == (*it).node1()) {
        n2 = (*it).node2();
      } else n2 = (*it).node1();
      double distance = norm(n.position() - n2.position());
      fSpring = fSpring + (-1) * K * (distance - L) * (n.position() - n2.position()) / distance;
      // std::cout << "t = " << t << ", i = " << i << ": (*it).node2() " << (*it).node2().position() << std::endl;
    }
    (void) t;
    return fSpring;    
  }
};
/** Damping Force function object. */
struct DampingForce {
  DampingForce(double c):c_(c) {}
  template<typename NODE>
  Point operator()(NODE n, double t) {
    Point fDamping = (-1) * c_ *n.value().vel;
    (void) t;
    return fDamping;
  }
  double c_ = 0.0;  
};

template <typename F1, typename F2>
struct TwoForce {
  TwoForce(F1 f1, F2 f2):f1_(f1), f2_(f2){}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return f1_(n, t) + f2_(n, t);
  }
  F1 f1_;
  F2 f2_;
};

template <typename F1, typename F2>
TwoForce<F1, F2> make_combined_force(F1 f1, F2 f2) {

  return TwoForce<F1, F2>(f1, f2);
}

template <typename F1, typename F2, typename F3>
TwoForce<TwoForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  TwoForce<F1, F2> twoForce = make_combined_force(f1, f2);
  return TwoForce<TwoForce<F1, F2>, F3>(twoForce, f3);
}

struct PlaneConstraint {
  PlaneConstraint(double z) : z_(z) {}
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if (dot((*it).position(), Point(0, 0, 1)) < z_) {
        (*it).position().elem[2] = z_;
        (*it).value().vel.elem[2] = 0;
      }
    }
    (void) t;
  }
  double z_;
};

struct SphereConstraint {
  SphereConstraint(Point center, double radius) : center_(center), radius_(radius) {}
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if (norm((*it).position() - center_) < radius_) {
        Point R = ((*it).position() - center_) / norm((*it).position() - center_);
        (*it).position() =  center_ + R * radius_;
        (*it).value().vel = (*it).value().vel - dot((*it).value().vel, R) * R;
      }
    }
    (void) t;
  }
  Point center_;
  double radius_;
};

struct RemoveSphereConstraint {
  RemoveSphereConstraint(Point center, double radius) : center_(center), radius_(radius) {}
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end();) {
      if (norm((*it).position() - center_) < radius_) {
        g.remove_node(*it);
      } else{
        ++it;
      }
    }
    (void) t;
  }


  Point center_;
  double radius_;  
};

struct ConstantConstraint {
  ConstantConstraint() {}
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).position() == Point(0, 0, 0) || (*it).position() == Point(1, 0, 0))
        (*it).value().vel = Point(0, 0, 0);
    }
    (void) t;
  }
};

template <typename C1, typename C2>
struct TwoConstraint {
  TwoConstraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    c1_(g, t);
    c2_(g, t);

  } 
  C1 c1_;
  C2 c2_;
};

template <typename C1, typename C2>
TwoConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  return TwoConstraint<C1, C2>(c1, c2);
}

template <typename C1, typename C2, typename C3>
TwoConstraint<TwoConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  TwoConstraint<C1, C2> twoConstraint = make_combined_constraint(c1, c2);
  return TwoConstraint<TwoConstraint<C1, C2>, C3>(twoConstraint, c3);
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
  while (CME212::getline_parsed(nodes_file, p)) {
    nodes.push_back(graph.add_node(p));
    // std::cout << "idx: " << nodes[nodes.size() - 1].index() << std::endl;
  }


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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().vel = Point();
    (*it).value().mass = 1.0/graph.num_nodes();
  }

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    // (*it).value().K = 100;
    (*it).value().length = (*it).length();
    // std::cout << "length is: " << (*it).length() << std::endl;
    // (*it).otherSide().value().length = (*it).length();
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the Viewer
  CME212::SFML_Viewer viewer;
  auto node_map = viewer.empty_node_map(graph);

  // for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
  //   std::cout << "node1: " << (*it).node1().index() << "; node2: " << (*it).node2().index() << std::endl;
  // }

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
        // symp_euler_step(graph, t, dt, Problem1Force(), PlaneConstraint(-0.75));
        // symp_euler_step(graph, t, dt, Problem1Force(), SphereConstraint(Point(0.5, 0.5, -0.5), 0.15));
        // symp_euler_step(graph, t, dt, Problem1Force(), RemoveSphereConstraint(Point(0.5, 0.5, -0.5), 0.15));
        auto totalForce = make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes()));
        symp_euler_step(graph, t, dt, totalForce, make_combined_constraint(ConstantConstraint(), RemoveSphereConstraint(Point(0.5, 0.5, -0.5), 0.15), PlaneConstraint(-0.75)));

        // Update viewer with nodes' new positions
        viewer.clear();
        node_map.clear();
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
