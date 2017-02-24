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

struct EdgeData {
  double length;
};

struct Problem1Constraint;

double init_length;

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
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
template <typename G, typename F, typename C = Problem1Constraint>
double symp_euler_step(G& g, double t, double dt, F force, C constraint = C() ) {
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
  constraint(g, t);

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
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    double K = 100;
    Point springForce = Point(0,0,0);
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) return springForce;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Point node2 = (*it).node2().position();
      Point diff = n.position() - node2;
      springForce += diff/norm(diff) * (norm(diff)- (*it).value().length);
    }
    return -1.0 * K * springForce + n.value().mass * Point(0,0,-1.0 * grav);
  }
};

struct GravityForce {
  double grav_;
  GravityForce(double g = grav) : grav_(g) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n; (void) t;     // silence compiler warnings
    return n.value().mass * Point(0,0,-1.0 * grav_);
  }
};

struct MassSpringForce {
  double K_;
  MassSpringForce(double K = 100) : K_(K) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    Point springForce = Point(0,0,0);
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) return springForce;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Point node2 = (*it).node2().position();
      Point diff = n.position() - node2;
      springForce += diff/norm(diff) * (norm(diff)- (*it).value().length);
    }
    return -1.0 * K_ * springForce;
  }
};

struct DampingForce {
  double c_;
  DampingForce(double c) : c_(c) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n; (void) t;     // silence compiler warnings
    return -1.0 * c_ * n.value().vel;
  }
};

template <typename Force1, typename Force2>
struct CombinedForces {
  Force1 f1_;
  Force2 f2_;
  CombinedForces(Force1 f1, Force2 f2) : f1_(f1), f2_(f2) {}

  Point operator() (Node n, double t){
    (void) n; (void) t;     // silence compiler warnings
    return f1_(n,t) + f2_(n,t);
  }
};

/** Combine two Forces
 * @param[in] f1 f2 Forces to be combined
 * @pre Forces @a f1, @a f2 are valid Force functions that take two inputs: a node n and time t as f(n,t)
 * @return CombinedForces object represents the combination of forces @a f1, @a f2.
 */
template <typename Force1, typename Force2>
CombinedForces<Force1, Force2> make_combined_force(Force1 f1, Force2 f2){
  return CombinedForces<Force1, Force2> ({f1, f2});
}

/** Combine three Forces
* @param[in] f1 f2 f3 Forces to be combined
* @pre Forces @a f1, @a f2, @a f3 are valid Force functions that take two inputs: a node n and time t as f(n,t)
* @return CombinedForces object represents the combination of forces @a f1, @a f2, @a f3.
*/
template <typename Force1, typename Force2, typename Force3>
CombinedForces<CombinedForces<Force1, Force2>, Force3> make_combined_force(Force1 f1, Force2 f2, Force3 f3){
  return CombinedForces<CombinedForces<Force1, Force2>, Force3> (CombinedForces<Force1, Force2>(f1,f2), f3);
}

struct Problem1Constraint {
  void operator()(GraphType& graph, double t) {
    (void) graph; (void) t;     // silence compiler warnings
    for(auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      if((*it).position() == Point(0,0,0) || (*it).position() == Point(1,0,0)) {
        (*it).value().vel =  Point(0,0,0);
      }
    }
  }
};

struct PlaneConstraint {
  double z_; //coordinate for the horizontal plane
  PlaneConstraint(double z) : z_(z) {}

  void operator()(GraphType& graph, double t) {
    (void) graph; (void) t;     // silence compiler warnings
    for(auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      if((*it).position().z < z_) {
        (*it).position().z = z_;
        (*it).value().vel.z = 0;
      }
    }
  }
};

struct SphereConstraint {
  Point center_;
  double radius_;
  SphereConstraint(double r, Point c) : center_(c), radius_(r) {}

  void operator()(GraphType& graph, double t) {
    (void) graph; (void) t;     // silence compiler warnings
    for(auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      double distance = norm((*it).position() - center_);
      if(distance < radius_){
        (*it).position() = center_ + (radius_ /distance) * ((*it).position() - center_) ;
        Point R = ((*it).position() - center_)/distance;
        (*it).value().vel -= dot((*it).value().vel, R) * R;
      }
    }
  }
};

struct RemoveConstraint {
  Point center_;
  double radius_;
  RemoveConstraint(double r, Point c) : center_(c), radius_(r) {}

  void operator()(GraphType& graph, double t) {
    (void) graph; (void) t;     // silence compiler warnings
    for(auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      double distance = norm((*it).position() - center_);
      if(distance < radius_){
        graph.remove_node(*it);
      }
    }
  }
};

template <typename Constraint1, typename Constraint2>
struct CombinedConstraints {
  Constraint1 c1_;
  Constraint2 c2_;
  CombinedConstraints(Constraint1 c1, Constraint2 c2): c1_(c1), c2_(c2) {}

  void operator() (GraphType& g, double t){
    c1_(g,t);
    c2_(g,t);
  }
};

template <typename Constraint1, typename Constraint2>
CombinedConstraints<Constraint1, Constraint2> make_combined_constraints(Constraint1 c1, Constraint2 c2){
  return CombinedConstraints<Constraint1, Constraint2> ({c1, c2});
}

template <typename Constraint1, typename Constraint2, typename Constraint3>
CombinedConstraints<CombinedConstraints<Constraint1, Constraint2>, Constraint3> make_combined_constraints(Constraint1 c1, Constraint2 c2, Constraint3 c3){
  return CombinedConstraints<CombinedConstraints<Constraint1, Constraint2>, Constraint3> ({{c1, c2}, c3});
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
//#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().vel = Point(0,0,0);
    n.value().mass = 1.0/graph.num_nodes();
  }

  init_length = (*(graph.edge_begin())).length();

  for (auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei ) {
    (*ei).value().length = (*ei).length();
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

    GravityForce gravity(grav);
    MassSpringForce massSpringForce(100);
    DampingForce damping(1.0/graph.num_nodes());
    //Final Force
    auto f = make_combined_force(gravity, massSpringForce, damping);

    //Setting Constraints
    PlaneConstraint plane(-0.75);
    SphereConstraint sphereConstraint(0.15, Point(0.5, 0.5, -0.5));
    RemoveConstraint removeSphere(0.15, Point(0.5, 0.5, -0.5));

    for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
      //std::cout << "t = " << t << std::endl;
      symp_euler_step(graph, t, dt, f, make_combined_constraints(plane, removeSphere, Problem1Constraint()));

      viewer.clear();
      node_map.clear();

      viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
      viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

      // Update viewer with nodes' new positions
      viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
