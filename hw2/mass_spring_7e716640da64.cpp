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
  double len;
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
double symp_euler_step(G& g, double t, double dt, F force) {
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

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C cons) {
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
  cons(g);
  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  Point operator()(Node n, double t) {
    // HW2 #1: YOUR CODE HERE
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
    (void) t;
    double K = 100;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    else{
      Point spring_force = Point(0,0,0);
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
        auto no = *it;
        Point adj = no.node2().position();
        double norm_distance = norm(n.position() - adj);
        spring_force += -K*(n.position() - adj) / norm_distance *(norm_distance - no.value().len);

      }
      return spring_force + n.value().mass * Point(0,0,-1.0*grav);

    }
    return Point(0);
  }
};
// gravity force
struct GravityForce{
  Point operator() (Node n, double t) {
      (void) t;
      return n.value().mass * Point(0,0,-1.0*grav);
  }

};
//spring force
//@ param @K spring constant
struct MassSprintForce {
  double K_;
  MassSprintForce(double K) : K_(K){

  }
  Point operator() (Node n, double t) {
    (void) t;
    Point spring_force = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      auto no = *it;
      Point adj = no.node2().position();
      double norm_distance = norm(n.position() - adj);
      spring_force += -K_*(n.position() - adj) / norm_distance *(norm_distance - no.value().len);

      }
    return spring_force;


  }

};

// dampling force
//@param @a c dampling constant
struct DamplingForce {
  double c_;
  DamplingForce(double c) : c_(c){}
  Point operator() (Node n, double t){
    (void) t;
    return -c_*n.value().vel;
  }

};


template<typename F1, typename F2>
struct CombineForce {
  F1 f1_;
  F2 f2_;
  CombineForce(F1 f1, F2 f2): f1_(f1), f2_(f2){

  }
  Point operator() (Node n, double t) {
    return f1_(n,t)+ f2_(n,t);
  }

};

template<typename F1, typename F2>
CombineForce<F1, F2> make_combined_force(F1 f1, F2 f2){
  return CombineForce<F1 ,F2>(f1, f2);
}

template<typename F1, typename F2, typename F3>
CombineForce<CombineForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3){
  return CombineForce<CombineForce<F1, F2>, F3>(CombineForce<F1 ,F2>(f1, f2), f3);
}


// constraint for a plane
struct PlaneConstraint {
    double c = -0.75;
    void operator()(GraphType &g){
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        auto n = *it;
        if (n.position().z < c){
          n.position().z = c;
          n.value().vel.z = 0;
        }
      }
    }

};

// constraint for a sphere
struct SphereConstraint {
  double r = 0.15;
  Point c = Point(0.5, 0.5, -0.5);
  void operator()(GraphType &g){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      auto p = n.position();
      double dis = norm(p-c);
      if ( dis < r){
        Point R = (p-c) / dis;
        n.position() = c + r*R;
        n.value().vel -= dot(n.value().vel, R)*R;
      }
    }
  }
};

struct SphereRemove {
  double r = 0.15;
  Point c = Point(0.5, 0.5, -0.5);
  void operator()(GraphType &g){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      double dis = norm(n.position()-c);
      if (dis < r){
        g.remove_node(n);
      }
    }
  }
};

// combine contraints

template<typename C1, typename C2>
struct CombineConstraint {
  C1 c1_;
  C2 c2_;
  CombineConstraint(C1 c1, C2 c2): c1_(c1), c2_(c2){

  }
  Point operator() (GraphType &g) {
    c1_(g); c2_(g);
  }

};

/* 
// @param two contraints @a c1, @a c2
// @return combined constraints
*/
template<typename C1, typename C2>
CombineConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2){
  return CombineConstraint<C1 ,C2>(c1, c2);
}

template<typename C1, typename C2, typename C3>
CombineConstraint<CombineConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3){
  return CombineConstraint<CombineConstraint<C1, C2>, C3>(CombineConstraint<C1 ,C2>(c1, c2), c3);
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
#if 1
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  // set length
  for (GraphType::EdgeIterator it = graph.edge_begin(); it != graph.edge_end(); ++it){
    (*it).value().len = (*it).length();
  }

  // set mass
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    (*it).value().vel = Point(0,0,0);
    (*it).value().mass = 1.0 / graph.num_nodes();
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
        symp_euler_step(graph, t, dt, Problem1Force());
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSprintForce(100)));
        //auto com_force = make_combined_force(GravityForce(), MassSprintForce(100),DamplingForce(1.0/graph.num_nodes()));
        //symp_euler_step(graph, t, dt, com_force);

        //auto com_cons = make_combined_constraint(PlaneConstraint(), SphereConstraint());
        //symp_euler_step(graph, t, dt, com_force, com_cons);

        //viewer.clear();
        //node_map.clear();

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
