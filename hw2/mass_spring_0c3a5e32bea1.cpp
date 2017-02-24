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
#include <limits>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// Damping constant
static double c = 0;

// Variables for HW2 problem 1
// static constexpr double K = 100;  
// static double L = 0;


// Set the position of nodes that will remain constant
static Point const_p1 = Point(0,0,0);
static Point const_p2 = Point(1,0,0);

// Unique id's of nodes that will remain constant (for HW2 problem 4)
static unsigned const_uid1 = 0;
static unsigned const_uid2 = 0;

/** Custom structures for Mass, Force, and Velocity */
struct Mass{
  double val_;
  Mass(double val) : val_(val){}
};
struct Force{
  Point val_;
  Force(double a) {
    val_ = Point(a,a,a);
  }
  Force(Point val): val_(val){}
};
struct Velocity{
  Point val_;
  Velocity(double a) {
    val_ = Point(a,a,a);
  }
  Velocity(Point val): val_(val){}
};

/*Some arithmetic operators for objects of class Force */

Force operator*(Force f, double a){
  return Force(a*f.val_);
}

Force operator+(Point a, Force b) {
  return Force(a + b.val_);
}
Force operator+(Force a, Force b) {
  return Force(a.val_ + b.val_);
}

/** Custom structure of data to store with Nodes */
struct NodeData {
  Velocity vel;       //< Node velocity
  Mass mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;        //< Spring constant
  double L;        //< Rest-length
  EdgeData() : K(100), L(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


// 
// CONSTRAINT FUNCTION OBJECTS:
//

/* Constraint function objects find nodes in the Graph that violate the 
 * constraint (after the position update) and reset this nodes' positions and
 * velocities (before the forces are calculated). 
 */

// Add constant node constraints
struct ConstantConstraint{
  template <typename GRAPH>
  void operator()(GRAPH& g, double t){
    (void) t; // To silence compiler warning
    g.node_uid(const_uid1).position() = const_p1;
    g.node_uid(const_uid2).position() = const_p2;
    g.node_uid(const_uid1).value().vel.val_ = Point(0,0,0);
    g.node_uid(const_uid2).value().vel.val_ = Point(0,0,0);
  }
};

// Add a constraint for a plain with set properties
struct PlaneConstraint{
  template <typename GRAPH>
  void operator()(GRAPH& g, double t){
    (void) t; // To silence compiler warning
    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      if((*ni).position().z < -0.75){
        // Set this node's position to the nearest point on the plane
        (*ni).position().z = -0.75; 
        // Set z component of Node's velocity to zero
        (*ni).value().vel.val_.z = 0; 
      }
    }
  }
};

// Add a constraint for a sphere with set properties.
// Fix by setting position to the nearest point on the sphere.
struct SphereConstraint{
  template <typename GRAPH>
  void operator()(GRAPH& g, double t){
    (void) t; // To silence compiler warning
    Point  c = Point(0.5, 0.5, -0.5);
    double r = 0.15; 
    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      Point p = (*ni).position()-c;
      if(norm(p) < r){
        // Set position to nearest point on the surface of the sphere.
        Point Ri = (p/norm(p));
        (*ni).position() += (r-norm(p))*Ri; 
        (*ni).value().vel.val_ -= dot((*ni).value().vel.val_, Ri)*Ri; 
      }
    }
  }
};

// Add a constraint for a sphere with set properties.
// Fix by removing the node and all of its edges.
struct SphereConstraintR{
  template <typename GRAPH>
  void operator()(GRAPH& g, double t){
    (void) t; // To silence compiler warning
    Point  c = Point(0.5, 0.5, -0.5);
    double r = 0.15; 
    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      Point p = (*ni).position()-c;
      if(norm(p) < r){
        // Set position to nearest point on the surface of the sphere.
        g.remove_node(*ni);
      }
    }
  }
};


/* 
 * Make combined constraint function objects
 */

/* Return type of combine constraint function objects. */
template<typename C1, typename C2> 
struct CombinedConstraint{
  C1 c1_;
  C2 c2_;
  template<typename GRAPH>
  void operator()(GRAPH& g, double t){
    c1_(g,t);
    c2_(g,t);
  }
  CombinedConstraint(C1 c1, C2 c2): c1_(c1), c2_(c2){};
};

// Combination constraint objects that combine two or three constraints together

// Takes two components as arguments
template <typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  return CombinedConstraint<C1, C2>(c1, c2);
};

// Takes three components as arguments
template <typename C1, typename C2, typename C3>
CombinedConstraint<C1, CombinedConstraint<C2, C3>> make_combined_constraint(
                                                        C1 c1, C2 c2, C3 c3){
  return CombinedConstraint<C1, CombinedConstraint<C2,C3>>(c1,
                                CombinedConstraint<C2,C3> (c2,c3));
};

//
// SYMPLECTIC EULER STEP
//

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
    n.position() += n.value().vel.val_ * dt;
    /* Keep specific nodes fixed for HW2 problem 1 */
    // if(n.position()!=Point(0,0,0) && n.position()!= Point(1,0,0)){
    //   n.position() += n.value().vel.val_ * dt;
    // }
  }
  
  /* Fix nodes that violated a constraint for HW2 problem 4 */
  ConstantConstraint constraint1;
  PlaneConstraint  constraint2;
  //SphereConstraint constraint3;
  SphereConstraintR constraint4;
  
  //constraint1(g,t);
  //constraint2(g,t);
  //constraint3(g,t);
  //constraint4(g,t);

   /*auto constraint = make_combined_constraint(constraint2, 
                                                constraint4); */

   auto constraint = make_combined_constraint(constraint1, 
                                              constraint2, 
                                              constraint4); 
   constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel.val_ += force(n, t) * (dt / n.value().mass.val_);
    /* Skip update step at constant nodes for HW2 problem 3 */
    // if(n.position()!=Point(0,0,0) && n.position()!= Point(1,0,0)){
    //   n.value().vel.val_ += force(n, t) * (dt / n.value().mass.val_);
    // }
  }

  return t + dt;
}

//
// FORCE FUNCTION OBJECTS
//


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; // To silence compiler warning
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    else{
      Point f_spring = Point(0,0,0);
      Point f_grav = (n.value().mass.val_)*Point(0,0,-grav);
      Point delta;
      double L, K;
      for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii){
        delta = n.position() - (*ii).node2().position();
        L = (*ii).value().L;
        K = (*ii).value().K;
        f_spring -= K*(delta/norm(delta))*(norm(delta) - L);
      }
      return f_spring + f_grav;
    }
  }
};

/** Force of gravity. */
struct GravityForce {
  template <typename NODE>
  Force operator()(NODE n, double t){
    (void) t; // To silence compiler warning
    return Force((n.value().mass.val_)*Point(0,0,-grav));
  }
};

/** Class of spring forces. */
struct MassSpringForce {
  template <typename NODE>
  Force operator()(NODE n, double t){
    (void) t; // To silence compiler warning
    Point p_spring = Point(0,0,0); // Initial spring 
    Point delta; // Distance between nodes
    double L, K;
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii){
      delta = n.position() - (*ii).node2().position();
      L = (*ii).value().L;
      K = (*ii).value().K;
      p_spring -= K*(delta/norm(delta))*(norm(delta) - L);
    }
    return Force(p_spring);
  }
};

/** Damping force. */
struct DampingForce {
  template <typename NODE>
  Force operator()(NODE n, double t){
    (void) t; // To silence compiler warning
    return Force(-c*(n.value().vel.val_));
  }
};

/* 
 * Make combined forces function objects
 */

/* Return type of combination force function objects. */
template<typename F1, typename F2> 
struct CombinedForce{
  F1 f1_;
  F2 f2_;
  template<typename NODE>
  Point operator()(NODE n, double t){
    return (f1_(n, t) + f2_(n, t)).val_;
  }
  CombinedForce(F1 f1, F2 f2): f1_(f1), f2_(f2){};
};

/* Combination force objects that find their's component forces' values and 
 * returns their sums.
 */

// Takes two components as arguments
template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return CombinedForce<F1, F2>(f1, f2);
};

// Takes three components as arguments
template <typename F1, typename F2, typename F3>
CombinedForce<F1,CombinedForce<F2, F3>>make_combined_force(F1 f1, F2 f2, F3 f3){
  return CombinedForce<F1,CombinedForce<F2, F3>>(f1,CombinedForce<F2,F3>(f2,f3));
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

  // Calculate the initial length of edges. (For HW2 Problem 1)
  // L = (*graph.edge_begin()).length();

  // Calculate the damping constant. (For HW2 Problem 3)
  c = (double)1/graph.size();

  // Set initial conditions for nodes and edges, if necessary; and
  // find the unique id's of nodes that will remain constant.
  NodeData initial_n_data;
  initial_n_data.vel.val_  = Point(0,0,0);
  initial_n_data.mass.val_ = (double)1/graph.size();

  EdgeData initial_e_data;
  initial_e_data.K = 100;

  for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
    // Set initial values for nodes
    (*ni).value() = initial_n_data;
    if((*ni).position() == const_p1){
      const_uid1 = (*ni).uid(); // For constant constraint
    }
    if((*ni).position() == const_p2){
      const_uid2 = (*ni).uid(); // For constant constraint
    }
    // Iterate through incident edges and set initial values
    for(auto ii = (*ni).edge_begin(); ii != (*ni).edge_end(); ++ii){
      initial_e_data.L = (*ii).length();
      (*ii).value() = initial_e_data; 
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
        //std::cout << "t = " << t << std::endl;
        //symp_euler_step(graph, t, dt, Problem1Force()); //For HW2 Problem 1
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), 
                                                          MassSpringForce(),
                                                          DampingForce()));
        // Clear the viewer's nodes and edges
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
