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
  NodeData(Point vel_, double mass_) : vel(vel_), mass(mass_) {}
};
/** Custom structure of data to store with Edges */
struct EdgeData {
  double nat_length;
  double constant;
  EdgeData(): nat_length(0), constant(0) {}
  EdgeData(double nat_, double constant_) : nat_length(nat_), constant(constant_) {}
};

/** Constraint functor
@param normal: normal vector to PlaneConstraint
@param level: height of PlaneConstraint*/
struct PlaneConstraint{
  Point normal;
  double level;
  PlaneConstraint(Point n_, double l_) : normal(n_), level(l_) {}
  template<typename V, typename E>
  void operator()(Graph<V,E> &g, double t){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if (inner_prod(n.position(),normal)<level){
        n.position() = n.position() - (inner_prod(n.position(), normal) - level) * normal;
        n.value().vel = n.value().vel - inner_prod(n.value().vel, normal) * normal;
      }
    }
  }
};

/** Constraint functor
@param center: center of sphere
@param radius: radius of sphere*/
struct SphereConstraint{
  Point center;
  double radius;
  SphereConstraint(Point c_, double r_) : center(c_), radius(r_) {}
  template<typename V, typename E>
  void operator()(Graph<V,E> &g, double t){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      Point h = n.position() - center;
      double length = norm(h);
      if (length<radius){
        Point normal =  h/length;
        n.position() = center + normal * radius;
        n.value().vel = n.value().vel - (inner_prod(n.value().vel, normal)) * normal;
      }
    }
  }
};

/** Constraint functor
@param center: center of sphere
@param radius: radius of sphere*/
struct SphereDisappear{
  Point center;
  double radius;
  SphereDisappear(Point c_, double r_) : center(c_), radius(r_) {}
  template<typename V, typename E>
  void operator()(Graph<V,E> &g, double t){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      Point h = n.position() - center;
      double length = norm(h);
      if (length<radius){
        g.remove_node(n);
      }
    }
  }
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
    if (n.position()!=Point(0,0,0) && n.position() != Point(1,0,0) ){
    n.position() += n.value().vel * dt;}
  }

  //PlaneConstraint constraint = PlaneConstraint(Point(0,0,1), -0.75);
  //SphereConstraint constraint = SphereConstraint(Point(0.5, 0.5, -0.5), 0.15);
  SphereDisappear constraint = SphereDisappear(Point(0.5, 0.5, -0.5), 0.15);
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position()!=Point(0,0,0) && n.position() != Point(1,0,0) ){
    n.value().vel += force(n, t) * (dt / n.value().mass);}
  }

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
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    } else {
      //0.25 0.04166666667 0.02040816327 0.0101010101
      return spring_force(n) + gravity_force(n);
    }
  }
  template <typename NODE>
  Point spring_force(NODE n) {
    Point force = Point(0,0,0);
    Point node_position = n.position();
    // HW2 #1: YOUR CODE HERE
    for (auto it = n.edge_begin(); it!=n.edge_end(); ++it){
      auto n = *it;
      Point distance = node_position - n.node2().position();
      force+=-n.value().constant*(distance/norm(distance))*(norm(distance)-n.value().nat_length);
    }
    return force;
  }
  template <typename NODE>
  Point gravity_force(NODE n) {
    return Point(0,0,-grav * n.value().mass);
  }
};

/* Gravity force class
  @param g: Gravity constant
*/
struct GravityForce{
  GravityForce(double g_) : g(g_) {}
  double g;
  template <typename NODE>
  Point operator()(NODE n){
    return Point(0,0,-g * n.value().mass);
  }
};

/* Spring force class
*/
struct MassSpringForce{
  template <typename NODE>
  Point operator()(NODE n){
    Point force = Point(0,0,0);
    Point node_position = n.position();
    // HW2 #1: YOUR CODE HERE
    for (auto it = n.edge_begin(); it!=n.edge_end(); ++it){
      auto n = *it;
      Point distance = node_position - n.node2().position();
      force+=-n.value().constant*(distance/norm(distance))*(norm(distance)-n.value().nat_length);
    }
    return force;
  }
};

/* Gravity force class
  @param g: Damping constant
*/
struct DampingForce{
  DampingForce(double c_) : c(c_) {}
  double c;
  template <typename NODE>
  Point operator()(NODE n){
    return -c* n.value().vel;
  }
};
/* A zero force*/
struct ZeroForce{
  template <typename NODE>
  Point operator()(NODE n){
    return Point(0,0,0);
  }
};

/* A  class to add forces up*/
template<typename F1, typename F2, typename F3>
struct combined_force{
  F1 force1;
  F2 force2;
  F3 force3;
  combined_force(F1 f1, F2 f2, F3 f3) : force1(f1), force2(f2), force3(f3) {}
  template<typename NODE>
  Point operator()(NODE n, double t){
    return force1(n) + force2(n) + force3(n);
  }
};

/*This function adds two forces
@param f1: functor that evaluates to the force at a given node
@param f1: functor that evaluates to the force at a given node
@return A functor that sums the forces at a given node
*/
template<typename F1, typename F2>
combined_force<F1, F2, ZeroForce> make_combined_force(F1 f1, F2 f2){
  return combined_force<F1, F2, ZeroForce> (f1, f2, ZeroForce());
}
/*This function adds three forces
@param f1: functor that evaluates to the force at a given node
@param f2: functor that evaluates to the force at a given node
@param f3: functor that evaluates to the force at a given node
@return A functor that sums the forces at a given node
*/
template<typename F1, typename F2, typename F3>
combined_force<F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3){
  return combined_force<F1, F2, F3> (f1, f2, f3);
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

    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.set_value(NodeData(Point(0,0,0),(1.0/((double)graph.size()))));
  }

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    auto n = *it;
    n.set_value(EdgeData(n.length(),100));
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;
  std::cout << graph.edge(0).length() <<std::endl;

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
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(grav), MassSpringForce(),DampingForce(1/(double)graph.size())));
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
