
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
  double K; // spring constant
  double init_length; //initial length of the edge
  EdgeData() : K(100), init_length(1) {}
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
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
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
    //fix the corner points

  }
  constraint(g, t);
 
  return t + dt;
}

/** constraint for the two corner points of the graph
 *  The constructor will find the nodes with coordinate
 *  (0, 0, 0) and (1, 0, 0) and will reset the position
 *  and velocity of these two points everytime this
 *  constraint is called. Since the nodes are stored
 *  as members, the operator() takes constant time.
 */
template<typename GRAPH>
struct CornerConstraint {
  void operator()(GRAPH& g, double t) {
    (void)g; (void)t;
    n1.position() = Point(0, 0, 0);
    n2.position() = Point(1, 0, 0);
    n1.value().vel = Point(0, 0, 0);
    n2.value().vel = Point(0, 0, 0);
  }

  CornerConstraint(GRAPH& g) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      if ((*it).position() == Point(0, 0, 0)) n1 = *it;
      if ((*it).position() == Point(1, 0, 0)) n2 = *it;
    }
  }
private:
  typename GRAPH::node_type n1;
  typename GRAPH::node_type n2;
};

/** Planar constraint z = -0.75
 *  Everytime the operator is called, it will iterate through
 *  the graph's nodes and reset the position and z-component
 *  of the velocity of all nodes if its z coordinate is less than
 *  -0.75.
 */

struct PlanarConstraint {
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).position().z < z_bound){
	(*it).position().z = z_bound;
	(*it).value().vel.z = 0;
      }
    }
  }
private:
  double z_bound = -0.75;
};

/** spherical constraints that makes the graph wrap around a sphere
 *  The operator()(Graph g, double t) will iterate through all the
 *  nodes of the graph and reset the position of any node that are in  *  in the sphere to be on the sphere and set their velocity to be
 *  tangent to be sphere
 */

struct SphereConstraint {
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      if (norm((*it).position() - c) < r){
	(*it).position() += (r - norm( (*it).position() - c)) * ((*it).position() - c);
	(*it).value().vel -= dot((*it).value().vel, (*it).position() - c) / r / r * ((*it).position() - c);
      }
    }
  }
private:
  double r = 0.15;
  Point c = Point(0.5, 0.5, -0.5);
};

/** sphereical constraint that removes nodes from the graph if they 
 *  touch the sphere
 */
struct SphereRemoveConstraint {
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void)t;
    auto it = g.node_begin();
    while(it != g.node_end()){
      if (norm((*it).position() - c) < r){
	it = g.remove_node(it);
      }
      else ++it;
    }
  }
private:
  double r = 0.15;
  Point c = Point(0.5, 0.5, -0.5);
};


/** functor that combines two constraints
 *  Assuming that the constraints in conjunction will not cause
 *  conflict (e.g., if a node satisfies both constraints but the
 *  constraints will set the node's position or velocity differently),
 *  This will apply both constraints to the graph
 */
template<typename C1, typename C2>
struct CombinedConstraint {
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void)t;
    _c1(g, t);
    _c2(g, t);
  }

/** constructor will take in two constraints and store them
 *  @pre The types of both arguments must implement 
      operator()(Graph& g, double t)
 */
  CombinedConstraint (C1 c1, C2 c2)
    : _c1(c1), _c2(c2) {}
private:
  C1 _c1;
  C2 _c2;
};

/** function for creating combined constraints
 *  Have the same specification as the constructor for
 *  CombinedConstraints
 */
template<typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraints(C1 c1, C2 c2){
  return CombinedConstraint<C1, C2>(c1, c2);
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
    (void)t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);
    else {
      // initialize force to be gravity
      Point force = Point(0, 0, -grav * n.value().mass);
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
	force += - (*it).value().K * (n.position() - (*it).node2().position()) / norm(n.position() - (*it).node2().position()) * (norm(n.position() - (*it).node2().position()) - (*it).value().init_length);
      }
      return force;
    }
  }

};

//functor for damping force
struct DampingForce {
public:
  Point operator()(Node n, double t) {
    (void)t;
    return -c * n.value().vel;
  }

  DampingForce (double coef)
    : c(coef) {}

private:
  double c;  
};

//functor for mass spring force on node
struct MassSpringForce {
public:
  Point operator()(Node n, double t) {
    (void)t;
    Point force = Point(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
	force += - (*it).value().K * (n.position() - (*it).node2().position()) / norm(n.position() - (*it).node2().position()) * (norm(n.position() - (*it).node2().position()) - (*it).value().init_length);
      }
   return force;
  }
};

//functor for gravity force
struct GravityForce {
public:
  Point operator()(Node n, double t) {
    (void)t;
    return Point(0, 0, -grav * n.value().mass);
  }
};


/** functor for combining two forces 
 *  Each force object will store two component forces;
 *  the operator () will return the sum of two component forces
 *  The two private force components must be of types that implements operator()(Node n, double t)
 *  Sice CombinedForce has operator()(Node n, double t), components could contain another functor
 *  of type CombinedForce */

template<typename F1, typename F2>
struct CombinedForce {
public:
  virtual Point operator()(Node n, double t) {
    return force1(n, t) + force2(n, t);
  }

  CombinedForce(F1 force_1, F2 force_2)
    : force1(force_1), force2(force_2) {}
private:
  F1 force1;
  F2 force2;
};

/** returns a CombinedForce object containing two component forces
 *  @pre all input tupes must implement Point operator()(Node n, double t)
 */

template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force (F1 force_1, F2 force_2) {
  return CombinedForce<F1, F2>(force_1, force_2);
}

/** overloads make_combined_force(), this takes in three types of forces
 *  @return a CombinedForce object with another CombinedForce object
       containing the second and third force as the second component
 *  @pre all input types must implement Point operator()(Node n, double t)
 */
template <typename F1, typename F2, typename F3>
CombinedForce<F1, CombinedForce<F2, F3>> make_combined_force (F1 force_1, F2 force_2, F3 force_3) {
  return make_combined_force(force_1, make_combined_force(force_2, force_3));
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
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
    (*it).value().mass = 1/double(graph.num_nodes());
    (*it).value().vel = {0, 0, 0};
  }


  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it){
    (*it).value().K = 100;
    (*it).value().init_length = norm((*it).node1().position() - (*it).node2().position());
  }

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it){
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
      double t_end = 5;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1/double(graph.num_nodes()))), make_combined_constraints(CornerConstraint<GraphType>(graph), SphereRemoveConstraint()));

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
