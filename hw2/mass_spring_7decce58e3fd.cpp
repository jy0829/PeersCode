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
 static constexpr double K = 100;
 static double L;
static double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;   //< Node velocity
  double mass; //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double L;
  double K;

  EdgeData() {
    L = 1.0;
    K = 100.0;
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
double symp_euler_step(G &g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // if ((n.position()==Point(0,0,0))||(n.position()==Point(1,0,0)))
    // continue;
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

/** Structure for calculating gravity force*/
struct GravityForce {

  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    Point force_grav = n.value().mass * Point(0, 0, -grav);
    return force_grav;
  }
};
/** Structure for calculating mass spring force*/
struct MassSpringForce {

  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    Point force_spring = Point(0, 0, 0);
    auto xi = n.position();
    
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      auto xj = (*iter).node2().position();
      double length = (*iter).value().L;
      force_spring = force_spring -
                     (*iter).value().K  * (xi - xj)/
                         norm(xi - xj)* (norm(xi - xj) - length) ;
    
    }
    
    return force_spring;
  }
};
/** Structure for calculating the Damping force*/
struct DampingForce {
  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    Point damp_force = -c * n.value().vel;
    return damp_force;
  }
};
/** Structure used for combining the forces*/
template <typename force1, typename force2> struct CombinedForce {
  force1 f1_;
  force2 f2_;

  template <typename NODE> Point operator()(NODE n, double t) {
    (void)t;
    return f1_(n, t) + f2_(n, t);
  }
};

/** Function for combining two forces*/
template <typename force1, typename force2, typename force3>
CombinedForce<CombinedForce<force1, force2>, force3>
make_combined_force(force1 f1, force2 f2, force3 f3) {
  return {{f1, f2}, f3};
}

/** Function for combining three forces*/
template <typename force1, typename force2>
CombinedForce<force1, force2> make_combined_force(force1 f1, force2 f2) {
  return {f1, f2};
}
/** Force function object for HW2 #1. */
struct Problem1Force {

 /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
     (void) t;  // silence compiler warnings
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    Point force_grav=n.value().mass*Point(0,0,-grav);
    Point force_spring=Point(0,0,0);
    auto xi=n.position();

    for (auto iter=n.edge_begin();iter!=n.edge_end();++iter){
      auto xj=(*iter).node2().position();
      double length=(*iter).value();
      force_spring=force_spring-K*(xi-xj)/norm(xi-xj)*(norm(xi-xj)-length);
    }
    return force_grav+force_spring;
  }
};

/** Struct for ensuring the fixed nodes*/
template <typename NODE> struct FixedConstraint {
  Point p;
  NODE n;
  template <typename Graph> void operator()(Graph &g, double t) {
    (void)t;
    g.node(n.index()).position() = p;
    g.node(n.index()).value().vel = Point(0, 0, 0);
  }
};

/** Struct for ensuring the Plane Constraint*/
struct PlaneConstraint {
  template <typename Graph> void operator()(Graph &g, double t) {
    (void)t;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      auto n = *iter;
      if (n.position().z < -0.75) {
        n.value().vel.z = 0.0;
        n.position().z = -0.75;
      }
    }
  }
};

/** Struct for ensuring the Sphere constraint with node remove */
struct SphereConstraint2 {
  template <typename Graph> void operator()(Graph &g, double t) {
    (void)t;
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      auto n = *iter;

      if (norm(n.position() - c) < r) {
        //std::cout<<n.index()<<" removed node\n";
        g.remove_node(n);
      }
    }
  }
};

/**Struct for ensuring the Sphere Constraint */
struct SphereConstraint {
  template <typename Graph> void operator()(Graph &g, double t) {
    (void)t;
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      auto n = *iter;

      if (norm(n.position() - c) < r) {
        Point unit_normal = (n.position() - c) / norm(n.position() - c);
        n.position() = r * unit_normal + c;
        n.value().vel =
            n.value().vel - (dot(n.value().vel, unit_normal)) * unit_normal;
      }
    }
  }
};

/** Struct for combining two constraints*/
template <typename C1, typename C2> struct CombinedConstraint {
  C1 c1;
  C2 c2;
  template <typename Graph> void operator()(Graph &g, double t) {
    
    c1(g, t);
    c2(g, t);
  }
};

/** function for combining two constraints */
template <typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraints(C1 c1, C2 c2) {
  return {c1,c2};
}

/** function for combining three constraints */
template <typename C1, typename C2, typename C3>
CombinedConstraint<CombinedConstraint<C1, C2>, C3>
make_combined_constraints(C1 c1, C2 c2, C3 c3) {
  return {{c1,c2},c3};
}

/** function for combining four constraints */
template <typename C1, typename C2, typename C3, typename C4>
CombinedConstraint<CombinedConstraint<CombinedConstraint<C1, C2>, C3>, C4>
make_combined_constraints(C1 c1, C2 c2, C3 c3, C4 c4) {
  return {{{c1,c2},c3},c4};
}

/** function for combining five constraints */
template <typename C1, typename C2, typename C3, typename C4, typename C5>
CombinedConstraint<
    CombinedConstraint<CombinedConstraint<CombinedConstraint<C1, C2>, C3>, C4>,
    C5>
make_combined_constraints(C1 c1, C2 c2, C3 c3, C4 c4, C5 c5) {

  return {{{{c1,c2},c3},c4},c5};
}


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
    // Diagonal edges: include as of HW2 #2
   graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    //#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  
  // Set initial conditions for your nodes, if necessary.
  unsigned int number_of_nodes = graph.num_nodes();
  c = 1.0 / number_of_nodes;
  FixedConstraint<Node> f1, f2;
  for (auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
    
    (*iter).value().mass = 1.0 / number_of_nodes;
    if ((*iter).position() == Point(0, 0, 0)) {
      f1.n = (*iter);
      f1.p = Point(0, 0, 0);
    }
    if ((*iter).position() == Point(1, 0, 0)) {
      
      f2.n = (*iter);
      f2.p = Point(1, 0, 0);
    }

    for (auto edgeiter = (*iter).edge_begin(); edgeiter != (*iter).edge_end();
         ++edgeiter) {
      (*edgeiter).value().L = (*edgeiter).length();
      
    }
  }
   //used in problem 1:
   L = graph.edge(0).length();
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
    double dt = 0.001;
    double t_start = 0.00;
    double t_end = 5.0; 

    for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
      // std::cout << "t = " << t << std::endl;
      //combine constraints
      auto mcc = make_combined_constraints(f1,f2, PlaneConstraint(),
                                           SphereConstraint2(),SphereConstraint()
                                          );
      //mcc(graph,t);
      symp_euler_step(graph, t, dt,
                      make_combined_force(GravityForce(), 
                                         MassSpringForce(),DampingForce()));

      //call constraints
      mcc(graph, t);

      viewer.clear();
      node_map.clear();

      // Update viewer with nodes’ new positions and new edges

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
