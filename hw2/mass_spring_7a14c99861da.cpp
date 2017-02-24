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
#include <iostream>

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
  NodeData(Point v, double m, Point org) : vel(v), mass(m), origin(org) {}
  Point origin;
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;
  double L;
  EdgeData() : K(100), L(1) {}
  EdgeData(double Kval, double Lval) : K(Kval), L(Lval) {}
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
 * @tparam G::node_value_type supports the struct NodeData: has mass and
 *           velocity as attributes.
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

  // Compute the t+dt position, skipping (0,0,0) and (1,0,0)
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
  }
 
  //reset any nodes that violate the constraint
  constraint(g, t);
  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Constraint that keeps point (0,0,0) and (1,0,0) fixed. */
struct Pins{
  void operator()(GraphType& Graph, double t) {
    (void) t;
    for (auto it = Graph.node_begin(); it != Graph.node_end(); ++it) {
      Node n = *it;
      if (n.value().origin == Point(0,0,0) || n.value().origin == Point(1,0,0)) {
        n.position() = n.value().origin;
      }
    }
  }
};

/** Constraint in the form of a plane with z = -0.75. */
struct Plane{
  void operator()(GraphType& Graph, double t) {
    (void) t;
    for (auto it = Graph.node_begin(); it != Graph.node_end(); ++it) {
      Node n = *it;
      if (dot(n.position(), Point(0,0,1)) < -0.75) {
        n.position().z = -0.75;
        n.value().vel.z = 0;
      }
    }
  }
};

/** Constraint in the form of a sphere with center at (0.5,0.5,-0.5)
    and radius = 0.15. */
struct Sphere{
  void operator()(GraphType& Graph, double t) {
    (void) t;
    Point c(0.5,0.5,-0.5);
    for (auto it = Graph.node_begin(); it != Graph.node_end(); ++it) {
      Node n = *it;
      if (norm(n.position() - c) < 0.15) {
        Point R = (n.position() - c) / norm(n.position() - c);
        n.position() = (n.position() - c)/(double)norm(n.position() - c) * 0.15 
          + c;
        n.value().vel = n.value().vel - dot(n.value().vel, R) * R;
      }
    }
  }
};

/** Constraint in the form of a sphere centered at (0.5,0.5,-0.5) with
    radius = 0.15 that destroys nodes when it comes into contact with them. */
struct SphereDestroy{
  void operator()(GraphType& Graph, double t) {
    (void) t;
    Point c(0.5,0.5,-0.5);
    for (auto it = Graph.node_begin(); it != Graph.node_end(); ++it) {
      Node n = *it;
      if (norm(n.position() - c) < 0.15) {
        Graph.remove_node(it);
      }
    }
  }
};

/** Constraint struct with functor that doesn't change the graph. */
struct NullConstraint{
  void operator()(GraphType& Graph, double t) {
    (void) t;
    (void) Graph;
  }
};

/** struct with a functor that combines two constraints into one.  */
template <typename C1, typename C2, typename C3>
struct CombConstraint {
  C1 c1;
  C2 c2;
  C3 c3;
  /** @pre fixing one constraint does not cause the other constraint
   *  to be violated.
   *  @post No node in the graph is violating either of the two constraints.
  */
  CombConstraint(C1 c1val) : c1(c1val) {
    c2 = NullConstraint();
    c3 = NullConstraint();
  }
  CombConstraint(C1 c1val, C2 c2val) : c1(c1val), c2(c2val) {
    c3 = NullConstraint();
  }
  CombConstraint(C1 c1val, C2 c2val, C3 c3val) : c1(c1val), c2(c2val),
    c3(c3val) {}
  void operator()(GraphType& Graph, double t) {
    c1(Graph, t);
    c2(Graph, t); 
    c3(Graph, t);
  }
};

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }
    Point spring_force(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Edge e = *it;
      spring_force += -100 * (n.position() - e.node2().position()) * 
        (norm(n.position() - e.node2().position()) - e.value().L) / 
        norm(n.position() - e.node2().position());
    }
    (void) t;
    return spring_force + n.value().mass * Point(0,0,-grav);
  }
};

/** struct with a functor that returns the force due to gravity on
    the given node. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().mass * Point(0,0,-grav);
  }
};

/** struct with a functor that returns the force due to the springs on
    the given node. */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point spring_force(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Edge e = *it;
      spring_force += -e.value().K * (n.position() - e.node2().position()) * 
        (norm(n.position() - e.node2().position()) - e.value().L) / 
        norm(n.position() - e.node2().position());
    }
    return spring_force;
  }
};

/** struct with a functor that returns the force due to damping
    the given node. */
struct DampingForce {
  double c;
  DampingForce(double cval) : c(cval) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -c * n.value().vel;  
  }
};

/** struct with a functor that returns (0,0,0) as the force on the given
    node. */
struct NullForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    (void) n;
    return Point(0);
  }
};

/** struct that can construct a null force object, or a combined force
    object with one, two, or three forces.  */
template <typename F1, typename F2, typename F3>
struct CombForce {
  F1 f1;
  F2 f2;
  F3 f3;
  CombForce() : f1(NullForce()), f2(NullForce()), f3(NullForce()) {}
  CombForce(F1 f1val) : f1(f1val), f2(NullForce()), f3(NullForce()) {}
  CombForce(F1 f1val, F2 f2val) : f1(f1val), f2(f2val), f3(NullForce()) {}
  CombForce(F1 f1val, F2 f2val, F3 f3val) : f1(f1val), f2(f2val), f3(f3val) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return(f1(n, t) + f2(n, t) + f3(n, t)); 
  } 
};

/** Function to return a functor that calculates the total force on a given node
 * of the given forces.
 * @param f1    A force functor
 * @param f2    A force functor
 * @type F1     structs with a () operation that takes in a Node n and 
 *                  double t and returns a point.
 * @type F2     structs with a () operation that takes in a Node n and 
 *                  double t and returns a point.
 * @return      functor that has a () operation that takes in a Node n and 
 *                  double t and returns f1(n, t) + f2(n, t). 
*/
template <typename F1, typename F2>
CombForce<F1, F2, NullForce> make_combined_force(F1 f1, F2 f2) {
  return CombForce<F1, F2, NullForce>(f1, f2);
}

/** Function to return a functor that calculates the total force on a given node
 * of the given forces.
 * @param f1    A force functor
 * @param f2    A force functor
 * @param f3    A force functor
 * @type F1     structs with a () operation that takes in a Node n and 
 *                  double t and returns a point.
 * @type F2     structs with a () operation that takes in a Node n and 
 *                  double t and returns a point.
 * @type F3     structs with a () operation that takes in a Node n and 
 *                  double t and returns a point.
 * @return      functor that has a () operation that takes in a Node n and 
 *                  double t and returns f1(n, t) + f2(n, t) + f3(n, t). 
*/
template <typename F1, typename F2, typename F3>
CombForce<F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return CombForce<F1, F2, F3>(f1, f2, f3);
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

  //set initial conditions
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    Node n = *ni;
    n.value() = NodeData(Point(0), (double)1/(double)graph.size(), n.position());
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
      Edge e = *ii;
      e.value() = EdgeData(100, norm(e.node1().position() - 
        e.node2().position()));
  
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
        //symp_euler_step(graph, t, dt, Problem1Force(), NullConstraint());
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), 
          DampingForce((double)1/(double)graph.size()), MassSpringForce()), 
          //comment out one of the following lines based on which constraints
          //you want
          CombConstraint<Pins, SphereDestroy, NullConstraint>(Pins(), SphereDestroy()));
          //CombConstraint<Pins, Plane, Sphere>(Pins(), Plane(), Sphere()));

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
