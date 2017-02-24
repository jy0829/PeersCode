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
 * @tparam G::node_value_type has the publicly visible fields 'vel' and 'mass', 
 *           where mass is some primitive numeric type and vel has type double.
 *
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


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings 
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0);
    Point force;
    force += n.value().mass * Point(0,0,-grav);
    for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
      Point other_node_pos = (*e).node2().position();
      double edge_value = (*e).value().len;
      double dist = norm(n.position() - other_node_pos);
      force += -K_ * (n.position() - other_node_pos) / dist
                   * (dist - edge_value);
    } 
    return force;
  }

  private:
    double K_ = 100;
};

/** Force function object corresponding to gravitational force. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings 
    return n.value().mass * Point(0,0,-grav);
  }
};

/** Force function object corresponding to the spring force on a particular
 *  node due to all the adjacent nodes (and corresponding edges).
 */
struct MassSpringForce {
  MassSpringForce(double K = 100) : K_(K) {
  }
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings 
    Point force;
    for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
      Point other_node_pos = (*e).node2().position();
      double edge_value = (*e).value().len;
      double dist = norm(n.position() - other_node_pos);
      force += -K_ * (n.position() - other_node_pos) / dist
                   * (dist - edge_value);
    }
    return force;
  }
   
  private:
    double K_;
};

/** Force function object corresponding to velocity damping. */
struct DampingForce {
  DampingForce(double c) : c_(c) {
  }

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -c_ * n.value().vel;
  }
  private:
    double c_;
};

/** Force function object whose () operator returns the sum of
 *  the forces output by the two input forces passed in to the constructor.
 */
template <typename F1, typename F2>
struct CombinedForce {
  CombinedForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {
  }

  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n, t);
  }
  
  private:
    F1 f1_;
    F2 f2_;
};

/** Returns a force functor that, given two force functors in the
 *  constructor, calculates the force to be the sum of the forces returned
 *  by the two input functors.
 *
 *  @param f1    The first force functor
 *  @param f2    The second force functor
 *  @type F1     struct or class with a () operation that takes in a Node n and
 *               a double t and returns a point.
 *  @type F2     struct or class with a () operation that takes in a Node n and
 *               a double t and returns a point.
 *  @return      struct that has a () operation that takes in a Node n and
 *               a double t and returns the point f1(n, t) + f2(n, t). 
*/
template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return CombinedForce<F1, F2>(f1, f2);
}

/** Returns a force functor that, given three force functors in the
 *  constructor, calculates the force to be the sum of the forces returned
 *  by the three input functors.
 *
 *  @param f1    The first force functor
 *  @param f2    The second force functor
 *  @param f3    The third force functor
 *  @type F1     struct or class with a () operation that takes in a Node n and
 *               a double t and returns a point.
 *  @type F2     struct or class with a () operation that takes in a Node n and
 *               a double t and returns a point.
 *  @type F3     struct or class with a () operation that takes in a Node n and
 *               a double t and returns a point.
 *  @return      struct that has a () operation that takes in a Node n and a
 *               double t and returns the point f1(n, t) + f2(n, t) + f3(n, t).
*/
template <typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return CombinedForce<CombinedForce<F1, F2>, F3>(CombinedForce<F1, F2>(f1, f2), f3);
}

/** Constraint that fixes all nodes at (0,0,0) or (1,0,0). */
struct FixedConstraint {
  void operator()(GraphType& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (violates_constraint(n, t))
        reset(n);
    } 
  }

  template <typename NODE>
  bool violates_constraint(const NODE n, double t) {
    (void) t;
    return n.position() == Point(0,0,0) || n.position() == Point(1,0,0);
  }

  template <typename NODE>
  void reset(NODE& n) {
    n.value().vel = Point(0,0,0);
  }
};

/** Constraint that blocks nodes from moving below the plane z = -0.75. */
struct PlaneConstraint {
  void operator()(GraphType& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (violates_constraint(n, t))
        reset(n);
    }
  }

  template <typename NODE>
  bool violates_constraint(const NODE n, double t) {
    (void) t;
    return n.position().z < height_;
  }

  template <typename NODE>
  void reset(NODE& n) {
    n.position().z = height_;
    n.value().vel.z = 0;
  }

  private:
    double height_ = -0.75;
};

/** Constraint that blocks nodes from moving inside the sphere of radius
 *  0.15 centered at (0.5, 0.5, -0.5).
 */
struct SphereConstraint {
  void operator()(GraphType& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (violates_constraint(n, t))
        reset(n);
    }
  }

  template <typename NODE>
  bool violates_constraint(const NODE n, double t) {
    (void) t;
    return norm(n.position() - center_) < r_;
  }

  template <typename NODE>
  void reset(NODE& n) {
    Point pos_vec = n.position() - center_;
    pos_vec = pos_vec / norm(pos_vec) * r_;
    n.position() = pos_vec + center_;

    Point ri = pos_vec / r_;
    n.value().vel -= dot(n.value().vel, ri) * ri;
  }

  private:
    Point center_ = Point(0.5, 0.5, -0.5);
    double r_ = 0.15;
};

/** Constraint that removes all nodes that enter the sphere of radius
 *  0.15 centered at (0.5, 0.5, -0.5).
 */
struct TearingSphereConstraint {
  void operator()(GraphType& g, double t) {
    auto it = g.node_begin();
    while (it != g.node_end()) {
      if (violates_constraint(*it, t))
        it = g.remove_node(it);
      else
        ++it;
    }
  }

  template <typename NODE>
  bool violates_constraint(const NODE n, double t) {
    (void) t;
    return norm(n.position() - center_) < r_;
  }

  private:
    Point center_ = Point(0.5, 0.5, -0.5);
    double r_ = 0.15;
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

  // Set initial conditions for nodes
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    (*it).value().mass = 1.0 / graph.num_nodes();
  }

  // Set initial conditions for edges
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    Edge e = *it;
    e.value().len = norm(e.node1().position() - e.node2().position());
  }

  FixedConstraint fixed_constraint;
  PlaneConstraint plane_constraint;
  SphereConstraint sphere_constraint;
  TearingSphereConstraint tearing_sphere_constraint;

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
        symp_euler_step(graph, t, dt,
                          make_combined_force(GravityForce(), MassSpringForce(),
                            DampingForce(1.0 / graph.num_nodes())));
        
        fixed_constraint(graph, t);
        //plane_constraint(graph, t);
        //sphere_constraint(graph, t);
        tearing_sphere_constraint(graph, t);
        
        // Problem 4 updates
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
