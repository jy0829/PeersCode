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

/** Custom structure of data t store with Edges */
struct EdgeData {
  double L;
  double K;
  EdgeData() : L(0), K(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = typename GraphType::size_type;

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
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (!(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))) {
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }
    //n.position() += n.value().vel * dt;
  }
 
  //run the constraits functor
  constraints(g, t);

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
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }    

    // gravitational 
    Point grav_force = n.value().mass * Point(0, 0, -grav); 

    // spring
    Point spring_force = Point(0, 0, 0);
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      spring_force += -(*it).value().K  * (n.position() - (*it).node2().position()) * 
                       ((*it).length() - (*it).value().L) / (*it).length();
    } 
    return grav_force + spring_force;
  }
};

/** Generalized force object that allows the forces to be combined.
 *  This can take up to 3 different forces.
 */
template <typename F1, typename F2, typename F3>
struct Force {
  std::tuple<F1,F2,F3> force3;
  std::tuple<F1, F2> force2;
  int num_forces;
  
  /** Make a tuple of two given forces. */
  void set_forces(F1 f1, F2 f2){
    force2 = std::make_tuple(f1, f2);
    num_forces = 2;
  }

  /** Make a tuple of three given forces. */
  void set_forces(F1 f1, F2 f2, F3 f3) {
    force3 = std::make_tuple(f1, f2, f3);
    num_forces = 3;
  }
   
  /** Return the sum of forces for the given node.  */ 
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (num_forces == 2) {
      return (std::get<0>(force2)(n) + std::get<1>(force2)(n));  
    } else {
      return (std::get<0>(force3)(n) + std::get<1>(force3)(n) + std::get<2>(force3)(n));
    }
  }
};


/** Functor that calculates the gravitaitonal Force. */
struct GravityForce {
  /** Return the gravitational force applied to the node @a n. */
  template <typename NODE>
  Point operator()(NODE n) {
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** Functor that calculates the spring force */
struct MassSpringForce {
  /** Return the spring force applied to the node @a n. */
  template <typename NODE>
  Point operator()(NODE n) {
    Point spring_force = Point(0, 0, 0);
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      spring_force += -(*it).value().K  * (n.position() - (*it).node2().position()) * 
                       ((*it).length() - (*it).value().L) / (*it).length();
    }
    return spring_force;
  }
};

/** Functor that calculates the damping force */
struct DampingForce {
  
  double c;

  DampingForce() : c(1) {
  }

  DampingForce(double c_val) : c(c_val) {
  }

  /**  Return the damping force applied to the node @a n. */
  template <typename NODE>
  Point operator()(NODE n) {
    return -c * n.value().vel;
  }
};


/** Combine the force value of two given forces. 
 * @param[in]     f1     First force functor 
 * @param[in]     f2     Second force functor
 * @return the force object (functor) that holds the combined forces.
 *
 * @tparam F1  Type of the first force functor
 * @tparam F2  Type of the second force functor
 */
template<typename F1, typename F2>
Force<F1, F2, F2> make_combined_force(F1 f1, F2 f2) {
  Force<F1, F2, F2> f;
  f.set_forces(f1, f2);
  return f;
}


/** Combine the force value of three given forces 
 * @param[in]     f1     First force functor 
 * @param[in]     f2     Second force functor
 * @param[in]     f3     Thrid force functor
 * @return the force object (functor) that holds the combined forces.
 *
 * @tparam F1  Type of the first force functor
 * @tparam F2  Type of the second force functor
 * @tparam F3  Type of the thrid force functor
 */
template<typename F1, typename F2, typename F3>
Force<F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  Force<F1, F2, F3> f;
  f.set_forces(f1, f2, f3);
  return f;
}


/** Check all constraints in for the nodes in the given graph. */
template <typename C>
struct Constraints { 
  std::vector<C*> v;

  /** Add a pointer to Constraint object to the contraint vector. */
  void add(C* c) {
    v.push_back(c);   
  } 
 
  /** Apply all constraints in the constraint vector to each node of the 
   *  graph @a g at time @a t. 
   */
  template <typename G>
  void operator()(G& g, double t) {
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      for (size_type i = 0; i < v.size(); i++) {
        Node n = *ni;
        v[i]->apply(g, n);
      }   
    } 
  }
};

/** Abstract class with member function condition and fix.
    Member function apply runs fix if the node meets the condition.    
 */
struct Constraint {

  /** Check if the condition applies to the node @a n. */
  virtual bool condition(GraphType& g, Node& n) = 0; 
  /** Fixes the node @a n. */
  virtual void fix(GraphType& g, Node& n) = 0;
  
  /** Fix the relevant values of node @a n if the condition is met. */
  void apply(GraphType& g, Node& n) {
    if (condition(g, n)) {
      fix(g, n); 
    }
  }
};

/** Constrant of a constant point near (0, 0, 0). */
struct CornerConstraint1 : public Constraint {
  bool condition(GraphType& g, Node& n) {
    return pow(n.position().x - 0, 2) + pow(n.position().y - 0, 2) + pow(n.position().z - 0, 2) < pow(.015, 2);
  }
 
  void fix(GraphType& g, Node& n) {
      n.position() = Point(0, 0, 0);
      n.value().vel = Point(0, 0, 0);
  }
};

/** Constrant of a constant point near (1, 0, 0). */
struct CornerConstraint2 : public Constraint {
  bool condition(GraphType& g, Node& n) {
    return pow(n.position().x - 1, 2) + pow(n.position().y - 0, 2) + pow(n.position().z - 0, 2) < pow(.015, 2);
  }
 
  void fix(GraphType& g, Node& n) {
      n.position() = Point(1, 0, 0);
      n.value().vel = Point(0, 0, 0);
  }
};

/** Constraint of a plane at z = -0.75 . */
struct PlaneConstraint : public Constraint {
  bool condition(GraphType& g, Node& n) {
    return (n.position().z < -0.75);
  }

  void fix(GraphType& g, Node& n) {
    n.position() = Point(n.position().x, n.position().y, -0.75);
    n.value().vel = Point(n.value().vel.x, n.value().vel.y, 0);
  }
};

/** Constraint of a sphere of radius 0.15 centered at (0.5, 0.5, -0.5). */
struct SphereConstraint : public Constraint {
  bool condition(GraphType& g, Node& n) {
    return pow(n.position().x - 0.5, 2) + pow(n.position().y - 0.5, 2) + pow(n.position().z + 0.5, 2) 
           < pow(.15, 2);
  }

  void fix(GraphType& g, Node& n) {
    n.position() = Point(0.5, 0.5, -0.5) + 0.15 * (n.position() - Point(0.5, 0.5, -0.5)) /
                   sqrt(pow(n.position().x - 0.5, 2) + 
                        pow(n.position().y - 0.5, 2) + 
                        pow(n.position().z + 0.5, 2)); 

    Point normal = (n.position() - Point(0.5, 0.5, -0.5)) / 
    sqrt(pow(n.position().x - 0.5, 2) + pow(n.position().y - 0.5, 2) + pow(n.position().z + 0.5, 2));

    n.value().vel -= -(n.value().vel.x * normal.x + n.value().vel.y * normal.y + n.value().vel.z * normal.z)
                     * normal;  
  }
};

/** Constraint of a sphere of radius 0.15 centered at (0.5, 0.5, -0.5)
    that remove the nodes that touches the sphere. 
 */
struct RemoveSphereConstraint : public Constraint {
  bool condition(GraphType& g, Node& n) {
    return pow(n.position().x - 0.5, 2) + pow(n.position().y - 0.5, 2) + pow(n.position().z + 0.5, 2) < pow(.15, 2);
  }

  void fix(GraphType& g, Node& n) {
    g.remove_node(n); 
  }
};



int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // /onstruct an empty graph
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

  // Set initial conditions for your nodes, if necessary.
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) { 
    (*ni).value().mass = 1.0 / graph.num_nodes();
    for (auto ii = (*ni).edge_begin(); ii != (*ni).edge_end(); ++ii) {
      (*ii).value().L = (*ii).length();
      (*ii).value().K = 100;
    }
  } 

  // Initialize constraint objects (pointers) and constraints object.
  Constraints<Constraint> cons;
  Constraint *c1, *c2, *c3, *c4, *c5;
  CornerConstraint1 corner1;
  CornerConstraint2 corner2;
  PlaneConstraint plane; 
  SphereConstraint sphere;
  RemoveSphereConstraint remove;

  c1 = &corner1;
  c2 = &corner2;
  c3 = &plane;
  c4 = &sphere;
  c5 = &remove;
  //cons.add(c1); cons.add(c2); 
  cons.add(c3); cons.add(c4); cons.add(c5);
 
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
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 2.0;

      double c = 1.0 / graph.num_nodes();
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        //P1 symp_euler_step(graph, t, dt, Problem1Force());
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce()));
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c)), cons);
   
        //remove constraint
        viewer.clear();
        node_map.clear();
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
 
        // Update viewer with nodes' new positions
        //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 30)
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
