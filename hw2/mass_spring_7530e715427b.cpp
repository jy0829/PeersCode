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

// Damping force
static double damp = 1.0;

// Uids of fixed points
static unsigned int uid1 = 0;
static unsigned int uid2 = 0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double len;
  double springk;
  EdgeData() : len(0), springk(100) {}

};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using NodeIterator = typename GraphType::node_iterator;
using EdgeIterator = typename GraphType::edge_iterator;
using IncidentIterator = typename GraphType::incident_iterator;
typedef Point (*fptr)(Node n, double t);

// CONSTRAINTS

/** Constant fixed points constraint object. */
struct ConstantConstraint {
  /* Apply the constraint that the two corner positions are held fixed
   */
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    Node n1 = g.node(g.uid2idx(uid1));
    Node n2 = g.node(g.uid2idx(uid2));
    n1.position() = Point(0,0,0);
    n1.value().vel = Point(0,0,0);

    n2.position() = Point(1,0,0);
    n2.value().vel = Point(0,0,0);
  }

};

/** Plane constraint object. */
struct PlaneConstraint {
  /* Apply the constraint that points cannot fall below the plane z = -0.75
   */
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    for (auto n : nodesof(g)) {
      if (n.position().z < -0.75) {
        n.position().z = -0.75;
        n.value().vel.z = 0;
      }
    }
  }
};

/** Sphere constraint object. */
struct SphereConstraint {
  /* Apply the constraint that points cannot go inside a sphere of radius 0.15
   * centred at (0.5, 0.5, -0.5).
   */
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    for (auto n : nodesof(g)) {
      Point vec = n.position() - Point(0.5,0.5,-0.5);
      if(norm(vec) < 0.15) {
        Point Ri = (1/norm(vec))*vec;
        n.position() += (0.15-norm(vec))*Ri;
        n.value().vel -= dot(n.value().vel, Ri)*Ri;
      }
    }
  }
};

/** Sphere constraint object that removes nodes. */
struct SphereRemoveConstraint {
  /* Apply the constraint that points cannot go inside a sphere of radius 0.15
   * centred at (0.5, 0.5, -0.5). Nodes are removed if they violate this.
   */
  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    (void) t;
    for (auto n : nodesof(g)) {
      Point vec = n.position() - Point(0.5,0.5,-0.5);
      if(norm(vec) < 0.15) {
        g.remove_node(n);
      }
    }
  }
};

/** Combined constraint object. */
template <typename T, typename T2>
struct CombinedConstraint {
  T c1_;
  T2 c2_;
  /** Return the combined force applying to @a n at time @a t.
  */
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    c1_(g,t);
    c2_(g,t);
  }

  CombinedConstraint(T c1, T2 c2) : c1_(c1), c2_(c2) {}
};

/** Function to combine two constraints into a CombinedConstraint object */
template <typename T, typename T2>
CombinedConstraint<T, T2> make_combined_constraint(T first, T2 second) {
  return CombinedConstraint<T, T2>(first, second);
}

/** Function to combine three constraints into a CombinedConstraint object */
template<typename T, typename T2, typename T3>
CombinedConstraint<CombinedConstraint<T, T2>, T3> make_combined_constraint(T first, T2 second, T3 third) {
  return CombinedConstraint<CombinedConstraint<T, T2>, T3>(CombinedConstraint<T, T2>(first, second), third);
}


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports velocity and mass
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto n : nodesof(g)) {
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  // Create constraint
  auto constraint = make_combined_constraint(ConstantConstraint(), PlaneConstraint(), SphereRemoveConstraint());
  //ConstantConstraint constraint;
  // Apply constraints
  constraint(g, t);
  
  // Compute the t+dt velocity
  for (auto n : nodesof(g)) {
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
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    } else {
      Point grav_force(0,0,-n.value().mass*grav);
      Point spring_force(0,0,0);
      for (auto e : n) {
        Node n2 = e.node2();
        Point p = n.position() - n2.position();
        double dist = norm(p);
        p /= dist;
        p *= (dist - e.value().len)*(-e.value().springk);
        spring_force += p;
      }

      return spring_force+grav_force;
    }
  }
};

/** Gravity force function object. */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point grav_force(0,0,-n.value().mass*grav);
    return grav_force;
  }

};

/** Mass force function object. */
struct MassSpringForce {
  /** Return the mass force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point spring_force(0,0,0);
    for (auto e : n) {
      Node n2 = e.node2();
      Point p = n.position() - n2.position();
      double dist = norm(p);
      p /= dist;
      p *= (dist - e.value().len)*(-e.value().springk);
      spring_force += p;
    }

      return spring_force;
  }

};

/** Damping force function object. */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point damping_force = -damp*n.value().vel;
    return damping_force;
  }
};

/** Combined force function object. */
template <typename T, typename T2>
struct CombinedForce {
  T f1_;
  T2 f2_;
  /** Return the combined force applying to @a n at time @a t.
  */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n,t);
  }

  CombinedForce(T f1, T2 f2) : f1_(f1), f2_(f2) {}
};

/** Function to combine two forces into a CombinedForce object */
template <typename T, typename T2>
CombinedForce<T, T2> make_combined_force(T first, T2 second) {
  return CombinedForce<T, T2>(first, second);
}

/** Function to combine three forces into a CombinedForce object */
template<typename T, typename T2, typename T3>
CombinedForce<CombinedForce<T, T2>, T3> make_combined_force(T first, T2 second, T3 third) {
  return CombinedForce<CombinedForce<T, T2>, T3>(CombinedForce<T, T2>(first, second), third);
}

// Attempt at variadic
/*template<typename T, typename T2>
CombinedForce<T, T2> make_combined_force(T first, T2 second) {
  return CombinedForce<T, T2>(first, second);
}

template<typename T, typename T2, typename... Args>
CombinedForce make_combined_force(T first, T2 second, Args... args) {
  return make_combined_force(CombinedForce<T,T2>(first,second), args...);
}*/

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
  damp = 1/graph.num_nodes();

  // Set initial conditions for your nodes, if necessary.
  for (auto n : nodesof(graph)) {
    NodeData node_data;
    node_data.vel = Point(0,0,0);
    node_data.mass =  1.0/(float)nodes.size();
    n.value() = node_data;
    if (n.position() == Point(0,0,0)) {
      uid1 = n.uid();
    } else if (n.position() == Point(1,0,0)) {
      uid2 = n.uid();
    }
    for (auto e : n) {
      EdgeData edge_data;
      edge_data.len = e.length();
      edge_data.springk = 100;
      e.value() = edge_data;
    }
  }

  // Set initial conditions for edges.
  /*for (auto e : edgesof(graph)) {
    EdgeData edge_data;
    edge_data.len = e.length();
    edge_data.springk = 100;
    e.value() = edge_data;
  }*/
  // Get the size of L, which for Problem 1 we set to the distance of the first edge
  //Edge firste = *graph.edge_begin();
  //double L = firste.length();

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
        // Problem 1 & 2
        //symp_euler_step(graph, t, dt, Problem1Force());
        // Problem 3-5
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()));

        // Clear the viewer's nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions and new edges
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

        // Update viewer with nodes' new positions
        //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
