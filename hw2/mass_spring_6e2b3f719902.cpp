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
#include <utility>

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

/*
 * Custom structure for Edge.
 */
struct EdgeData {
  double K;
  double L;
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
 * @param[in]     constraint Function object for the constaints on the graph.
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

  // HW2Q4: constraint processing.
  constraint(g, t);

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    #if 0
    // Skip updates for two fixed points. Deprecated.
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      continue;
    }
    #endif

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
    // HW2 #1: YOUR CODE HERE
    (void)t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }
    Point force(0);
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      Edge edge = *iter;
      Node n1 = edge.node1();
      Node n2 = edge.node2();
      if (n2 == n) {
        std::swap(n1, n2);
      }
      force += -edge.value().K * (n1.position() - n2.position()) /
      edge.length() * (edge.length() - edge.value().L);
    }
    force += n.value().mass * Point(0, 0, -grav);
    return force;
  }
};

/*
 * HW2 #3: GravityForce.
 */
struct GravityForce {
  template<typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    return n.value().mass * Point(0, 0, -grav);
  }
};

/*
 * HW2 #3: MassSpringForce.
 */
struct MassSpringForce {
  template<typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    Point force(0);
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      Edge edge = *iter;
      Node n1 = edge.node1();
      Node n2 = edge.node2();
      if (n2 == n) {
        std::swap(n1, n2);
      }
      force += -edge.value().K * (n1.position() - n2.position()) /
      edge.length() * (edge.length() - edge.value().L);
    }
    return force;
  }
};

/*
 * HW2 #3: DampingForce.
 */
struct DampingForce {
  template<typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    return -c * n.value().vel;
  }
  double c;
  DampingForce(double c = 0.0) : c(c) {}
};

/*
 * HW2 #3: Combine forces.
 */
template<typename Force1, typename Force2>
struct CombinedForce {
  template<typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n, t) + f2(n, t);
  }
  Force1 f1;
  Force2 f2;
  CombinedForce(Force1 f1, Force2 f2) : f1(f1), f2(f2) {}
};

template<typename Force1, typename Force2>
CombinedForce<Force1, Force2> make_combined_force(
    Force1 f1, Force2 f2) {
  return CombinedForce<Force1, Force2>(f1, f2);
}

template<typename Force1, typename Force2, typename Force3>
CombinedForce<CombinedForce<Force1, Force2>, Force3> make_combined_force(
    Force1 f1, Force2 f2, Force3 f3) {
  return make_combined_force(make_combined_force(f1, f2), f3);
}

//HW2Q4: Constant node constraint to fix certain nodes.
struct ConstantNodeConstraint {
  template<typename G>
  void operator()(G& g, double t) {
    (void)t;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      Node node = *iter;
      if (find(constant_points.begin(), constant_points.end(),
          node.position()) != constant_points.end()) {
        node.value().vel = Point(0);
      }
    }
  }
  std::vector<Point> constant_points;
  ConstantNodeConstraint(std::vector<Point> constant_points = 
      std::vector<Point>({Point(0, 0, 0), Point(1, 0, 0)})) :
      constant_points(constant_points) {}
};

//HW2Q4: z-plane constraint.
struct ZPlaneConstraint {
  template<typename G>
  void operator()(G& g, double t) {
    (void)t;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      Node node = *iter;
      if (node.position().elem[2] < z) {
        node.position().elem[2] = z;
        node.value().vel[2] = 0;
      }
    }
  }
  double z;
  ZPlaneConstraint(double z = -0.75) : z(z) {}
};

//HW2Q4: sphere constraint.
struct SphereConstraint {
  template<typename G>
  void operator()(G& g, double t) {
    (void)t;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
      Node node = *iter;
      double dis = norm(node.position() - c);
      if (dis < r) {
        Point R = (node.position() - c) / dis;
        node.position() = c + r * R;
        node.value().vel -= dot(node.value().vel, R) * R;
      }
    }
  }
  Point c;
  double r;
  SphereConstraint(const Point& c = Point(0.5, 0.5, -0.5),
      double r = 0.15) : c(c), r(r) {}
};

//HW2Q4: sphere constraint with removing nodes.
struct SphereConstraintRemoveNodes {
  template<typename G>
  void operator()(G& g, double t) {
    (void)t;
    for (auto iter = g.node_begin(); iter != g.node_end();) {
      Node node = *iter;
      double dis = norm(node.position() - c);
      if (dis < r) {
        g.remove_node(node);
      } else {
        ++iter;
      }
    }
  }
  Point c;
  double r;
  SphereConstraintRemoveNodes(
      const Point& c = Point(0.5, 0.5, -0.5), double r = 0.15) :
      c(c), r(r) {}
};

struct EmptyConstraint {
  template<typename G>
  void operator()(G& g, double t) {
    (void)g;
    (void)t;
  }
};

/*
 * HW2#4: Combined constraint.
 */
template<typename Constraint1, typename Constraint2>
struct CombinedConstraint {
  template<typename G>
  void operator()(G& g, double t) {
    c1(g, t);
    c2(g, t);
  }
  Constraint1 c1;
  Constraint2 c2;
  CombinedConstraint(Constraint1 c1, Constraint2 c2) : c1(c1), c2(c2) {}
};

template<typename Constraint1, typename Constraint2>
CombinedConstraint<Constraint1, Constraint2> make_combined_constraint(
    Constraint1 c1, Constraint2 c2) {
  return CombinedConstraint<Constraint1, Constraint2>(c1, c2);
}

template<typename Constraint1, typename Constraint2, typename Constraint3>
CombinedConstraint<CombinedConstraint<Constraint1, Constraint2>, Constraint3>
    make_combined_constraint(Constraint1 c1, Constraint2 c2, Constraint3 c3) {
  return make_combined_constraint(make_combined_constraint(c1, c2), c3);
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
  double N = graph.num_nodes();
  for (auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
    Node n = *iter;
    n.value().vel = Point(0);
    n.value().mass = 1.0 / N;
  }

  // HW2 #2 add initial conditions for edges.
  for (auto iter = graph.edge_begin(); iter != graph.edge_end(); ++iter) {
    Edge edge = *iter;
    edge.value().K = 100.0;
    edge.value().L = edge.length();
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
        #if 0
        symp_euler_step(graph, t, dt, Problem1Force(), EmptyConstraint());
        #endif

        #if 0
        symp_euler_step(graph, t, dt,
            make_combined_force(GravityForce(), MassSpringForce()),
            EmptyConstraint());
        #endif

        #if 0
        symp_euler_step(graph, t, dt,
            make_combined_force(GravityForce(), MassSpringForce(),
            DampingForce(1.0 / N)), EmptyConstraint());
        #endif

        #if 0
        symp_euler_step(graph, t, dt,
            make_combined_force(GravityForce(), MassSpringForce()),
            ZPlaneConstraint());
        #endif

        #if 0
        symp_euler_step(graph, t, dt,
            make_combined_force(GravityForce(), MassSpringForce()),
            SphereConstraint());
        #endif

        #if 0
        symp_euler_step(graph, t, dt,
            make_combined_force(GravityForce(), MassSpringForce()),
            make_combined_constraint(ZPlaneConstraint(),
            SphereConstraint()));
        #endif

        #if 0
        symp_euler_step(graph, t, dt,
            make_combined_force(GravityForce(), MassSpringForce()),
            SphereConstraintRemoveNodes());
        #endif

        #if 0
        symp_euler_step(graph, t, dt,
            make_combined_force(GravityForce(), MassSpringForce()),
            make_combined_constraint(ConstantNodeConstraint(),
            ZPlaneConstraint(), SphereConstraintRemoveNodes()));
        #endif

        #if 1
        symp_euler_step(graph, t, dt,
            make_combined_force(GravityForce(), MassSpringForce(),
            DampingForce(1.0 / N)),
            make_combined_constraint(ConstantNodeConstraint(),
            ZPlaneConstraint(), SphereConstraintRemoveNodes()));
        #endif

        // Update viewer with nodes' new positions
        viewer.clear();
        node_map.clear();
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
