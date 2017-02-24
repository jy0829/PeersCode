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
#include <cassert>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// Damping constant
double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Custom structure of data to store with Edges
struct EdgeData {
  double K;
  double length;
  EdgeData() : K(0), length(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
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
 * @tparam G::node_value_type supports the NodeData struct
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
    if (!(n.position() == Point(0,0,0)
          || n.position() == Point(1,0,0))) {
      n.position() += n.value().vel * dt;
    }
  }

  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (!(n.position() == Point(0,0,0)
          || n.position() == Point(1,0,0))) {
      n.value().vel += force(n,t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}

// Gravity object
struct GravityForce {
  // Operator to assign gravity to a node
  template <typename NODE>
  Point operator()(NODE n, double) const {
    return Point(0,0,-grav*n.value().mass);
  }
};


// Spring force object
struct MassSpringForce {
  /** Return the force applying to @a n at time @a t.  */
  template <typename NODE>
  Point operator()(NODE n, double) const {
    // Calculate spring force for the node
    Point spring_force;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      auto node1 = e.node1();
      auto node2 = e.node2();
      assert(node1 == n);
      Point diff = node1.position() - node2.position();
      double distance = norm(diff);
      spring_force += (-1)*(e.value().K)*(diff)*(1-(e.value().length)/distance);
    }
    return spring_force;
  }
};


// Damping force object
struct DampingForce {
  // Operator to assign damping force to node based on damping constance c
  template <typename NODE>
  Point operator()(NODE n, double) const {
    return (n.value().vel)*c;
  }
};


// Creating type Force2 to combine 2 forces
template <typename F1, typename F2>
struct Force2 {
  F1 f1_;
  F2 f2_;
  Force2<F1,F2>(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n,t) + f2_(n,t);
  }
};


// Function to combine 2 forces
template <typename F1, typename F2>
Force2<F1,F2> make_combined_force(F1 f1, F2 f2) {
  return Force2<F1,F2>(f1, f2);
}

// Function to combine 3 forces
template <typename F1, typename F2, typename F3>
Force2<F1,Force2<F2,F3>> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return Force2<F1,Force2<F2,F3>>(f1, Force2<F2,F3>(f2,f3));
}



// Constraint for plane z
struct planeConstraint {
  double z_;
  planeConstraint(double z) : z_(z) {}
  void operator()(GraphType& g, double) {
    for (auto it = g.node_begin(); it!=g.node_end(); ++it) {
      auto n = *it;
      if (inner_prod(n.position(), Point(0,0,1)) < z_) {
        n.position().z = z_;
        n.value().vel.z = 0;
      }
    }
  }
};


// Constraint for the sphere
struct sphereConstraint {
  Point c_;
  double r_;
  sphereConstraint(Point c, double r) : c_(c), r_(r) {}
  void operator()(GraphType& g, double) {
    for (auto it = g.node_begin(); it!=g.node_end(); ++it) {
      auto n = *it;
      Point initpos = n.position();
      Point velocity = n.value().vel;
      double distance = norm(n.position() - c_);
      Point R = (initpos-c_)/distance;
      if (distance < r_) {
        n.position() = (initpos-c_)/distance*r_ + c_;
        n.value().vel = velocity - inner_prod(velocity,R)*R;
      }
    }
  }
};

// Remove node constraint
struct sphereRemove {
  Point c_;
  double r_;
  sphereRemove(Point c, double r) : c_(c), r_(r) {}
  void operator()(GraphType& g, double) {
    for (auto it = g.node_begin(); it!=g.node_end(); ++it) {
      auto n = *it;
      double distance = norm(n.position() - c_);
      if (distance < r_) {
        g.remove_node(n);
      }
    }
  }
};

// Combine 2 constraints struct
template <typename C1, typename C2>
struct constraint2 {
  C1 c1_;
  C2 c2_;
  constraint2<C1,C2>(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}
  void operator()(GraphType& g, double t) {
    c1_(g,t);
    c2_(g,t);
  }
};

// Function to combine 2 constraints
template <typename C1, typename C2>
constraint2<C1,C2> combine_constraint(C1 c1, C2 c2) {
  return constraint2<C1,C2>(c1,c2);
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

  // Set initial velocity and mass of the node
  c = 1/double(graph.num_nodes());
  for (auto it = graph.node_begin(); it!=graph.node_end(); ++it) {
    auto n = *it;
    n.value().vel = Point(0,0,0);
    n.value().mass = c;
  }

  // Set the spring constant and length of the edge
  for (auto eit = graph.edge_begin(); eit!=graph.edge_end(); ++eit) {
    auto e = *eit;
    e.value().K = 100;
    e.value().length = norm(e.node1().position()-e.node2().position());
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
        //std::cout << "t = 90" << t << std::endl;
        symp_euler_step(graph, t, dt, make_combined_force(MassSpringForce(),
                        DampingForce(),GravityForce()),
                        combine_constraint(sphereRemove(Point(0.5,0.5,-0.5),0.15),
                        planeConstraint(-0.75)));

        // Clear the viewer's nodes and edges
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
