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
double c = 0.0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edge */
struct EdgeData {
  double K;   //< Spring constant
  double L;   //< Resting length
  EdgeData() : K(100.0), L(0.0) {}
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

  // Apply the constraints
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/* Structure for defining constant constraints on specific nodes */
struct node_constraint {
  Point pos;
  Point vel;
  Node node;
};

/** Constant node constraint */
struct ConstConstraint {
  ConstConstraint(std::vector<node_constraint> n_in) : n(n_in) {}

  template<typename GRAPH>
  void operator()(GRAPH&, double) {
    for (unsigned int i=0; i < n.size(); ++i) {
      n[i].node.position() = n[i].pos;
      n[i].node.value().vel = n[i].vel;
    }
  }

  std::vector<node_constraint> n;
};        

/** Plane constraint */
struct PlaneConstraint {
  PlaneConstraint() : z(-0.75) {}

  template<typename GRAPH>
  void operator()(GRAPH& g, double) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if ( dot(n.position(),Point(0,0,1)) < z) {
        // Project point onto plane
        n.position().z = z;
        // Remove velocity in z
        n.value().vel.z = 0.0;
        //std::cout << n.position() << std::endl;
      }
    }
  }

  double z;
};

/** Sphere constraint */
struct SphereConstraint {
  SphereConstraint() : c(Point(0.5,0.5,-0.5)), r(0.15) {}

  template<typename GRAPH>
  void operator()(GRAPH& g, double) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point dir = n.position() - c;
      double dist = norm(dir);
      if (dist < r) {
        // Project point onto sphere
        dir *= r/dist;
        n.position() = c + dir;
        // Set velocity normal to sphere equal to zero
        n.value().vel -= dot(n.value().vel,dir)*dir/(r*r);
      }
    }
  } 

  Point c;
  double r;
};

/* Sphere Removal Constraint */
struct SphereRmConstraint {
  SphereRmConstraint() : c(Point(0.5,0.5,-0.5)), r(0.15) {}

  template<typename GRAPH>
  void operator()(GRAPH& g, double) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      double dist = norm(n.position() - c);
      if (dist < r) {
        g.remove_node(n);
      }
    }
  }

  Point c;
  double r;
};

// Combine Constraints
template<typename C1, typename C2, typename C3>
struct Constraint {
  Constraint(const C1& a,const C2& b,const C3& c) : a_(a),b_(b),c_(c) {}

  template<typename GRAPH>
  void operator()(GRAPH& g, double t) {
    a_(g,t);
    b_(g,t);
    c_(g,t);
  }
  C1 a_;
  C2 b_;
  C3 c_;

};

template<typename C1, typename C2, typename C3>
Constraint<C1,C2,C3> make_combined_constraint(const C1& a, const C2& b, const C3& c) {
  return Constraint<C1,C2,C3>(a,b,c);
}

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double) {
    Point f = {0.0,0.0,0.0};
   
    // skip for points (0,0,0) and (1,0,0)
    if ((n.position() != Point(0,0,0)) and (n.position() != Point(1,0,0))) {

      // Accumulate Spring force from adjacent nodes
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        auto e = *it;
        auto n2 = e.node2();
        f +=  -1.0*e.value().K*(n.position()-n2.position())/e.length()*(e.length() -
            e.value().L);
      }
      // Add graviational force
      f += Point(0.0,0.0,-1.0*grav)*n.value().mass;
    }

    return f;
  }
};

/* Gravity Force Function */
struct GravityForce {
  /* Return the gravitational force applying to @a n */
  template <typename NODE>
  Point operator()(NODE n, double) {
    return {0.0,0.0,n.value().mass*-1.0*grav};
  }
};

/* Spring Force function */
struct MassSpringForce {
  /* Return the mass spring force applying to @a n at time @a t */
  template<typename NODE>
  Point operator()(NODE n, double) {
    Point f = {0.0,0.0,0.0};
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      auto n2 = e.node2();
       f +=  -1.0*e.value().K*(n.position()-n2.position())/e.length()*(e.length() -             e.value().L);
    }
    return f;
  }
};

/* Damping force function */
struct DampingForce {
  /* Return the damping force applying to @a n at time @a t */
  template<typename NODE>
  Point operator()(NODE n, double) {
    return n.value().vel*-1.0*c;
  }
};

// Combine forces
template<typename F1, typename F2>
struct Force2 {
  template<typename NODE>
  Point operator()(NODE n, double t) {
    return a(n,t)+b(n,t);
  }
  F1 a;
  F2 b;
};
template<typename F1, typename F2, typename F3>
struct Force3 {
  template<typename NODE>
  Point operator()(NODE n, double t) {
    return a(n,t)+b(n,t)+c(n,t);
  }
  F1 a;
  F2 b;
  F3 c;
};
 
// Create combined force functors
template<typename F1, typename F2>
Force2<F1,F2> make_combined_force(const F1&, const F2&) {
  return Force2<F1,F2>();
}
template<typename F1, typename F2, typename F3>
Force3<F1,F2,F3> make_combined_force(const F1&, const F2&, const F3&) {
  return Force3<F1,F2,F3>();
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
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set node masses and resting edge lengths (also constant constraints)
  std::vector<node_constraint> nc;
  double N = graph.size();
  for (auto it=graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().mass = 1.0/N;
    for (auto it2=n.edge_begin(); it2 != n.edge_end(); ++it2) {
      (*it2).value().L = (*it2).length();
    }
    // constant constraints
    if (n.position() == Point(0,0,0)) {
      nc.push_back({Point(0,0,0),Point(0,0,0),n});
    }
    if (n.position() == Point(1,0,0)) {
      nc.push_back({Point(1,0,0),Point(0,0,0),n});
    } 
  }
 
  // Set the damping constant
  c = 1.0/N;

  // Create combined force
  auto f = make_combined_force(GravityForce(),MassSpringForce(),DampingForce());

  // Create combined constraint
  //auto constraint = make_combined_constraint(ConstConstraint(nc),PlaneConstraint(),SphereConstraint());
  auto constraint = make_combined_constraint(ConstConstraint(nc),PlaneConstraint(),SphereRmConstraint());

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
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, f, constraint);
        
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
