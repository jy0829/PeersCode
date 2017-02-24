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
  NodeData() : vel(Point(0,0,0)), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;
  double L;
  EdgeData() : K(100), L(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

struct ConstConstraint {
  // void operator()(Node& n, double t) {
  //   if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))  {
  //      n.value().vel = Point(0,0,0);
  //   } 
  //   (void) t;
  // }

  void operator()(GraphType& g, double t) {
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
      Node n = *i;
      if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))  {
        n.value().vel = Point(0,0,0);
      }
    }
  (void) t;
  }
};

struct PlaneConstraint {
  //   void operator()(Node& n, double t) {
  //     if (dot(n.position(), Point(0,0,1)) < -0.75) {
  //       n.position().elem[2] = -0.75;
  //       n.value().vel.elem[2] = 0;
  //     }
  //   (void) t;
  // }
  void operator()(GraphType& g, double t) {
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
      Node n = *i;
      if (dot(n.position(), Point(0,0,1)) < -0.75) {
        n.position().elem[2] = -0.75;
        n.value().vel.elem[2] = 0;
      }
    }
  (void) t;
  }
};

struct SphereConstraint {
  Point c = Point(0.5,0.5,-0.5);
  double r = 0.15;
  //template<typename NODE>
  void operator()(GraphType& g, double t) {
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
      Node n = *i;
      if (norm(n.position()-c) < r) {
        Point ri = (n.position()-c)/norm(n.position()-c);
        n.position() = c + ri*r;
        n.value().vel += -dot(n.value().vel, ri)*ri;
      }
    }
  (void) t;
  }
};

struct SphereRemove {
  Point c = Point(0.5,0.5,-0.5);
  double r = 0.15;
  //template<typename NODE>
  void operator()(GraphType& g, double t) {
    auto iter = g.node_begin();
    while(iter != g.node_end()) {
      if (norm((*iter).position()-c) < r) {
        iter = g.remove_node(iter);
      } else { ++iter;}

    }
    (void) t;
  }
};

template<typename C1,typename C2>
struct CombineConstraints {
  C1 c1;
  C2 c2;
  //template<typename NODE>
  void operator()(GraphType& g, double t) {
    (void) t;
    c1(g,0);
    c2(g,0);
    //return (c1(n,0)+c2(n,0));
  }
  CombineConstraints(C1 c1, C2 c2) : c1(c1), c2(c2) {}
};

template<typename C1,typename C2>
CombineConstraints<C1,C2> make_combined_constraint(C1 c1, C2 c2) {
  return CombineConstraints<C1,C2>(c1, c2);
}

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
    if (!(n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) ){
      n.position() += n.value().vel*dt;
    }
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    //n.position() = n.position() + n.value().vel * dt;
    //auto con = make_combined_constraint(make_combined_constraint(ConstConstraint(), PlaneConstraint()),SphereConstraint());
    //con(n,0);
  }
  // auto con = make_combined_constraint(ConstConstraint(), PlaneConstraint());
  auto con = make_combined_constraint(make_combined_constraint(ConstConstraint(), PlaneConstraint()),SphereConstraint());
  con(g,0);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    // if (!(n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) )
    //   n.value().vel += force(n, t) * (dt / n.value().mass);
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
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    Point spring_force = Point(0,0,0);
    Point gravity_force = Point(0,0,-grav)*n.value().mass;

    for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
      Edge edge = *e; 
      spring_force += (-edge.value().K)*(edge.node1().position()-edge.node2().position())/edge.length()*(edge.length()-edge.value().L);
    }
    (void) t;
    Point final_force = spring_force+gravity_force;
    return final_force;
  }
};

/** Force function object for gravity */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return (Point(0,0,-grav)*n.value().mass);
  }

};

/** Force function object the mass spring force */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point spring_force = Point(0,0,0);

    for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
      Edge edge = *e; 
      spring_force += (-edge.value().K)*(edge.node1().position()-edge.node2().position())/edge.length()*(edge.length()-edge.value().L);
    }
    (void) t;
    return spring_force;
  }

};

/** Force function object for the damping force */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point damp_force = (-c)*(n).value().vel;
    (void) t;
    return damp_force;
  }
  static double c;
};
double DampingForce::c = 0;

template<typename F1,typename F2>
struct Combine {
  F1 f1;
  F2 f2;
  template<typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return (f1(n,0)+f2(n,0));
  }
  Combine(F1 f1, F2 f2) : f1(f1), f2(f2) {}

};

template<typename F1,typename F2>
Combine<F1,F2> make_combined_force(F1 f1, F2 f2) {
  return Combine<F1,F2>(f1, f2);
}

template<typename F1,typename F2,typename F3>
Combine<Combine<F1,F2>,F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return Combine<Combine<F1,F2>,F3>(make_combined_force(f1, f2), f3);
}

//Point make_combined_force(Point f1, Point f2) {
//  return f1+f2;  
//}

// Point make_combined_force(F1 f1, F2 f2, F3 f3) {
//   Point operator()(NODE n, double t) {
//     (void) t;
//     return (f1(n,0)+f2(n,0)+f3(n,0)); 
//   } 
// }


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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  DampingForce::c = 1.0/graph.size();

  for (auto n = graph.node_begin(); n != graph.node_end(); ++n) {
    (*n).value().mass = 1.0/graph.size();
    (*n).value().vel = Point(0,0,0);
  }

  for (auto e = graph.edge_begin(); e != graph.edge_end(); ++e) {
    (*e).value().K = 100;
    (*e).value().L = (*e).length();
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
        //symp_euler_step(graph, t, dt, Problem1Force());
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce(),DampingForce() ));

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