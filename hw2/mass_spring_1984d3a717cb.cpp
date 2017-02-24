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
  double K;       //< Edge spring constant
  double L;     //< Edge rest length
  EdgeData() : K(100), L(0) {}
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

/** the original euler_step function
 */
#if 0
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
#endif

/** the euler_step function for Problem 3 which skip two points
 */
#if 0
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    // skip two points
    if (n.position() != Point(0,0,0) and n.position() != Point(1,0,0)) {
      n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    // skip two points
    if (n.position() != Point(0,0,0) and n.position() != Point(1,0,0)) {
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}
#endif

/** the symp_euler_step function for Problem 4 which has constraints
 */
#if 1
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
  }
  constraint(g,t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
   
      n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}
#endif


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
    // it seems that t is not needed, since f_spring is not explicitly 
    // dependent on t
    (void)t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0);
    } 

    Point f_grav = Point(0,0,-grav*n.value().mass);
    Point f_spring = Point(0);
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      f_spring += -K_*(e.length()-L_)/e.length() * 
        (e.node1().position()-e.node2().position());
    }
    return f_spring + f_grav;
  }

  Problem1Force (double K, double L) {
    K_ = K;
    L_ = L;
  }
  double K_;
  double L_;
};

/** Force function object for HW2 #2. */
struct Problem2Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #2, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    // it seems that t is not needed, since f_spring is not explicitly 
    // dependent on t
    (void)t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {      
      return Point(0);
    }

    Point f_grav = Point(0,0,-grav*n.value().mass);      
    Point f_spring = Point(0);  
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      f_spring += -e.value().K * (e.length()-e.value().L)/e.length() * 
        (e.node1().position()-e.node2().position());
    }
    return f_spring + f_grav;
  }
};

/** Force function object for HW2 #3. */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // it seems that t is not needed, since f_spring is not explicitly 
    // dependent on t
    (void)t;
    Point f_grav = Point(0,0,-grav*n.value().mass);
    return f_grav;
  }
};

/** Force function object for HW2 #3. */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // it seems that t is not needed, since f_spring is not explicitly 
    // dependent on t
    (void)t;
    Point f_spring = Point(0);
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      f_spring += -e.value().K * (e.length()-e.value().L)/e.length() * 
        (e.node1().position()-e.node2().position());
    }
    return f_spring;
  }
};

/** Force function object for HW2 #3. */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // it seems that t is not needed, since f_spring is not explicitly 
    // dependent on t
    (void)t;
    Point f_damping = -c_ * n.value().vel;
    return f_damping;
  }

  DampingForce (double c) {
    c_ = c;
  }

  DampingForce(DampingForce& df) : c_(df.c_){}

  double c_;
};

struct ZeroForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t; (void)n;
    return Point(0);
  }
};


// /** just for self-learning purpose
// struct Problem3Force {
//   /** Return the force applying to @a n at time @a t.
//    */
//   template <typename NODE>
//   Point operator()(NODE n, double t) {
//     if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) return Point(0);
//     GravityForce gf;
//     MassSpringForce mf;
//     DampingForce df(c_);
//     return gf(n,t) + mf(n,t) + df(n,t);
//   }
//   Problem3Force (double c) {
//     c_ = c;
//   }
//   double c_;
// };

/** combined force function object for HW2 #3. */
template <typename F1, typename F2>       
struct make_combined_force {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n,t) + f2_(n,t);
  }
  make_combined_force (const F1& f1, const F2& f2) {
    f1_ = f1;
    f2_ = f2;
  }  
  F1 f1_;
  F2 f2_;
};

/** Force function object for HW2 #1. */
// template <typename F1, typename F2, typename F3=ZeroForce>       
// struct make_combined_force {
//   template <typename NODE>
//   Point operator()(NODE n, double t) {
//     return f1_(n,t) + f2_(n,t) + f3_(n,t);
//   }
//   make_combined_force (const F1& f1, const F2& f2, const F3& f3=ZeroForce()) {
//     f1_ = f1;
//     f2_ = f2;
//     f3_ = f3;
//   }  
//   F1 f1_;
//   F2 f2_;
//   F3 f3_;  
// };

/** combined constraints function object for HW2 #4. */
template <typename F1, typename F2>       
struct make_combined_constraint {

  template <typename G>
  void operator()(G& g, double t) {
    f1_(g,t);
    f2_(g,t);
    return;
  }
  make_combined_constraint (const F1& f1, const F2& f2) {
    f1_ = f1;
    f2_ = f2;
  }  
  F1 f1_;
  F2 f2_;
};

struct PlaneConstraint {
  /** constraint to @a g at time @a t.
   */
  template <typename G>
  void operator()(G& g, double t) {
    (void)t;
    Point zn = Point(0,0,1);
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (n.position()[2] < z_) {
        n.position() *= Point(1,1,0);
        n.position() += z_ * zn;
        n.value().vel *= Point(1,1,0);
      }
    }
    return;
  }
  PlaneConstraint (double z) {
    z_ = z;
  }

  PlaneConstraint(PlaneConstraint& p) { z_ = p.z_; }
  double z_;
};

struct SphereConstraint {
  /** constraint to @a g at time @a t.
   */  
  template <typename G>
  void operator()(G& g, double t) {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (norm(n.position()-c_) < r_) {
        Point xd = n.position() - c_;
        Point R = xd / norm(xd);
        n.position() = c_ + r_ * R;
        n.value().vel = n.value().vel - dot(n.value().vel,R) * R;
      }
    }
    return;
  }
  SphereConstraint (Point c, double r) {
    c_ = c;
    r_ = r;
  }
  Point c_;
  double r_;
};

struct SphereConstraintRemove {
  /** constraint to @a g at time @a t.
   */
  template <typename G>
  void operator()(G& g, double t) {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (norm(n.position()-c_) < r_) {
        g.remove_node(n);
      }
    }
    return;
  }
  SphereConstraintRemove (Point c, double r) {
    c_ = c;
    r_ = r;
  }
  Point c_;
  double r_;
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

    // Diagonal edges: include as of HW2 #2 and beyond
    #if 1
        graph.add_edge(nodes[t[0]], nodes[t[3]]);
        graph.add_edge(nodes[t[1]], nodes[t[2]]);      
    #endif

    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

// set spring constant and rest length for HW2 #1
#if 0
  double K = 100;
  double L = (*graph.edge_begin()).length();    
#endif

// set to be true if use damping force
#if 1
  double c = 1.0/graph.size();
#endif

  // Set initial conditions for nodes
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;
    n.value().vel = Point(0);
    n.value().mass = 1.0/graph.size();
  }

// Set initial conditions for edges, as of HW2 #2 and beyond, eg #3, 4, 5
#if 1
  /* iterate over edge like following is problematic since
   * here we need to change value of every edge, both (node_a, node_b)
   * and (node_b, node_a)
  // for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
  //   auto e = *it;
  //   e.value().K = 100;
  //   e.value().L = e.length();
  // }
  */
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    auto n = *ni;
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      auto e = *ei;
      e.value().K = 100;
      e.value().L = e.length();
    }
  }
#endif        

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
      // reduce dt to half to improve stability
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        
        // for Problem 1

      #if 0
        symp_euler_step(graph, t, dt, Problem1Force(K, L));
      #endif

        // for Problem 2
      #if 0
        symp_euler_step(graph, t, dt, Problem2Force());
      #endif

        // for Problem 3 without damping force
      #if 0
        symp_euler_step(graph, t, dt, 
          make_combined_force<GravityForce, DampingForce>(GravityForce(),DampingForce(c)));
      #endif

      #if 0
        symp_euler_step(graph, t, dt, 
          make_combined_force<GravityForce, MassSpringForce, DampingForce>
          (GravityForce(),MassSpringForce(),DampingForce(c)));
      #endif


      // for Problem 4 PlaneConstraint
      #if 0
        symp_euler_step(graph, t, dt, 
          make_combined_force<GravityForce, MassSpringForce>(GravityForce(),MassSpringForce()),
          PlaneConstraint(-0.75));
      #endif

      // for Problem 4 SphereConstraint
      #if 0
        symp_euler_step(graph, t, dt, 
          make_combined_force<GravityForce, MassSpringForce>(GravityForce(),MassSpringForce()),
          SphereConstraint(Point(0.5,0.5,-0.5),0.15));
      #endif

      // for Problem 4 SphereConstraintRemove
      #if 1
        viewer.clear();
        node_map.clear();
        symp_euler_step(graph, t, dt, 
          make_combined_force<GravityForce, MassSpringForce>(GravityForce(),MassSpringForce()),
          SphereConstraintRemove(Point(0.5,0.5,-0.5),0.15));
          viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
          viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

      #endif

      #if 0
        symp_euler_step(graph, t, dt, Problem2Force());
      #endif                

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
