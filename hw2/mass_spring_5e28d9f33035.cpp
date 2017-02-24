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

/** HW2: Custom structure of data to store with Edges */
struct EdgeData {
  double length;   //< Edge length
  double K;        //< Spring constant
  EdgeData() : length(0), K(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using idx_type = unsigned;

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
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      if ( n.position() != Point(0,0,0) && n.position() != Point(1,0,0) ) {
          n.position() += n.value().vel * dt;
      }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/* Indices of corner nodes - those with position() == Point(0,0,0)
   or Point(1,0,0) */
idx_type corner1_idx;
idx_type corner2_idx;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g           Graph
 * @param[in]     t           The current time (useful for time-dependent forces)
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @param[in]     constraint  Constraint object which checks for violations and
 *                            updates node position/velocity if necessary.
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a constraint object called as @a constraint(@a g, @a t),
 *           where @a g is the graph and @a t is the current time.
 *           @a constraint checks if a condition is violated at time @a t
 *           and if so, updates a node's position and velocity at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      if ( n.position() == Point(0,0,0) ) {
          corner1_idx = n.index();
      }
      else if ( n.position() == Point(1,0,0) ) {
          corner2_idx = n.index();
      }
      n.position() += n.value().vel * dt;
  }

  /* Check if constraints are violated */
  constraint(g, t);

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
   * model that by returning a zero-valued force.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {

    /** Prevent cloth from falling to infinity */
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        return Point(0,0,0);

    /** Calculate spring force */
    Point xi = n.position();
    Point F_spring(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        Point xj = (*it).node2().position();
        F_spring += -K_ * (xi - xj) * ( (*it).length() - L_ ) / (*it).length();
    }

    (void) t;

    return F_spring + Point(0,0, n.value().mass * -grav);
  }

   /** HW2: Constructor */
   Problem1Force(double K, double L)  : K_(K), L_(L) {}

   /** HW2: Spring constant */
   double K_;

   /** HW2: Spring rest-length */
   double L_;

};


/** Force function object for HW2 #2. */
struct Problem2Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE

    /** Prevent cloth from falling to infinity */
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        return Point(0,0,0);
    }

    /** Calculate spring force */
    Point xi = n.position();
    Point F_spring(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        Point xj = (*it).node2().position();
        double Lij = (*it).value().length;
        double Kij = (*it).value().K;
        F_spring += -Kij * (xi - xj) * ( (*it).length() - Lij ) / (*it).length();
    }

    (void) t;

    return F_spring + Point(0,0, n.value().mass * -grav);
  }

   /** HW2: Constructor */
   Problem2Force() {}

};

/* Gravity force function object for HW2 #3. */
struct GravityForce {

  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    return Point(0,0, n.value().mass * -grav);
  }

   /* Constructor */
   GravityForce() {}

};

/** Spring force function object for HW2 #3. */
struct MassSpringForce {

  /** Return the spring force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {

      Point xi = n.position();

      Point F_spring(0,0,0);

      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
          Point xj = (*it).node2().position();
          double Lij = (*it).value().length;
          double Kij = (*it).value().K;
          F_spring += -Kij * (xi - xj) * ( (*it).length() - Lij ) / (*it).length();
      }

      (void) t;
      return F_spring;
  }

  /* Constructor */
  MassSpringForce() {}

};

/** Damping force function object for HW2 #3. */
struct DampingForce {

  /** Return the damping force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
      (void) t;
      return -c_ * n.value().vel;
  }

  /* Constructor */
  DampingForce(const double c) : c_(c) {}

  private:
      /* Damping coefficient */
      const double c_;

};

/** Struct to combine two forces */
template <typename F1, typename F2>
struct combined_force {

  template <typename NODE>
  Point operator()(NODE n, double t) {
      return f1_(n, t) + f2_(n, t);
  }

  /* Constructor */
  combined_force(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}

  /* Force 1 and Force 2 */
  F1 f1_;
  F2 f2_;

};

/** Combine two forces */
template <typename F1, typename F2>
combined_force<F1, F2> make_combined_force(F1 f1, F2 f2) {
    return combined_force<F1,F2>(f1, f2);
}

/** Combine three forces */
template <typename F1, typename F2, typename F3>
combined_force<combined_force<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
    return combined_force<combined_force<F1,F2>,F3>(combined_force<F1,F2>(f1, f2), f3 );
}

/** Constant constraint for HW2 #4 */
struct ConstantConstraint {

  /* Updates the position of two globally defined corner nodes */
  void operator()(GraphType& g, double t) {
      g.node(corner1_idx).position() = Point(0,0,0);
      g.node(corner2_idx).position() = Point(1,0,0);
      (void) t;
  }

  /* Constructor */
  ConstantConstraint() {}

};

/** Plane constraint for HW2 #4 */
struct PlaneConstraint {

  void operator()(GraphType& g, double t) {
      double z = -0.75;
      for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
          if ( (*ni).position().z < z ) {
              /* Set z-component of position to 0 */
              (*ni).position().z = z;
              /* Set z-component of velocity to 0 */
              (*ni).value().vel.z = 0;
          }
      }
      (void) t;
  }

  /* Constructor */
  PlaneConstraint() {}

};

/** Sphere constraint for HW2 #4 */
struct SphereConstraint4 {

  void operator()(GraphType& g, double t) {
      for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
          Point pos = (*ni).position();
          double dist_from_center = norm(pos - center_);
          Point Ri = (pos - center_) / dist_from_center;
          if ( dist_from_center < radius_ ) {
              /* Set position to the nearest point on sphere's surface */
              (*ni).position() = center_ + radius_*Ri;
              /* Set component of velocity normal to the sphere's surface to 0 */
              (*ni).value().vel -= ( dot( (*ni).value().vel, Ri) * Ri );
          }
     }
    (void) t;
  }

  /* Constructor */
  SphereConstraint4(const Point center, double radius) :
      center_(center), radius_(radius) {}

  /* Sphere center */
  const Point center_;
  /* Sphere radius */
  double radius_;

};

/** Sphere constraint for HW2 #5 */
struct SphereConstraint5 {

  void operator()(GraphType& g, double t) {
      auto ni = g.node_begin();
      auto end = g.node_end();
      while ( ni != end ) {
          if ( norm( (*ni).position() - center_ ) < radius_ ) {
              g.remove_node(ni);
          }
          else {
              ++ni;
          }
      }


      /*}
      for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
          if ( norm( (*ni).position() - center_ ) < radius_ ) {
              g.remove_node(ni);
          }
      }*/
      (void) t;
  }

  /* Constructor */
  SphereConstraint5(const Point center, double radius) :
      center_(center), radius_(radius) {}

  /* Sphere center */
  const Point center_;
  /* Sphere radius */
  double radius_;
};

/** Struct to combine three constraints */
template <typename C1, typename C2, typename C3>
struct combined_constraint {

  void operator()(GraphType& g, double t) {
      c1_(g, t);
      c2_(g, t);
      c3_(g, t);
  }

  /* Constructor */
  combined_constraint(C1 c1, C2 c2, C3 c3) :
    c1_(c1), c2_(c2), c3_(c3) {}

  /* Constraints */
  C1 c1_;
  C2 c2_;
  C3 c3_;

};

/** Combine three constraints */
template <typename C1, typename C2, typename C3>
combined_constraint<C1, C2, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
    return combined_constraint<C1,C2,C3>(c1, c2, c3);
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.

  /* Initiliaze node mass, node velocity, edge rest lengths,
  and edge spring constants */
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {

      (*ni).value().mass = 1./graph.num_nodes();
      (*ni).value().vel = Point(0,0,0);

      for (auto ii = (*ni).edge_begin(); ii != (*ni).edge_end(); ++ii) {
          /* Set resting length to initial edge length */
          (*ii).value().length = (*ii).length();
          /* Set spring constant */
          (*ii).value().K = 100.;
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

      /** HW2: Define initial conditions for Problem1Force */
      double K = 100.;
      double L = (*graph.edge_begin()).length();
      (void) K; (void) L;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;

        /* HW2: Problem 1 */
        //symp_euler_step(graph, t, dt, Problem1Force(K, L));

        /* HW2: Problem 2 */
        //symp_euler_step(graph, t, dt, Problem2Force());

        /* HW2: Problem 3 */
        /*double c = 1./graph.num_nodes();
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),
          MassSpringForce(), DampingForce(c)) );*/

        /* HW2: Problem 4 & 5 */
        double c = 1./graph.num_nodes();
        const Point center(0.5,0.5,-0.5);
        double radius = 0.15;
        symp_euler_step(graph, t, dt,
          make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c)),
          make_combined_constraint(ConstantConstraint(), PlaneConstraint(),
            SphereConstraint5(center, radius)) );

        /* HW2 #5: Clear the viewer's nodes and edges */
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
