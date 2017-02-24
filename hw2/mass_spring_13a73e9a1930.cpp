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
// Spring constant
static constexpr double K = 100;
// Initial length
// static double L;

// Damping coefficient
static double C;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,double>;
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
    if (n.position() != Point(0,0,0) && n.position()!=Point(1,0,0)) {
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0,0,0) && n.position()!=Point(1,0,0)) {
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  // Call the passed constraint
  constraint(g);

  return t + dt;
}

/** struct GravityForce
 * @brief Returns the force from gravity applying to @a n at time @a t.
 * @pre @a n must be a valid node of graph object
 * @param Node @a n, double @a t - input node and time
 * @return Point representing force
 */

struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    double mass_n = n.value().mass;
    Point f_gravity = mass_n*Point(0,0,-grav);
    return f_gravity;
  }
};

/** struct MassSpringForce
 * @brief Returns the spring force applying to @a n at time @a t.
 * @pre @a n must be a valid node of graph object
 * @param Node @a n, double @a t - input node and time
 * @return Point representing force
 */

struct MassSpringForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point f_spring;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      /* Initialize edge value to resting edge length */
      if (t==0) {
        (*it).value()=(*it).length();
      }

      /* Calculate damping force */
      double eucDiff;
      Point p1 = n.position();
      Point p2 = (*it).node2().position();
      eucDiff = norm(p1 - p2);
      double l = (*it).value();
      Point foo = -K*(p1 - p2)/eucDiff * (eucDiff - l);
      f_spring = f_spring + foo;
    }
    return f_spring;
  }
};

/** struct GravityForce
 * @brief Returns the damping force applying to @a n at time @a t.
 * @pre @a n must be a valid node of graph object @a g
 * @param Node @a n, double @a t - input node and time
 * @return Point representing force
 *
 */

struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point vel_n = n.value().vel;
    Point f_damping = -C*vel_n;
    return f_damping;
  }
};

/** struct NullForce
 * @brief Returns a force of zero
 * @param Node @a n, double @a t - input node and time
 * @return Point representing force
 */
struct NullForce {
   template <typename NODE>
   Point operator()(NODE n, double t) {
     (void) t;
     (void) n;
     return Point(0,0,0);
   }
};

/** @struct make_combined_force
    @brief returns a point which is the sum of the input forces
    @param F1, F2, F3 - the input forces (F3 is null if not specified)
    @param Node @a n, double @a t - input node and time
    @return object of type Point which is the resulting force
*/
template <typename F1 = NullForce, typename F2 = NullForce, typename F3 = NullForce>
struct make_combined_force {
  /** Return a combination of forces
   * */
  F1 f1_;
  F2 f2_;
  F3 f3_;

  make_combined_force(F1 f1 = NullForce(), F2 f2 = NullForce(), F3 f3 = NullForce()) {
     f1_ = f1;
     f2_ = f2;
     f3_ = f3;
  }

  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point p1 = f1_(n,t);
    Point p2 = f2_(n,t);
    Point p3 = f3_(n,t);
    return ( p1 + p2 + p3 );
  }
};

 /** Implements a planar constraint on the positions and velocities of the
  * graph nodes
  * @param Graph @a g - input graph
  * @post Changes position and velocities of all nodes which violate constraint
  */
struct PlaneConstraint {
  template <typename G>
  void operator()(G& g) {
    double z = -0.75;
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto n = *ni;
      if (inner_prod(n.position(), Point(0,0,1)) < z) {
        // Fix the position
        Point p = n.position();
        n.position() = Point(p.x,p.y,z);

        // Fix the velocity
        Point v = n.value().vel;
        n.value().vel = Point(v.x,v.y,0);
      }
    }
  }
};

 /** Implements a spherical constraint on the positions and velocities of the
  * graph nodes using fixed sphere center and radius
  * @param Graph @a g - input graph
  * @post Changes position and velocities of all nodes which violate constraint
  */
struct SphereConstraint {
  double r = 0.15;
  Point C = Point(0.5,0.5,-0.5);

  template <typename G>
  void operator()(G& g) {

    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto n = *ni;
      Point P = n.position();
      if (norm(P-C) < r) {

        // Fix the position
        Point d = P - C;
        Point E = C + d*(r/norm(d));
        n.position() = E;

        // Fix the velocity
        Point v0 = n.value().vel;
        Point R_i = (E-C) / norm(E-C);
        Point v1 = v0 - inner_prod(v0,R_i)*R_i;
        n.value().vel = v1;
      }

    }
  }
};

 /** Implements a spherical constraint on the positions and velocities of the
  * graph nodes
  * @param Graph @a g - input graph
  * @post Removes all nodes which violate constraint
  */
struct RemoveSphereConstraint {
  double r = 0.15;
  Point C = Point(0.5,0.5,-0.5);

  template <typename G>
  void operator()(G& g) {
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      auto n = *ni;
      Point P = n.position();
      if (norm(P-C) < r) {
        g.remove_node(n);
      }

    }
  }
};

/** Imposes no constraints on graph @a g
 * @param - input graph @a g
 */
struct NullConstraint {

   template <typename G>
   void operator()(G& g) {
     (void) g;
   }
};

/** Return a combination of constraints
 */
template <typename C1 = NullConstraint, typename C2 = NullConstraint>
struct make_combined_constraint {
   C1 c1_;
   C2 c2_;

   make_combined_constraint(C1 c1 = NullConstraint(), C2 c2 = NullConstraint()) {
     c1_ = c1;
     c2_ = c2;
   }

   template <typename G>
   void operator()(G& g) {
     c1_(g);
     c2_(g);
   }
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
//#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  /* Initialize all node value data */
  NodeData nd;
  nd.mass = (double)1/graph.num_nodes();

  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    (*ni).value() = nd;
  }

  /* Set the damping constant */
  C = (double)1/graph.num_nodes();


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
        // symp_euler_step(graph, t, dt,
        //   make_combined_force<GravityForce,MassSpringForce,DampingForce>
        //   (GravityForce(),MassSpringForce(),DampingForce() ),
        //   make_combined_constraint<PlaneConstraint,SphereConstraint>
        //   (PlaneConstraint(),SphereConstraint()) );
        symp_euler_step(graph, t, dt,
          make_combined_force<GravityForce,MassSpringForce,DampingForce>
          (GravityForce(),MassSpringForce(),DampingForce() ),
          make_combined_constraint<RemoveSphereConstraint>
          (RemoveSphereConstraint()) );

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
