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
double c; //damping coefficient

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double Kij;       //< Edge spring constant
  double Lij;       //< Edge rest length
  EdgeData() : Kij(100), Lij(0) {}
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

    if(!(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))){
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    }
  }

  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    if(!(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))){
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    }
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
    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }

    Point xi = n.position();

    Point fgrav = n.value().mass*Point(0,0,-grav);
    Point fspring = Point(0,0,0);
    Point ftotal = Point(0,0,0);

    // Calculate fspring
    for(auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){
      Edge e = *ei;
      Node neighbor = e.node2();                               // adjacent node to n
      Point xj = neighbor.position();

      fspring += (-e.value().Kij)*((xi-xj)/norm(xi-xj))*(norm(xi-xj)-e.value().Lij);
      //fspring += (-K)*((xi-xj)/norm(xi-xj))*(norm(xi-xj)-L);
    }

    ftotal = fspring + fgrav;

    (void) t;
    return ftotal;
  }
};

/** A struct which acts as a functor to impose the force of gravity on
 * on a node.
 * @param n     Graph::Node on which the force should be applied to
 * @param t     current time
 * @return      a point representing the vector force imposed by gravity on @a n
 *
 * @post        no elements within @a n are modified
 *
 */
struct GravityForce{
  /** Return the force applying to @a n at time @a t by gravity */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().mass*Point(0,0,-grav);
  }
};

/** A struct which acts as a functor to impose the force of the mass spring on
 * on a node.
 * @param n     Graph::Node on which the force should be applied to
 * @param t     current time
 * @return      a point representing the vector force imposed by the mass spring on @a n
 *
 * @post        no elements within @a n are modified
 *
 */
struct MassSpringForce{
  /** Return the force applying to @a n at time @a t by the mass-spring.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point xi = n.position();
    Point fspring = Point(0,0,0);

    for(auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){
      Edge e = *ei;
      Node neighbor = e.node2();                               // adjacent node to n
      Point xj = neighbor.position();

      fspring += (-e.value().Kij)*((xi-xj)/norm(xi-xj))*(norm(xi-xj)-e.value().Lij);
    }

    (void) t;
    return fspring;
  }
};

/** A struct which acts as a functor to impose the force of damping on
 * on a node.
 * @param n     Graph::Node on which the force should be applied to
 * @param t     current time
 * @return      a point representing the vector force imposed by damping on @a n
 *
 * @post        no elements within @a n are modified
 *
 */
struct DampingForce{
  /** Return the force applying to @a n at time @a t by damping.*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -c*n.value().vel;
  }
};

/** A struct which acts as a default force functor and returns a zero valued force (Point).
 * @param n     Graph::Node
 * @param t     current time
 * @return      a Point representing a zero vector or force which does not
 *              affect a node when added to it.
 *
 * @pre         a valid node or time is NOT required
 * @post        no elements within @a n are modified
 *
 */
struct nullForce{
  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    (void) n;
    return Point(0,0,0);
  }
};

/** A struct which acts as a 'combination functor', containing multiple functors.
  * @tparam F1, F2, F3   typename for forces @a f1, @a f2, @a f3
  * @param  f1, f2, f3   objects representing the forces to be applied
  */
template <typename F1, typename F2, typename F3>
struct combined_force{
  F1 f1;
  F2 f2;
  F3 f3;

  /** A constructor which creates a combination functor consisting of multiple
   * functors.
   * @tparam F1, F2, F3   typename for forces @a f1, @a f2, @a f3
   * @param  first, second, third   objects representing the forces to be applied
   * @return      a combined_force object containing the desired functors
   *
   */
  combined_force(F1 first, F2 second, F3 third) : f1(first), f2(second), f3(third){}

  /** operator() function which sums the outputs of the three force functors.
   * @param  n   Graph::Node object for which the forces are to be applied to
   * @param  t   current time
   * @return     a Point representing the combined force acting on the node
   */
  template <typename NODE>
  Point operator()(NODE n, double t){
    return f1(n, t) + f2(n, t) + f3(n, t);
  }
};

/** A function which creates and returns a combined_force object. The purpose is to
 * create a single functor representing multiple functors.
 * @param  f1, f2, f3   constructors for the desired forces
 * @tparam F1, F2, F3   typename for forces @a f1, @a f2, @a f3. The default
                        typename is set to nullForce so all three arguments do not
                        have to be passed into the function.
 * @return              an initialized combined_force object with the desired functors
 */
template <typename F1 = nullForce, typename F2 = nullForce, typename F3 = nullForce>
combined_force<F1, F2, F3> make_combined_force(F1 f1 = nullForce(), F2 f2 = nullForce(), F3 f3 = nullForce()){
  return combined_force<F1, F2, F3>(f1, f2, f3);
}


/** A struct which acts as a functor to impose a constraint of a plane at z = -0.75
 * on the entire simulation. If n.position() < z, the value of the node is 'corrected'.
 * @param g     Graph on which the constraint should be applied to
 * @param t     current time
 *
 * @post        If n.position() < z (where z = -0.75),
 *                1) the z component of n.position() is set to -0.75
 *                2) the z component of n.value().vel is set to 0
 */
struct PlaneConstraint{
  template <typename G>
  void operator()(G& g, double t){
    double z = -0.75;

    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      auto n = *ni;
      double zpos = dot(n.position(), Point(0,0,1));

      if(zpos < z){
        n.position() += Point(0,0,z-zpos);
        double zval = dot(n.value().vel, Point(0,0,1));
        n.value().vel += Point (0,0,-zval);
      }
    }

    (void) t;
    return;
  }
};

/** A struct which acts as a functor to impose a constraint of a sphere centered at
 * c = (0.5, 0.5, -0.5) and radius r = 0.15. If norm(n.position() - c) < r, the
 * value of the node is 'corrected'.
 * @param g     Graph on which the constraint should be applied to
 * @param t     current time
 *
 * @post        If norm(n.position() - c) < r:
 *                1) n.position() = c + r*(n.position() - c)/norm(n.position() - c)
 *                2) n.value().vel = n.value().vel - (dot(n.value().vel, Ri)*Ri), where
 *                   Ri = (n.position()-c)/norm(n.position()-c)
 *
 */
int check = 0;
struct SphereConstraint{
  template <typename G>
  void operator()(G& g, double t){
    Point c = Point(0.5,0.5,-0.5);    // center of sphere
    double r = 0.15;                  // sphere radius

    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      auto n = *ni;
      Point xi = n.position();


      if(norm(xi - c) < r){           //violates constraint
        // Fix position
        n.position() = c + r*(xi - c)/norm(xi - c);

        // Fix velocity
        Point Ri = (xi-c)/norm(xi-c);
        n.value().vel = n.value().vel - (dot(n.value().vel, Ri)*Ri);
      }

    }

    (void) t;
    return;
  }
};

/** A struct which acts as a functor to impose a constraint of a sphere centered at
 * c = (0.5, 0.5, -0.5) and radius r = 0.15. If norm(n.position() - c) < r, the
 * value of the node is 'corrected'.
 * @param g     Graph on which the constraint should be applied to
 * @param t     current time
 *
 * @post        If norm(n.position() - c) < r, remove node
 *
 */
struct SphereRemoveConstraint{
  template <typename G>
  void operator()(G& g, double t){
    Point c = Point(0.5,0.5,-0.5);    // center of sphere
    double r = 0.15;                  // sphere radius

    for(auto ni = g.node_begin(); ni != g.node_end(); ++ni){
      auto n = *ni;
      Point xi = n.position();

      if(norm(xi - c) < r){           //violates constraint
        g.remove_node(ni);
      }
    }

    (void) t;
    return;
  }
};

/** A struct which acts as a default constraint functor. The purpose is to act as
 * a placeholder for the make_combined_constraint function
 * @param g     Graph
 * @param t     current time
 *
 * @post        no elements within @a g are modified
 *
 */
struct nullConstraint{
  template <typename G>
  void operator()(G& g, double t){
    (void) t;
    (void) g;
    return;
  }
};


/** A struct which acts as a 'combination functor', containing multiple functors.
  * @tparam C1, C2, C3   typename for constraints @a c1, @a c2, @a c3
  * @param  c1, c2, c3   objects representing the constraints to be applied
  */
template <typename C1, typename C2, typename C3>
struct combined_constraint{
  C1 c1;
  C2 c2;
  C3 c3;

  /** A constructor which creates a combination constraint consisting of multiple
   * functors.
   * @tparam C1, C2, C3   typename for constraints @a c1, @a c2, @a c3
   * @param  first, second, third   objects representing the forces to be applied
   * @return      a combined_constraint object containing the desired functors
   *
   */
  combined_constraint(C1 first, C2 second, C3 third) : c1(first), c2(second), c3(third){}

  /** operator() function which enforces each of the three constraint functors.
   * @param  g   Graph object for which the constraints are to be applied to
   * @param  t   current time
   */
  template <typename G>
  void operator()(G& g, double t){
    c1(g, t);
    c2(g, t);
    c3(g, t);
    return;
  }
};

/** A function which creates and returns a combined_constraint object. The purpose is to
 * create a single functor representing multiple functors.
 * @param  c1, c2, c3   constructors for the desired constraint functors
 * @tparam C1, C2, C3   typename for constraint functors @a c1, @a c2, @a c3. The default
                        typename is set to nullConstraint so all three arguments do not
                        have to be passed into the function.
 * @return              an initialized combined_constraint object with the desired functors
 */
template <typename C1 = nullConstraint, typename C2 = nullConstraint, typename C3 = nullConstraint>
combined_constraint<C1, C2, C3> make_combined_constraint(C1 c1 = nullConstraint(), C2 c2 = nullConstraint(), C3 c3 = nullConstraint()){
  return combined_constraint<C1, C2, C3>(c1, c2, c3);
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
  c = 1.0/double(graph.num_nodes());

  // set mass of each node
  for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
    Node n = *ni;
    n.value().mass = 1.0/double(graph.num_nodes());
  }

  // set Lij of each edge
  for(auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei){
    Edge e = *ei;
    e.value().Lij = norm(e.node1().position() - e.node2().position());
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
        //symp_euler_step(graph, t, dt, Problem1Force());
        GravityForce f1 = GravityForce();
        MassSpringForce f2 = MassSpringForce();
        DampingForce f3 = DampingForce();
        PlaneConstraint c1 = PlaneConstraint();
        SphereConstraint c2 = SphereConstraint();
        SphereRemoveConstraint c3 = SphereRemoveConstraint();
        symp_euler_step(graph, t, dt, make_combined_force(f1,f2,f3), make_combined_constraint(c3, c1));

        //clear the viewer's nodes and edges
        viewer.clear(); // added for HW 2 #5
        node_map.clear(); // added for HW 2 #5

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map); // added for HW 2 #5
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
