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

struct EdgeData {
  double K; 
  double L;
  EdgeData(double K, double L): K(K), L(L){}
  EdgeData() : K(100.0), L(1.0){}
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
//    if(n.position() != Point(0,0,0) and n.position() != Point(1,0,0)){
       n.value().vel += force(n, t) * (dt / n.value().mass);
 //   }
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

    if(n.position() == Point(0,0,0) or n.position() == Point(1,0,0))
      return Point(0,0,0);


    Point f_grav = n.value().mass* Point(0,0,-grav);
    Point f_spring = Point(0);
    for( auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){
      Edge e = *ei;
      Node n2 = e.node2(); //Adjacent node 
      Point diff = n.position() - n2.position();
      f_spring -=  e.value().K* (diff/(double) norm(diff))* (norm(diff) - e.value().L);

    }
    return f_grav + f_spring; 


  }
};

/** Functor that represents gravity
 */

struct GravityForce{
  template <typename NODE>
 /* Return the gravity force applying to node @a n at time @a t
 */
  Point operator()(NODE n, double t) {
    return n.value().mass* Point(0,0,-grav);
  }
};

struct MassSpringForce{
  template <typename NODE>
 /* Return the mass-spring force applying to node @a n at time @a t
 */

  Point operator()(NODE n, double t) {
    Point f_spring = Point(0);
    for( auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){
      Edge e = *ei;
      Node n2 = e.node2(); //Adjacent node 
      Point diff = n.position() - n2.position();
      f_spring -=  e.value().K* (diff/(double) norm(diff))* (norm(diff) - e.value().L);

    }
    return f_spring; 
  }
};

struct DampingForce{
  double c; //Constant c for the force
  DampingForce(double c_) : c(c_){}
  template <typename NODE>
 /* Return the damping force applying to node @a n at time @a t
 */

  Point operator()(NODE n, double t) {
    return -n.value().vel*c;
  }
};
  

template <typename F_1, typename F_2>
struct pair_force{
/* Structure to take into account a combination of two forces
*/
  F_1 f1;
  F_2 f2;
  pair_force(F_1 f1_, F_2 f2_) : f1(f1_), f2(f2_){} 
/* Force of @a f1 + @a f2 applied to node @a n at time @a t
*/
  template <typename NODE>
  Point operator()(NODE n, double t){
    return f1(n,t) + f2(n,t);
  }

};

/* Return the combination of two forces 
*/ 
template <typename F_1, typename F_2>
pair_force<F_1, F_2> make_combined_force(F_1 f1, F_2 f2){
  return pair_force<F_1, F_2>(f1,f2);
}


template <typename F_1, typename F_2, typename F_3>
struct three_force{
/* Structure to take into account a combination of three forces
*/
  F_1 f1;
  F_2 f2; 
  F_3 f3;
  three_force(F_1 f1_, F_2 f2_, F_3 f3_) : f1(f1_), f2(f2_), f3(f3_){} 
/* Force of @a f1, @a f2 and @a f3 applied to node @a n at time @a t
*/
  template <typename NODE>
  Point operator()(NODE n, double t){
    return f1(n,t) + f2(n,t) + f3(n,t);
  }

};

/* Return the combination of three forces 
*/

template <typename F_1, typename F_2, typename F_3>
three_force<F_1, F_2, F_3> make_combined_force(F_1 f1, F_2 f2, F_3 f3){
  return three_force<F_1, F_2, F_3>(f1,f2,f3);
}

/*Set the constraint regarding to a certain point
 */


struct Constraint_Point{
  Point p;
  Constraint_Point(Point p_) : p(p_){}

 /*For all nodes that are in the position @a p we set the velocity to 0
  */ 
  void operator()(GraphType& g, double t){
    for(auto n_it = g.node_begin(); n_it != g.node_end(); ++n_it){
      Node n = *n_it;
      if(n.position() == p){
        n.value().vel = Point(0);
      }
    }

  } 
  

};

struct Constraint_Plane{

 /*For all nodes that are in the plane  we change them
  */ 

  void operator()(GraphType& g, double t){
    for(auto n_it = g.node_begin(); n_it != g.node_end(); ++n_it){
      Node n = *n_it; 
      if(n.position().z < -0.75){
        n.value().vel.z = 0; 
        n.position().z = -0.75; // Nearest point in the plane
      }

    }
  }
};

struct Constraint_Sphere1{

 /*For all nodes that are in the sphere we change them
  */ 

  void operator()(GraphType& g, double t){
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
    for(auto n_it = g.node_begin(); n_it != g.node_end(); ++n_it){
      Node n = *n_it; 
      Point Ri = (n.position() - c)/norm(n.position() -c);
      if(norm(n.position() - c) < r){
        n.value().vel -= dot(n.value().vel, Ri)*Ri; 
        n.position() = c + r*Ri; // Nearest point in the plane
      }

    }
  }
};

struct Constraint_Sphere2{

 /*For all nodes that are in the sphere we change them
  */ 

  void operator()(GraphType& g, double t){
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
    for(auto n_it = g.node_begin(); n_it != g.node_end(); ++n_it){
      Node n = *n_it; 
      if(norm(n.position() - c) < r){
        g.remove_node(n);
      }

    }
  }
};

template <typename C_1, typename C_2>
struct pair_constraint{
/* Structure to take into account a combination of two constraints
*/
  C_1 cst1;
  C_2 cst2;
  pair_constraint(C_1 c1, C_2 c2) : cst1(c1), cst2(c2){} 
/* Constraint of @a cst1 + @a cst2 applied to graph @a g at time @a t
*/
  void operator()(GraphType& g, double t){
    cst1(g,t); 
    cst2(g,t);
  }

};

/* Return the combination of two constraints 
*/ 
template <typename C_1, typename C_2>
pair_constraint<C_1, C_2> make_combined_constraint(C_1 cst1, C_2 cst2){
  return pair_constraint<C_1, C_2>(cst1,cst2);
}


template <typename C_1, typename C_2, typename C_3>
struct three_constraint{
/* Structure to take into account a combination of three forces
*/
  C_1 cst1;
  C_2 cst2; 
  C_3 cst3;
  three_constraint(C_1 c1, C_2 c2, C_3 c3) : cst1(c1), cst2(c2), cst3(c3){} 
/* Constraint of @a cst1, @a cst2, @a cst3 applied to graph @a g at time @a t
*/
  void operator()(GraphType& g, double t){
    cst1(g,t); 
    cst2(g,t);
    cst3(g,t);
  }

};
template <typename C_1, typename C_2, typename C_3>
three_constraint<C_1, C_2,C_3> make_combined_constraint(C_1 cst1, C_2 cst2, C_3 cst3){
  return three_constraint<C_1, C_2,C_3>(cst1,cst2,cst3);
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

  int N = graph.num_nodes();
//  int cmpt_n = 0;
//  int cmpt_e = 0;
  std::cout << 1.0/N<< std::endl;
//Initialize the node data 
  for( auto n_it = graph.node_begin(); n_it!= graph.node_end(); ++n_it){
    Node n = *n_it;
    n.value().mass = 1.0/N;
    n.value().vel = Point(0);
    for(auto e_it = n.edge_begin(); e_it != n.edge_end(); ++e_it){
      Edge e = *e_it;
      e.value() = EdgeData(100.0, e.length());
    } 
  }

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;
//  std::cout << cmpt_n << " " << cmpt_e << std::endl;

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
      auto total_force = make_combined_force(GravityForce(),MassSpringForce(), DampingForce(1.0/graph.num_nodes()));
      auto const1 = Constraint_Point(Point(0,0,0));
      auto const2 = Constraint_Point(Point(1,0,0));
      auto cst = make_combined_constraint(const1, const2); 
      auto cst_plane = Constraint_Plane();
      auto cst_sphere = Constraint_Sphere1();

      auto total_constraint = make_combined_constraint(cst,cst_plane, cst_sphere);

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, total_force);
//        cst(graph,t);
        total_constraint(graph,t);
       // symp_euler_step(graph, t, dt, Problem1Force());

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
