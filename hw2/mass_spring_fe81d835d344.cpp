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

//for HW2 PART 1
//double L;
//double K=100;
double c;
/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** custom structure of data to store with Edges */
struct EdgeData{

  double L;
  double K;
  EdgeData() : L(0), K(100) {}
};


// Define the Graph type
//using GraphType = Graph<NodeData,EdgeData>;
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


//HW2 PART 3
struct GravityForce{

   //f_grav at current node n-a Point object
  template<typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    Point f_grav(0,0,-n.value().mass*grav);
    return f_grav;
  }

};

struct MassSpringForce{

  //f_spring-a Point object
  template<typename NODE>
  Point operator()(NODE n,double t){
      Point f_spring(0);
      for(auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){
        Edge e = *ei;
        Point p1 = e.node1().position();
        Point p2 = e.node2().position();

        f_spring = f_spring + (-e.value().K*(norm(p1-p2)-e.value().L)*(p1-p2))/norm(p1-p2);

      }
   (void) t;
   return f_spring;
  }
};

struct DampingForce{

  template<typename NODE>
  Point operator()(NODE n,double t){

    Point force_damp;
    force_damp = -c*n.value().vel;
    (void) t;
    return force_damp;

  }

};



template<typename F1, typename F2>
struct Combined_Force{

  F1 f1;
  F2 f2;
  template<typename NODE>
  Point operator()(NODE n, double t) {

    return f1(n,t) + f2(n,t);

  }



};

/**Return the combined forces of Spring and Gravity
 *@param[in] f1 Force one.
 *#param[in] f2 Force two.
 *@return the summed forces of @a f1 and @a f2.

 */

template<typename F1, typename F2,typename F3>
Combined_Force<Combined_Force<F1,F2>,F3> make_combined_force(F1 f1, F2 f2){

  return {f1,f2};

}

template<typename F1, typename F2,typename F3>
Combined_Force<Combined_Force<F1,F2>,F3>  make_combined_force(F1 f1,F2 f2, F3 f3){

   return {{f1,f2},f3};
}


//H5 PART 4

struct FixedConst{

  template<typename G>
  void operator()(G& g, double t){

    for(auto it = g.node_begin(); it != g.node_end(); ++it){

      auto n =*it;
      if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        n.position() = Point(0,0,0);
    }
  }

};

struct PlaneConst{

  template<typename G>
  void operator()(G& g, double t){


    for(auto it = g.node_begin(); it != g.node_end(); ++it){

      auto n =*it;
      Point p(0,0,-0.75);

      if(n.position().z <  p.z){
        n.position().z = p.z;
        n.value().vel.z = 0;
      }

    }
  }

};

struct SphereConst{

  template<typename G>
  void operator()(G& g, double t){

    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      Point c(0.5,0.5,-0.5);
      double r = 0.15;
      if(norm(n.position()-c) < r){
        double d1 = r - n.position().x;
        double d2 = r - n.position().y;
        double d3 = r - n.position().z;
        n.position().x = n.position().x + d1;
        n.position().y = n.position().y + d2;
        n.position().z = n.position().z + d3;
        Point R = (n.position()-c)/norm(n.position()-c);
        Point perp = n.value().vel - inner_prod(n.value().vel,R)*R;
        if(perp.x==0) n.value().vel.x=0;
        if(perp.y==0) n.value().vel.y=0;
        if(perp.z==0) n.value().vel.z=0;

      }
    }
  }

};

struct SphereConst2{

 template<typename G>
  void operator()(G& g, double t){

    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      Point c(0.5,0.5,-0.5);
      double r = 0.15;
      if(norm(n.position()-c) < r){
         g.remove_node(n);
      }
    }
  }



};


template<typename G, typename C1, typename C2>
void combine_const(G& g, double t, C1 c1, C2 c2){

c1(g,t);
c2(g,t);

}


template<typename G, typename C1, typename C2, typename C3>
void combine_const(G& g, double t, C1 c1, C2 c2, C3 c3){

c1(g,t);
c2(g,t);
c3(g,t);

}


template<typename G, typename C1, typename C2, typename C3, typename C4>
void combine_const(G& g, double t, C1 c1, C2 c2, C3 c3, C4 c4){

c1(g,t);
c2(g,t);
c3(g,t);
c4(g,t);
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
    Point f_grav(0,0,-n.value().mass*grav);

     Point f_spring(0);
      for(auto ei = n.edge_begin(); ei != n.edge_end(); ++ei){
        Edge e = *ei;
        Point p1 = e.node1().position();
        Point p2 = e.node2().position();

        f_spring = f_spring + (-e.value().K*(norm(p1-p2)-e.value().L)*(p1-p2))/norm(p1-p2);
      }
    return f_spring + f_grav;
  }
};

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

    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      continue;
    
    

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
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


  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
/*
  //set L to 1 edge length in graph g-HW2 PART 1
  Node first_node = *graph.node_begin();
  Edge first_edge = *first_node.edge_begin();
  L = first_edge.length();

*/
  //set the mass to all nodes and vel
  for(auto it = graph.node_begin(); it != graph.node_end();++it){
    auto n  = *it;
    n.value().mass = 1.0/graph.num_nodes();
    n.value().vel = Point(0,0,0);
    for(auto ei = n.edge_begin(); ei != n.edge_end();++ei){
      Edge e = *ei;
      e.value().L = e.length();
    }

  }

  c = 1.0/graph.num_nodes();

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
        
       symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce(),DampingForce()));
       // symp_euler_step(graph, t, dt, Problem1Force());
       // combine_const(graph,t,PlaneConst(),SphereConst());
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
