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
  double c;
  NodeData() : vel(0), mass(1),c(0.001) {}
};

struct EdgeData {
  double K;
  double L;
  EdgeData() : K(100.0),L(0.01){}
};



// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

// 
double K = 100.0; 

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
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
    if(n.position() == Point(0) || n.position() == Point(1,0,0)){
        n.value().vel = Point(0);}
    else{
       n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }
  return t + dt;
}

template <typename G, typename C>
void updating(G& g, double t, C cons) {
   cons(g,t);
}

struct Problem1Force{
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    if(n.position() == Point(0) || n.position() == Point(1,0,0)) return Point(0);
    Point force = Point(0,0,-1)* grav * n.value().mass;
    for(auto iter = n.edge_begin() ; iter != n.edge_end(); ++iter){
       auto node_j = (*iter).node2();    
       double dist_i_j = norm(node_j.position() - n.position());
       force +=  K * (node_j.position() - n.position() ) / dist_i_j * (dist_i_j - (*iter).value().L);
    } 
    return force;
  }
};

struct MassSpringForce{
   template <typename NODE>
   Point operator()(NODE n, double t) {
    Point force = Point(0,0,0);
    for(auto iter = n.edge_begin() ; iter != n.edge_end(); ++iter){
       auto node_j = (*iter).node2();    
       double dist_i_j = norm(node_j.position() - n.position());
       force +=  K * (node_j.position() - n.position() ) / dist_i_j * (dist_i_j - (*iter).value().L);
    } 
    return force;
  }
};

struct GravityForce{
 template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force = Point(0,0,-1)* grav * n.value().mass;
    return force;
  } 
};

struct DampingForce{
 template <typename NODE>
  Point operator()(NODE n, double t) {
    Point force =  - 1 * n.value().c * n.value().vel;
    return force;
  }
};

/*class make_combined_force{
   std::vector<Force> forces;
public:
   template<typename... Args>
   make_combined_force(Args&&... args):forces(std::forward<Args>(args)...){ }
   template<typename NODE>
   Point operator()(NODE n,double t){
     return Point(0);  
   }
};
*/

struct constraint_1{
   template<typename GRAPH>
   void operator()(GRAPH& g,double t){
       for(auto iter = g.node_begin() ; iter != g.node_end(); ++iter){
          if((*iter).position().z < plane_z){
                (*iter).value().vel.z = 0;
                (*iter).position().z = plane_z;
          }
      }
   }
   double plane_z = -0.75;
};

struct constraint_2{
   template<typename GRAPH>
   void operator()(GRAPH& g,double t){
       for(auto iter = g.node_begin() ; iter != g.node_end(); ++iter){
          Point R = ((*iter).position() - c);
          double rn = norm(R);    
          if(rn < r){
             R  /= rn;
             (*iter).position() = R * r + c;
             (*iter).value().vel -=  inner_prod((*iter).value().vel,R)* R; 
          }
      }
   }
   Point c = Point(0.5,0.5,-0.5);
   double r = 0.15;
};

struct constraint_3{
   template<typename GRAPH>
   void operator()(GRAPH& g,double t){
       for(auto iter = g.node_begin() ; iter != g.node_end(); ++iter){
          Point R = ((*iter).position() - c);
          double rn = norm(R);    
          if(rn < r){
           g.remove_node((*iter));
          }
      }
   }
   Point c = Point(0.5,0.5,-0.5);
   double r = 0.15;
};

/*template<typename U1, typename U2>
struct make_combined_force{
   template<class T1,class T2>
   make_combined_force(T1&& t1_, T2&& t2_):t1(std::forward<T1>(t1_)),t2(std::forward<T2>(t2_)){}
   template<typename NODE>
   Point operator()(NODE n,double t){
     return t1(n,t) + t2(n,t); 
   }
   U1 t1;
   U2 t2;
};
*/

template<typename U1, typename U2, typename U3>
struct make_combined_force{
   template<class T1,class T2,class T3>
   make_combined_force(T1&& t1_, T2&& t2_, T3&& t3_):t1(std::forward<T1>(t1_)),t2(std::forward<T2>(t2_)),t3(std::forward<T3>(t3_)){}
   template<typename NODE>
   Point operator()(NODE n,double t){
     return t1(n,t) + t2(n,t) + t3(n,t); 
   }
   U1 t1;
   U2 t2;
   U3 t3;
};



//template<typename T1, typename T2>
//comb<T1,T2>& make_combined_force(T1 t1, T2 t2){
//   return comb<T1,T2>(t1,t2);
//};

//template<typename T1, typename T2, typename T3>
//comb<comb<T1,T2>,T3>& make_combined_force(T1 t1, T2 t2, T3 t3){
//     auto it = make_combined_force<T1,T2>(t1,t2);
//     return comb<comb<T1,T2>,T3>(it,t3);
//};

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

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
   for(auto iter = graph.node_begin(); iter != graph.node_end(); ++iter){
     (*iter).value().vel = Point(0,0,0);
     (*iter).value().mass =  1.0 / graph.num_nodes();
     (*iter).value().c = 1.0 / graph.num_nodes();
  }
  
 
 // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;
  for(auto iter = graph.edge_begin(); iter != graph.edge_end(); ++iter){
       (*iter).value().K = 100.0;
       (*iter).value().L = (*iter).length();
  }
  
   
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
        symp_euler_step(graph, t, dt, make_combined_force<GravityForce,MassSpringForce,DampingForce>(GravityForce(), MassSpringForce(), DampingForce()));
        //symp_euler_step(graph, t, dt, make_combined_force<GravityForce,MassSpringForce>(GravityForce(), MassSpringForce()));


        viewer.clear();
        node_map.clear();
        updating(graph,t,constraint_3());
        
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
