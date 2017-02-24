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
  double L;
  double K;
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
double symp_euler_step(G& g, double t, double dt, F force, C con = C()) {
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
  con(g,t);
  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
 // template <typename NODE>
  Point operator()(Node n, double t) {
    
    // HW2 #1: YOUR CODE HERE
    (void) t;  // silence compiler warnings
    
    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
        return Point(0,0,0);
    }
    
    double mass = n.value().mass;
    double K = 100;

    Point grav_F = Point(0,0,-grav*mass);
    Point spring_F = Point(0,0,0);
    
    Point dist;

    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        
        Point dist;
        if((*it).node1() == n){
            dist = (*it).node1().position() - (*it).node2().position();
        } else{
            dist = (*it).node2().position() - (*it).node1().position();
        }

        double length = (*it).length();
        double in_length = (*it).value().L;

        spring_F = spring_F + (-K)*(dist/length)*(length-in_length);
    } 
    
    return (spring_F + grav_F);
  }
};

struct GravityForce{
   Point operator()(Node n, double t) {
      (void) t;
      if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
        return Point(0,0,0);
      }

 
      double mass = n.value().mass;
      Point grav_F = Point(0,0,-grav*mass);
      
      return grav_F;
   }
 
};

struct MassSpringForce{
    Point operator()(Node n, double t) {
   
     (void) t;

        double K = 100;
        Point spring_F = Point(0,0,0);
    
        Point dist;

        for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        
            Point dist;
            if((*it).node1() == n){
                dist = (*it).node1().position() - (*it).node2().position();
            } else{
                dist = (*it).node2().position() - (*it).node1().position();
            }

            double length = (*it).length();
            double in_length = (*it).value().L;

            spring_F = spring_F + (-K)*(dist/length)*(length-in_length);
        }   
    
    return spring_F;        
    }
};

struct DampingForce{
    double con = 0;
    
    DampingForce(){
    }
    
    DampingForce(double con_input){
        con = con_input;
    }

    Point operator()(Node n, double t) {
        (void) t;
        Point damping_F = (-con)*n.value().vel;
        return damping_F;
    }
};

template<typename F1, typename F2>
struct TotalForce{
    F1 force_1;
    F2 force_2;

    TotalForce(){
    }
    
    TotalForce(F1 input1, F2 input2){
        force_1 = input1;
        force_2 = input2;
    }
    
    Point operator()(Node n, double t){
        (void) t;
        return (force_1(n,t)+force_2(n,t));
    }
};

template<typename F1, typename F2>
TotalForce<F1, F2> make_combined_force(F1 force1, F2 force2){
    return TotalForce<F1, F2>(force1, force2);
}

template<typename F1, typename F2, typename F3>
TotalForce<F1, TotalForce<F2, F3>> make_combined_force(F1 force1, F2 force2, F3 force3){
    return make_combined_force(force1,make_combined_force(force2, force3));
}

struct ZConstraint{
    double z = 0.75;
    void operator()(GraphType& g, double t){
        (void) t;
        for(auto it = g.node_begin(); it != g.node_end(); ++it){
        
            Node n = (*it);
            if(n.position().z < -z){
               n.position().z = -z;
               n.value().vel.z = 0;
            }
        } 
    }
    
};

struct SphereConstraint{
    double r = 0.15;
    Point center = Point(0.5, 0.5, -0.5);
    void operator()(GraphType& g, double t){
    
        (void) t;
        for(auto it = g.node_begin(); it != g.node_end(); ++it){

            Node n = (*it);
            Point dist = n.position()-center;
            double dis = norm(dist);

            if(dis == 0){
                n.position() = n.position() + Point(0,0,0.15);
                n.value().vel.z = 0;
            } else if(dis < r){
                n.position() = center+(dist*(r/dis));
                n.value().vel = n.value().vel - dot(n.value().vel, dist/dis)*(dist/dis);
            }
        }
    }
};

struct SphereConstraint2{
    double r = 0.15;
    Point center = Point(0.5, 0.5, -0.5);
    
    
    void operator()(GraphType& g, double t){
         (void) t;
         for (auto it = g.node_begin(); it != g.node_end(); ++it){
             Node n = (*it);
             Point dist = n.position()-center;
             double dis = norm(dist);

             if (dis < r){
                 g.remove_node(n);
             }
        }
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
#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  for(auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit){
      (*eit).value().L = (*eit).length();
  }
  
  for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit){
      (*nit).value().mass = (1.0/graph.num_nodes());
      (*nit).value().vel = Point(0,0,0);
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
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;
      //double c = 1.0/graph.num_nodes();

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        SphereConstraint con;
        //SphereConstraint2 con2;

        //auto totalF =  make_combined_force(GravityForce(), MassSpringForce(), DampingForce(c));
        symp_euler_step(graph, t, dt, Problem1Force(), con);
        
        // viewer.clear();
        //node_map.clear();
        //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);     
        //viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map); 
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
