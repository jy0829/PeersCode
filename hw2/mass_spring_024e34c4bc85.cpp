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
    
    
// Citations/Collaboration:
// Worked with Ian Shaw, Jonathon Roth, and Amery Martin (HW0) directly
// Attended office hours with Amy Shoemaker (TA) (HW0)
// Referenced (and followed the structure closely) of proxy_example.cpp
// Also referenced public github reposititories for previous students in 
// Harvard CS207: Lan, Tian; Piasevoli, Filip; Tran, Dustin; Zacarias, Lisa; 
// **Code was not copied direcly, only referenced for ideas when I was stuck


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Add addtional EdgeData structure (ref. CS207)
struct EdgeData{ 
  // Resting length
    double L;
  // Spring Constant
    double K;
};



// Define the Graph type with EdgeData included
using GraphType = Graph<NodeData, EdgeData>; 
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

// Generate plane constraint so nodes do not fall to infinity
struct PlaneConstraint {
    void operator() (GraphType& g, double){
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            Node node = (*it);
            if (dot(node.position(), Point(0, 0, 1)) < -0.75){
                node.position().elem[2] = -0.75;
                node.value().vel.elem[2] = 0;
            }
        }
    }
};

// Generate sphere constraint 
struct SphereConstraint { 
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    void operator()(GraphType& g, double){
        for (auto it = g.node_begin();  it != g.node_end(); ++it){
            Node node = (*it);
            if (norm(node.position()-c) < r){
                Point R = (node.position()-c)/norm(node.position()-c);
                node.position() = c + R*r;
                node.value().vel = node.value().vel - dot(node.value().vel, R)*R;
            }
        }
    }
};

// Generate sphere removal constraint
struct SphereRemoveConstraint{
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    void operator()(GraphType& g, double){
        auto it = g.node_begin();
        while(it != g.node_end()){
            if (norm((*it).position() - c) < r){
                it = g.remove_node(it);
            }
            else{
                ++it;
            }
        }
    }
};

// Combine two constraints
template<typename C1, typename C2>
struct CombinedConstraint{
    C1 cons1;
    C2 cons2;
    CombinedConstraint(C1 c1=C1(), C2 c2=C2()):cons1(c1), cons2(c2){}
    void operator()(GraphType& g, double){
        cons1(g,0);
        cons2(g,0);
    }
};

template<typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1=C1(), C2 c2=C2()){
    return CombinedConstraint<C1, C2>(c1, c2);
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
    CombinedConstraint<SphereRemoveConstraint, PlaneConstraint> constraint;
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0, 0, 0) && n.position() != Point(1,0,0)){
        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }
  }
  
  constraint(g,0);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)){
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
  //template <typename NODE> //<-------------------------------------------------------------
  Point operator()(Node n, double) {
    // HW2 #1: YOUR CODE HERE
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
    //return Point(0);
      if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
          return Point(0,0,0);
      Point spring = Point(0,0,0);
      Point gravity;
      gravity = Point (0,0,-grav)*n.value().mass;
    // Sum the spring forces
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
          //std::cout<< "Rest Length: " <<(*it).value().L << std::endl;
          spring += (*it).value().K * ((*it).node2().position()-(*it).node1().position()) / (*it).length() * ((*it).length()-(*it).value().L);  
      }
      return(spring + gravity);
  }
};
// Gravity force function object
struct GravityForce{
    Point operator()(Node n, double ){
        return (Point(0,0,-grav)*n.value().mass);
    }
};
// Spring force function object following format in Problem1Force
struct MassSpringForce{
    Point operator()(Node n, double){
        Point spring = Point(0,0,0);
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
            spring += (*it).value().K * ((*it).node2().position()-(*it).node1().position()) / (*it).length() * ((*it).length()-(*it).value().L); 
        }
        return spring;
    }
};

// Damping force function object
struct DampingForce{
    Point operator()(Node n, double ){
        return(-(c*n.value().vel));
    }
    static double c;
};
double DampingForce::c = 0;
// Function to combined two forces
template<typename F1, typename F2>
struct CombinedForce{
    F1 force1;
    F2 force2;
    CombinedForce(F1 f1=F1(), F2 f2=F2()): force1(f1), force2(f2){}
    Point operator()(Node n, double t){
        (void) t;
        return (force1(n,0)+force2(n,0));
    }
};

template<typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1 = F1(), F2 f2=F2()){
    return CombinedForce<F1, F2>(f1, f2);
}
// Function to combined three forces
template<typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>, F3> make_combined_force(F1 force1, F2 force2, F3 force3){
    return make_combined_force(make_combined_force(force1, force2), force3);
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
    //std::vector<Node> nodes;
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
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
        (*it).value().mass = float(1)/graph.size();
        (*it).value().vel = Point(0,0,0);
    }
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
        for (auto j = (*it).edge_begin(); j!= (*it).edge_end(); ++j){
            (*j).value().L = (*j).length();
            (*j).value().K = 100;
            //std::cout << "Rest Length: " << (*j).value().L << std::endl;
        }
    }
    DampingForce::c = float(1)/graph.num_nodes();
    

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
        symp_euler_step(graph, t, dt, Problem1Force());

        // Update viewer with nodes' new positions
//        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.set_label(t);
          
          viewer.clear();
          node_map.clear();
          viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
          viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

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
