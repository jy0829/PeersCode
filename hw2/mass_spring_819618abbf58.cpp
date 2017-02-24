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
#include <iostream>

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
  double K_;       //< Spring const
  double L_;     // Rest length
  EdgeData() : K_(100), L_(1) {}
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force,C constraint ) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    //std::cout << n.value().vel<<std::endl;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    //std::cout << t << " " << n.position() <<std::endl;
  }

  constraint(g,t);

  return t + dt;
}

/* Gravity Force 
*/
struct GravityForce{
  template <typename NODE>
  Point operator()(NODE n,double t){
    return Point(0,0,-n.value().mass*grav);
    (void) t;
  }
  
};

/* Spring force
*/
struct MassSpringForce{
  template <typename NODE>
  Point operator()(NODE n,double t){
    Point spring(0,0,0); // For spring force
    for(auto iit=n.edge_begin();iit!=n.edge_end();++iit){
      auto K = (*iit).value().K_;
      auto L = (*iit).value().L_;
      //std::cout << "L_" << " "<< L << "legn" << " " << (*iit).length() << std::endl;
      double force_const = -(K/(*iit).length())*((*iit).length()-L);
      spring += force_const*((*iit).node1().position()-(*iit).node2().position());
    }
    return spring;
    (void) t;    
  }
  
};

/* Damping force
*/
struct DampingForce{
  DampingForce(): c_(0.0){}; // Default constructor
  DampingForce(double c): c_(c){};

  template <typename NODE>
  Point operator()(NODE n,double t){
    return (-c_*n.value().vel);
    (void) t;
  }
  
  double c_; // Damping constant
};

// Constraints
// Plane constraint
struct PlaneConstraint{
  PlaneConstraint(): z_(-1e5){}; // Default constructor 
  // z value set a large number that no point will ever violate

  PlaneConstraint(double z): z_(z){};
  template<typename Graph>
  void operator()(Graph& g,double t){
    for (auto n=g.node_begin(); n!=g.node_end();++n){
      if (dot((*n).position(),Point(0,0,1))<z_){
        (*n).position().z = z_;
        (*n).value().vel.z = 0;
      }
    }
    (void) t;
  }
  double z_; 
};

// Sphere constraint
struct SphereConstraint{
  SphereConstraint() :c_(Point(1e5,1e5,1e5)), radius_(0){}; // Default constructor
  SphereConstraint(Point c,double radius): c_(c), radius_(radius){};

  template<typename Graph>
  void operator()(Graph& g,double t){
    for (auto n=g.node_begin(); n!=g.node_end();++n){
      auto dist = norm((*n).position()-c_);
      auto Ri = ((*n).position()-c_)/dist;
      if (dist<radius_){
        (*n).position() = c_+radius_*Ri;
        (*n).value().vel -= dot((*n).value().vel,Ri)*Ri;
      }
    }
    (void) t;
  } 
  Point c_;
  double radius_;
};

// Sphere constraint with remove nodes
struct SphereConstraint2{
  SphereConstraint2() :c_(Point(1e5,1e5,1e5)), radius_(0){}; // Default constructor
  SphereConstraint2(Point c,double radius): c_(c), radius_(radius){};

  template<typename Graph>
  void operator()(Graph& g,double t){
    for (auto n=g.node_begin(); n!=g.node_end();++n){
      auto dist = norm((*n).position()-c_);
      if (dist<radius_){
        g.remove_node(n);
      }
    }
    (void) t;
  }
  Point c_;
  double radius_; 
};

// Constant position constraint
struct FixedPosition{
  template<typename Graph>
  void operator()(Graph& g, double t){
    for (auto ni= g.node_begin(); ni!=g.node_end(); ++ni){
      if (((*ni).position()== Point(0,0,0)) ||
          ((*ni).position()== Point(1,0,0))){
        (*ni).value().vel = Point(0,0,0);
      }
    }
    (void) t;
  }
};

// Make combined constraint
/* Given any two constraints (including FixedPoint), it uses the default constructor of the other constraints 
* defined and applies all the constraints to the graph
*
* Given any three (of the 4 constraints), it again applies the default constructor
* of the other constraint and applies all the constraints to the graph
* 
* @pre The constraints have to specified in the order defined in the template. 
* The order is PlaneConstraint, SphereConstraint, SphereConstraint2, FixedPoint
*/
template <typename PC, typename SC1, typename SC2, typename FP>
struct make_combined_constraint{
  make_combined_constraint(PC cons1, FP cons4)
  : cons1_(cons1), cons2_(SphereConstraint())
   ,cons3_(SphereConstraint2()), cons4_(cons4) {};

  make_combined_constraint(SC1 cons2, FP cons4)
  : cons1_(PlaneConstraint()), cons2_(cons2)
   ,cons3_(cons3_), cons4_(cons4) {};

  make_combined_constraint(SC2 cons3, FP cons4)
  : cons1_(PlaneConstraint()), cons2_(SphereConstraint())
   ,cons3_(cons3), cons4_(cons4) {};

  make_combined_constraint(PC cons1, SC1 cons2, FP cons4)
  : cons1_(cons1), cons2_(cons2), cons3_(SphereConstraint2()),cons4_(cons4) {};

  make_combined_constraint(PC cons1, SC2 cons3, FP cons4)
  : cons1_(cons1), cons2_(SphereConstraint()), cons3_(cons3),cons4_(cons4) {};  
  

  template<typename Graph>
  void operator()(Graph& g,double t){
    cons1_(g,t);
    cons2_(g,t);
    cons3_(g,t);
    cons4_(g,t);
  }

  PC cons1_;
  SC1 cons2_;
  SC2 cons3_;
  FP cons4_;
};


/* Make combined force struct
* Given two forces, gravity force and mass spring force it sets 
* the damping force to zero and combines the three forces
*
* Given three forces, gravity force, mass spring force, damping force, i
* it returns the combination of the three forces
*/
template <typename GF, typename MSF, typename DF>
struct make_combined_force{  
  make_combined_force(GF force1, MSF force2)
  : force1_(force1), force2_(force2), force3_(DampingForce()){}; 

  make_combined_force(GF force1, MSF force2, DF force3)
  : force1_(force1), force2_(force2), force3_(force3){};

  template <typename NODE>
  Point operator()(NODE n, double t){
    return (force1_(n,t)+force2_(n,t)+force3_(n,t));
  }
  
  GF force1_;
  MSF force2_;
  DF force3_;
};

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  //Problem1Force(double K, double L) : K_(K), L_(L) {};
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    if ( n.position () == Point(0 ,0 ,0) || n.position ()== Point(1 ,0 ,0))
      return Point(0 ,0 ,0);

    Point gravity(0,0,-n.value().mass*grav); // For gravity force
    Point spring(0,0,0); // For spring force

    // Calculating spring force
    for(auto iit=n.edge_begin();iit!=n.edge_end();++iit){
      auto K_ = (*iit).value().K_;
      auto L_ = (*iit).value().L_;
      double force_const = -(K_/(*iit).length())*((*iit).length()-L_);
      spring += force_const*((*iit).node1().position()-(*iit).node2().position());
    }
  
    return (spring+gravity);
    (void) t;
  }
  //double K_,L_; // Spring force constants
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
// #if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
// #endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  // Initial positions,velocity already fixed

  //Fix mass of each node
  for (auto n=graph.node_begin(); n!=graph.node_end();++n){
    (*n).value().mass =  1.0/(double) graph.num_nodes();
  }

  // Fix Spring constant and Length of spring
  //double K =100.0;
  // for (auto e=graph.edge_begin(); e!=graph.edge_end();++e){
  //   std::cout<< (*e).length()<< std::endl;
  // }
  //auto L = (*graph.edge_begin()).length();

  //Fix K, L for each edge
  for (auto ni=graph.node_begin();ni!=graph.node_end();++ni){
    for (auto iit=(*ni).edge_begin(); iit!= (*ni).edge_end();++iit){
      (*iit).value().K_ =  100;
      (*iit).value().L_ = (*iit).length();
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
      double t_end = 5;
      using mcf_type = make_combined_force<GravityForce,MassSpringForce
                                                        ,DampingForce> ;
      using mcc_type = make_combined_constraint<PlaneConstraint,SphereConstraint
                                                ,SphereConstraint2,FixedPosition>;
      
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        // symp_euler_step(graph, t, dt, Problem1Force());
        
        mcf_type mcf((GravityForce()),(MassSpringForce()),
                      (DampingForce(1.0/(double) graph.num_nodes()))); 

        mcc_type mcc((PlaneConstraint(-0.75)),
                      (SphereConstraint2(Point(0.5,0.5,-0.5),0.15)),
                       FixedPosition());

        symp_euler_step(graph, t, dt, mcf,mcc);

        // // Update viewer with nodes' new positions
        // viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        // viewer.set_label(t);

        // Clear the viewer ’s nodes and edges
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes ’ new positions and new edges
        viewer.add_nodes(graph.node_begin(),graph.node_end(),node_map);
        viewer.add_edges(graph.edge_begin(),graph.edge_end(),node_map);

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
