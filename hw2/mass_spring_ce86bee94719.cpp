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
  double SpringConst;       //< Edge Spring constant 
  double Length;     //< Edge Resting Length 
  EdgeData() : SpringConst(100), Length(0.25) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using NodeIter  = typename GraphType::node_iterator;


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
template <typename G, typename F,typename C>
double symp_euler_step(G& g, double t, double dt, F force,C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  //Impose our constraints on the system
  constraint(g,t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
   
  return t + dt;
}


/* Gravity Force
*/
struct GravityForce{
   template <typename NODE>
   Point operator()(NODE n, double t) {
      (void) t; // silence compiler warning
      return Point(0.0,0.0,-1.0*n.value().mass*grav); 
   }
};

/* Damping Force
*/
struct DampingForce{
   template <typename NODE>
   Point operator()(NODE n, double t) {
      (void) t; // silence compiler warning
      return -1.0*n.value().mass*n.value().vel; 
   }
};

/* Mass Spring Force
*/
struct MassSpringForce{
   template <typename NODE>
   //@pre:templated node has a position() function
   Point operator()(NODE n, double t) {
    Point Force(0.0,0.0,0.0);
    Point Distance(0.0,0.0,0.0);
    double K;
    double L;
    for(auto it=n.edge_begin();it!=n.edge_end();++it){
       K=(*it).value().SpringConst;
       L=(*it).value().Length;
       Distance=(n.position()-(*it).node2().position()); 
       Force+=(-K*(norm(Distance)-L)/(norm(Distance)))*(Distance); 
    } 
      (void) t; // silence compiler warning
      return Force; 
   }
};
/* "Problem 1" Force
*@ Usage: Problem1Force(SpringConstant,Length)  
*/

struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  double K=20;
  double L=0.25;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position()==Point(0,0,0) || n.position()==Point(1,0,0)){
       (void) t; //silence compiler warning
       return Point(0,0,0);
    }
    //Find mass spring force then add gravity
    Point Force(0,0,0);
    // norm(Point a-Point b) should be euclidean distance
    for(auto it=n.edge_begin();it!=n.edge_end();++it){
       Point Distance;
       K=(*it).value().SpringConst;
       L=(*it).value().Length;
       Distance=(n.position()-(*it).node2().position()); 
       Force+=(-K*(norm(Distance)-L)/(norm(Distance)))*(Distance); 
    } 
    //Gravity = num_nodes*Point(0,0,-g)
    Force+=n.value().mass*Point(0,0,-grav);
    (void) t; //silence compiler warning
    return Force;
  }
};

//Constraint or force that does nothing
struct Empty{
   template <typename N>
   Point operator()(N& n,double t){
     (void) n; (void) t; // silence compiler warning
   return Point(0,0,0);
   }
   Empty(){}
};

/*Combined Force class
* @pre:Forces must have operator() that takes a templated node and time
* @return: operator() takes a templated node and time and returns a force as a Point
*/
template<typename f1,typename f2,typename f3=void>
class combined_force{
public: 
  combined_force(f1 first,f2 second,f3 third):f1_(first),f2_(second),f3_(third){} 
  template <typename NODE>
  Point operator()(NODE n, double t) {
   return f1_(n,t)+f2_(n,t)+f3_(n,t);
  }
private:  
  f1 f1_;
  f2 f2_;
  f3 f3_;
};

//Helper function to create a combined force with two forces
template<typename f1,typename f2>
combined_force<f1,f2,Empty> make_combined_force( f1 first,f2 second){
   return combined_force<f1,f2,Empty>(first,second,Empty());
}

//Helper function to create a combined force with three forces
template<typename f1,typename f2,typename f3>
combined_force<f1,f2,f3> make_combined_force( f1 first,f2 second, f3 third){
   return combined_force<f1,f2,f3>(first,second,third);
}

//Stationary point constraint that is constructed with two nodes 
//holds two nodes in the same position
struct StationaryPoints{
   Node n1_;
   Node n2_;
   Point firstPt_;
   Point secondPt_; 
   template <typename G>
   void operator()(G& g,double t){
     //cassert all graphs same
     n1_.position()=firstPt_; 
     n2_.position()=secondPt_; 
     n1_.value().vel=Point(0.0,0.0,0.0); 
     n2_.value().vel=Point(0.0,0.0,0.0); 
     (void) g; (void) t; // silence compiler warning
   }
   StationaryPoints(Node first,Node second):n1_(first),n2_(second),firstPt_(first.position()),secondPt_(second.position()){}
};

//Nearest Node helper function to initialize StationaryPoints
NodeIter nearest_node(const GraphType& g, const Point& point)
{
   /** Declaring comparator between two nodes*/ 
  //use lambda function
  auto dist=[&](const Node& n1,const Node & n2){return norm(n1.position()-point)<norm(n2.position()-point);};
  return std::min_element(g.node_begin(),g.node_end(),dist);
}


//Implements a planar boundary at z=-0.75
//Not required for hw but could implement for arbitray planar boundary
//with minor modifications
struct PlanarBoundary{
   Point plane_=Point(0.0,0.0,1.0);
   double val_=-0.75;
   template <typename G>
   void operator()(G& g,double t){
      for(auto it=g.node_begin();it!=g.node_end();++it){
         auto n=*it;
         if(n.position().z<val_){
           n.position().z=val_; 
           n.value().vel.z=0.0;
         }
      }
      (void) t; // silence compiler warning
      return;
   }
   PlanarBoundary(Point face,double position):plane_(face),val_(position){}
   PlanarBoundary(){}
};

// Spherical Constraint 
struct SphereConstraint{
   Point center_=Point(0.5,0.5,-0.5);
   double radius_=0.15;
   template <typename G>
   void operator()(G& g,double t){
      for(auto it=g.node_begin();it!=g.node_end();++it){
         auto n=*it;
         Point difference=n.position()-center_;
         if(norm(difference)<radius_){
           n.position()+=(radius_-norm(difference))*(difference/norm(difference)); 
           n.value().vel-=(dot(n.value().vel,difference/norm(difference)))*(difference/norm(difference));
         }
      }
      (void) t; // silence compiler warning
   return;
   }
};

// Spherical Constraint that removes nodes 
struct SphereHole{
   Point center_=Point(0.5,0.5,-0.5);
   double radius_=0.15;
   template <typename G>
   void operator()(G& g,double t){
      for(auto it=g.node_begin();it!=g.node_end();++it){
         auto n=*it;
         Point difference=n.position()-center_;
         if(norm(difference)<radius_){
            g.remove_node(n);
         }
      }
      (void) t; // silence compiler warning
   return;
   }
};

/*
template<typename f1,typename f2>
class combined_constraint{
public: 
  combined_constraint(f1 first,f2 second):f1_(first),f2_(second){} 
  template <typename G>
  void operator()(G& g, double t) {
   f1_(g,t);
   f2_(g,t);
   return ;
  }
private:  
  f1 f1_;
  f2 f2_;
};
*/

//A class that combines the effect of 3 other constraints
template<typename f1,typename f2,typename f3=void>
class combined_constraint{
public: 
  combined_constraint(f1 first,f2 second,f3 third):f1_(first),f2_(second),f3_(third){} 
  template <typename G>
  void operator()(G& g, double t) {
   f1_(g,t);
   f2_(g,t);
   f3_(g,t);
   return ;
  }
private:  
  f1 f1_;
  f2 f2_;
  f3 f3_;
};

//Helper function to create a combined constraint from 2 constraints
template<typename f1,typename f2>
combined_constraint<f1,f2,Empty> make_combined_constraint( f1 first,f2 second){
   return combined_constraint<f1,f2,Empty>(first,second,Empty());
}

//Helper function to create a combined constraint from 3 constraints
template<typename f1,typename f2,typename f3>
combined_constraint<f1,f2,f3> make_combined_constraint( f1 first,f2 second,f3 third){
   return combined_constraint<f1,f2,f3>(first,second,third);
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
  /* Sets mass of all nodes because velocity is initialized to 0 */
  for(auto it=graph.node_begin();it!=graph.node_end();++it){
     (*it).value().mass=(1.0/(double)graph.num_nodes());
  }
  //Determine stationary pts 
  Node Pt1=*(nearest_node(graph,Point(0,0,0)));
  Node Pt2=*(nearest_node(graph,Point(1,0,0)));
  // Set initial conditions for your graph, if necessary.
  /*Set initial edge length and spring constant for all edges */
  for(auto it=graph.edge_begin();it!=graph.edge_end();++it){
     ((*it).value()).Length=norm((*it).node1().position()-(*it).node2().position());
     ((*it).value()).SpringConst=100;
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
      //Create combined constraint and forces 
      auto InForce=make_combined_force(GravityForce(),MassSpringForce(),DampingForce());
      auto InConstraint=make_combined_constraint(StationaryPoints(Pt1,Pt2),PlanarBoundary(),SphereHole());
      // Begin the mass-spring simulation
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;
      //double Length=norm(graph.edge(1).node1().position()-graph.edge(1).node2().position());
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        symp_euler_step(graph, t, dt, InForce,InConstraint);
        //Clear Viewer's Nodes and Edges
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
