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
  NodeData() {} //invalid constructor
  NodeData(const Point& vel_, const double& mass_): vel(vel_),mass(mass_){}

};

struct EdgeData{
  double K;
  double L;
  EdgeData() {}
  EdgeData(const double& K_ij_, const double& L_ij_): K(K_ij_),L(L_ij_){}
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
double symp_euler_step(G& g, double t, double dt, F force, C combined_constr) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  combined_constr(g,t);

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
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    //(void) n; (void) grav;    // silence compiler warnings
    if (n.position()==Point(0,0,0) || n.position() == Point(1,0,0)){return Point(0,0,0);}

    Point F_spring(0,0,-n.value().mass*grav);
    for(auto iter =n.edge_begin();iter!=n.edge_end();++iter){
     F_spring-=K*((*iter).node1().position()-(*iter).node2().position())*(1-L/(*iter).length());
    }
    Point F_grav = Point(0,0,-n.value().mass*grav);
    (void) t;
    Point F=F_grav+F_spring;
    return F;
  }
double K;
double L;


Problem1Force(const double& K_, const double& L_):K(K_),L(L_){}


};

/**GravityForce defines the force relative to gravity*/

struct GravityForce{
    template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;//not unused variable
    return Point(0,0,-n.value().mass*grav);
  }
};

/**MassSpringForce defines the force relative to the spring forces of every other springs*/

struct MassSpringForce{
template <typename NODE>
Point operator()(NODE& n, double t){
  (void) t;
  Point F_spr(0,0,0);
 for (auto iter =n.edge_begin();iter!=n.edge_end();++iter){
//    F_spr-= (*iter).value().K_ij*(1.-(*iter).value().L_ij/(*iter).length())*((*iter).node1().position()-(*iter).node2().position());
 		Edge edge = *iter;
 		Point substract = n.position()=edge.node2().position();
 		double mynorm = norm(substract);
 		F_spr+=-(edge.value().K)*substract*(mynorm-edge.value().L)/mynorm;
                  }
return F_spr;
   
  }
};

/**DampingForce defines the force relative to damping*/

struct DampingForce {
  double damp_const;
  DampingForce(){}
  DampingForce(const double& damp_const_): damp_const(damp_const_){}
  template <typename NODE>
  Point operator()(NODE n, double t){
  (void) t;
  return -n.value().vel*damp_const;
 }
  
};

/**Combine different number of forces of different nature*/
template <typename F1, typename F2>
struct CombinedForce{
  F1 f1;
  F2 f2;

  CombinedForce( F1 f1_,F2 f2_) : f1(f1_),f2(f2_){}


  /**summation of forces f1 and f2 applied at node @a n at time @a t
  *@param[n] node to consider
  *@param[t] time at which the force is applied
     */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n,t) + f2(n,t);
  }
};

/**combine 2 forces of the structure CombinedForce*/
template <typename F1, typename F2>
CombinedForce<F1,F2> make_combined_force(F1 f1, F2 f2){
  return CombinedForce<F1,F2>(f1,f2);
}

/**combine 3 forces of the structure CombinedForce*/

template <typename F1, typename F2, typename F3>
CombinedForce< CombinedForce<F1,F2>, F3 > make_combined_force(F1 f1, F2 f2, F3 f3){
  return CombinedForce< CombinedForce<F1,F2>,F3 > (make_combined_force(f1,f2),f3);
}

//CONSTRAINTS

/**Constraint at a specific node with velocity ==0*/

struct FixedNodeConstraint{
  Point point;
  FixedNodeConstraint(Point point_): point(point_){
  }
  template <typename NODE>
  void operator()(GraphType& g, double t){
    for (auto iter =g.node_begin(); iter!=g.node_end();++iter){
      if((*iter).position()==Point(0,0,0) || (*iter).position()==Point(1,0,0))//boundary conditions
      {
        (*iter).value().vel =Point(0,0,0);
      }
    }
    (void) t;
  }
};

struct PlaneConstraint{
  Point point;
  Point origin_z;
  PlaneConstraint(Point point_ = Point(0,0,1), Point origin_z_ =Point(0,0,-0.75)) : point(point_),
  origin_z(origin_z_){}
/**The nodes of the graph have to be above the plan z=-0.75
*@param[g] the graph
*@param[t] the time
*/
void operator()(GraphType& g, double t){
  (void) t;
  for (auto iter =g.node_begin();iter!= g.node_end();++iter){
    if ((*iter).position().z<-0.75){
      (*iter).position().z=-0.75;
      (*iter).value().vel.z=0;//use the .z from the Point object
    }
  }
}
};

/**constraint from the sphere*/
struct SphereConstraint{
/**The nodes can't stay inside a sphere S((0.5,0.5,-0.5),0.15)*/
  void operator()(GraphType& g, double t){
    (void) t;
    Point center(0.5,0.5,-0.5);
    double radius = 0.15;
    for (auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
        Point& point =(*iter).position();
        if(norm(point-center)<radius){
          Point normalize = (point-center)/norm(point-center);
          point = (normalize*radius)+center;
          Point updatedVel = dot(normalize,(*iter).value().vel)*normalize;
          (*iter).value().vel = (*iter).value().vel-updatedVel;
        }
    }
  }

};

/**RemoveSphereConstraint is embedding a method to
remove nodes that are on the surface of a sphere. Center and radius are fixed
No constructor is provided here
*/ 

struct RemoveSphereConstraint{

  void operator()(GraphType& g, double t){
    (void) t;
    double radius =0.15;
    Point center(0.5,0.5,-0.5);
    for (auto iter =g.node_begin();iter !=g.node_end();++iter){
      if (norm((*iter).position()-center)<radius){
        g.remove_node((*iter));
      }
    }
  }


};


/**A way to combine constraints*/
template<typename C1, typename C2>
struct CombineConstraint{
  C1 c1;
  C2 c2;

  void operator()(GraphType& g, double t){
    c1(g,t);
    c2(g,t);
  }
  CombineConstraint(){} //non valid constructor
  CombineConstraint(C1 c1_, C2 c2_) : c1(c1_),c2(c2_){}
};
/**function that encode 2 combined constraints
*@param[in] two constraints to be combined
*@returns CombineConstraint member
*/
template<typename C1,typename C2>
CombineConstraint<C1,C2> make_combined_constraint(C1 c1,C2 c2){
  return CombineConstraint<C1,C2>(c1,c2);
}
/**function that encode 3 combined constraints
*@param[in] 3 constraints to be combined
*@returns CombineConstraint member
*/
template<typename C1,typename C2,typename C3>
CombineConstraint<CombineConstraint<C1,C2>,C3> make_combined_constraint(C1 c1,C2 c2,C3 c3){
CombineConstraint<C1,C2> combined_first_two(c1,c2);//use constructor
return make_combined_constraint(combined_first_two,c3);
}

struct Corner_points{
  void operator()(GraphType& g, double t){
    (void) t;
    Point three_zero(0,0,0);
    Point two_zero_1_one(1,0,0);
    for(auto iter =g.node_begin();iter!=g.node_end();++iter){
      if((*iter).position()==three_zero || (*iter).position()==two_zero_1_one){
        (*iter).value().vel = three_zero;
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

  //double K=100;

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p,NodeData(Point(0,0,0),0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
#if 1
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  //initialization for problem 1
  for (auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
    (*iter).value().mass = 1. / graph.num_nodes();
    (*iter).value().vel = Point(0, 0, 0);
  } 
  //initialization for problem 2
   for(auto iter2 =graph.edge_begin();iter2!=graph.edge_end();++iter2){
     EdgeData& edgedata = (*iter2).value();
     edgedata.L=(*iter2).length();
     edgedata.K=100;
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

      //definition of the forces
      GravityForce grav_F;
      MassSpringForce mass_spring_F;
      DampingForce damping_F(1/graph.num_nodes());
      auto combined_F = make_combined_force(grav_F,mass_spring_F,damping_F);
      
      //definition of the constraints
      RemoveSphereConstraint rmvesphereconstrt;
      Corner_points cornerpoint;
      PlaneConstraint planeconstrt;
      auto combined_C =make_combined_constraint(cornerpoint, planeconstrt,rmvesphereconstrt);

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {



        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, combined_F, combined_C);
        viewer.clear();
        node_map.clear();
        // Update viewer with nodes ’ new positions and new edges
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
    
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
