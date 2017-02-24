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

//#define REFRESH

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {} // default constructor
};

struct EdgeData {
  double ori_len;
  double sp_const;
  EdgeData():ori_len(1), sp_const(1){}
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
    ///////////////////ignore nodes /////////////////////
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) continue;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  constraint(g);  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    //////////////////igone nodes////////////////////////
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) continue;
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
    Point result = Point(0,0,0);
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return result;
    }
    for (auto edge_ptr = n.edge_begin(); edge_ptr != n.edge_end(); ++ edge_ptr){
      // add the spring force
      auto e = *edge_ptr;
      NODE n2 = (*edge_ptr).node2();
      double distance = e.length();
      //assert(distance >= 0);
      Point direction = (n.position() - n2.position()) / distance;
      result += -e.value().sp_const * direction * (distance - e.value().ori_len);
    }
    result += n.value().mass * grav * Point(0,0,-1);
    (void) t;    // silence compiler warnings
    return result;
  }
};

struct GravityForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point result = Point(0,0,0);
    result += n.value().mass * grav * Point(0,0,-1);
    (void) t;;    // silence compiler warnings
    return result;
  }
};

struct  MassSpringForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point result = Point(0,0,0);
    for (auto edge_ptr = n.edge_begin(); edge_ptr != n.edge_end(); ++ edge_ptr){
      // add the spring force
      auto e = *edge_ptr;
      NODE n2 = (*edge_ptr).node2();
      double distance = e.length();
      
      //assert(distance >= 0);
      Point direction = (n.position() - n2.position()) / distance;
      result += -e.value().sp_const * direction * (distance - e.value().ori_len);
    }
    (void) t;    // silence compiler warnings
    return result;
  }
};


struct DampingForce{
  DampingForce(double num):c(num) {}
  double c = 0.0;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point result = Point(0,0,0);
    result += - c * n.value().vel;
    (void) t;  // silence compiler warnings
    return result;
  }
};



template <class T1, class T2>
struct TwoForces{
  T1 f1;
  T2 f2;
  template <typename NODE>
  Point operator()(NODE n, double t){
    return f1(n,t) + f2(n,t);
  }
};

template <class T1, class T2>
TwoForces<T1, T2> make_combined_force(T1 f1, T2 f2){
  (void) f1; (void) f2; // silence compiler warning;
  return TwoForces<T1, T2>();
}


template <class T1, class T2, class T3>
struct ThreeForces{
  ThreeForces(T1 t1, T2 t2, T3 t3):f1(t1), f2(t2), f3(t3) {}
  T1 f1;
  T2 f2;
  T3 f3;
  template <typename NODE>
  Point operator()(NODE n, double t){
    return f1(n,t) + f2(n,t) + f3(n,t);
  }
};

template <class T1, class T2, class T3>
ThreeForces<T1, T2, T3> make_combined_force(T1 f1, T2 f2, T3 f3){
  return ThreeForces<T1, T2, T3>(f1, f2, f3);
}


struct FloorConstraint{
  FloorConstraint(double z):z_(z){};
  double z_;
  void operator()(GraphType& g){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      Point& p = (*it).position();
      if (p[2] < z_){
        p[2] = z_;
        (*it).value().vel[2] = 0;
      }
    }
  }
};

struct SphereConstraint{
  SphereConstraint(Point p, double r): p_(p),r_(r){};
  Point p_;
  double r_;
  void operator()(GraphType& g){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      Point& p = (*it).position();
      if (norm(p - p_) < r_){
        Point direction = (p - p_) /(norm(p - p_));
        p = p_ + direction * r_;
        Point& velocity = (*it).value().vel;
        velocity = velocity - (velocity * direction) * direction;
      }
    }
  }
};


struct SphereRemove{
  SphereRemove(Point p, double r): p_(p),r_(r){};
  Point p_;
  double r_;
  void operator()(GraphType& g){
    for (auto it = g.node_begin(); it != g.node_end();){
      Point& p = (*it).position();
      if (norm(p - p_) < r_){
        it =  g.remove_node(it);
      } else{
        ++it;
      }
    }
  }
};

template <class F1, class F2, class F3>
struct GeneralizedConstraint{
  GeneralizedConstraint(F1 f1, F2 f2, F3 f3): f1_(f1), f2_(f2), f3_(f3) {}
  F1 f1_;
  F2 f2_;
  F3 f3_;
  void operator()(GraphType& g){
    f1_(g);
    f2_(g);
    f3_(g);
  }
};

template < class F1, class F2, class F3>
GeneralizedConstraint<F1, F2, F3> make_constraint(F1 f1, F2 f2, F3 f3){
  return GeneralizedConstraint<F1, F2, F3>(f1, f2, f3);
}

struct Dummy{
  Point operator()(GraphType& g){
    (void)g;
    return Point();
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
  double numnodes = graph.num_nodes();
  for (auto currnode : nodes){
    currnode.value().vel = Point();
    currnode.value().mass = 1/numnodes;
  }

  for (GraphType::EdgeIterator edgeitr = graph.edge_begin(); edgeitr != graph.edge_end(); ++ edgeitr){
    Edge e = *edgeitr;
    e.value().sp_const = 100;
    e.value().ori_len = e.length();
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
      auto totalForce = make_combined_force(GravityForce(), MassSpringForce(), DampingForce{double(1)/graph.num_nodes()});
      auto totalConstraint = make_constraint(FloorConstraint(-0.75), SphereConstraint(Point(0.5,0.5,-0.5), 0.15), Dummy());
      //auto totalConstraint = make_constraint(FloorConstraint(-0.75), SphereRemove(Point(0.5,0.5,-0.5), 0.15), Dummy());
      //auto totalConstraint = make_constraint(FloorConstraint(-0.75), Dummy(), Dummy());
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;

        symp_euler_step(graph, t, dt, totalForce, totalConstraint);

        // clear
#ifdef REFRESH        
        viewer.clear();
        node_map.clear();
#endif
        //
        //        std::cout << "num nodes and edges " << graph.num_nodes() << "  " << graph.num_edges() << std::endl;
        // Update viewer with nodes' new positions        
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
#ifdef REFRESH
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
#endif        
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
