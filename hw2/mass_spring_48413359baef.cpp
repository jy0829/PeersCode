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
  //NodeData(Point vel, double mass) : vel(vel), mass(mass) {}
};

struct EdgeData {
  double K;
  double L;
  EdgeData() : K(0), L(0) {}
  EdgeData(double K,double L) : K(K), L(L) {}
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
 * @tparam G::node_value_type supports a NodeData struct 
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

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraint(g,t);

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
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
       return Point(0,0,0);
    }
    else{
      Point gravity= Point(0,0,-n.value().mass*grav);
      Point spring = Point(0,0,0);
      for (auto it=n.edge_begin(); it!=n.edge_end(); ++it){
        Point node_pos=(*it).node2().position();
        spring-=K*(n.position()-node_pos)*(norm(n.position()-node_pos)-L)/norm(n.position()-node_pos);
      }
      return (gravity+spring);
    }
   (void) t;
  }
  Problem1Force(double K, double L) : K(K), L(L) {}
  double K;
  double L;
};

struct Problem2Force{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    else{
      Point gravity= Point(0,0,-n.value().mass*grav);
      Point spring = Point(0,0,0);
      for (auto it=n.edge_begin(); it!=n.edge_end(); ++it){
        double K=(*it).value().K;
        double L=(*it).value().L;
        Point node_pos=(*it).node2().position();
        spring-=K*(n.position()-node_pos)*(norm(n.position()-node_pos)-L)/norm(n.position()-node_pos);
      }
      return (gravity+spring);
    }
   (void) t;
  }
};

struct GravityForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return Point(0,0,-n.value().mass*grav);
  }
};

struct MassSpringForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point spring = Point(0,0,0);
    for (auto it=n.edge_begin(); it!=n.edge_end(); ++it){
      double K=(*it).value().K;
      double L=(*it).value().L;
      Point node_pos=(*it).node2().position();
      spring-=K*(n.position()-node_pos)*(norm(n.position()-node_pos)-L)/norm(n.position()-node_pos);
    }
    (void) t;
    return spring;
  }
};

struct DampingForce{
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -1/(double)num_nodes_*n.value().vel;
  }
  DampingForce(int num_nodes) : num_nodes_(num_nodes) {}
  int num_nodes_;
};


template<typename F1, typename F2>
struct Combine2{
  F1 f1;
  F2 f2;
  template<typename NODE>
  Point operator()(NODE n, double t){
    return f1(n,t)+f2(n,t);
  }
  Combine2(F1 f1, F2 f2) : f1(f1), f2(f2) {}
};

template<typename F1, typename F2>
Combine2<F1,F2> make_combined_force(F1 f1, F2 f2){
  return Combine2<F1,F2>(f1,f2);
}

template<typename F1, typename F2, typename F3>
Combine2<Combine2<F1,F2>,F3> make_combined_force(F1 f1, F2 f2, F3 f3){
  return Combine2<Combine2<F1,F2>,F3>(Combine2<F1,F2>(f1,f2),f3);
}


struct ConstConstraint{
  template<typename G>
  void operator()(G& g, double t){
    (void) t;
    auto it1 = std::next(g.node_begin(),g.u2i(id1_));
    auto it2 = std::next(g.node_begin(),g.u2i(id2_));
    (*it1).position()=Point(0,0,0);
    (*it2).position()=Point(1,0,0);
  }
  ConstConstraint(unsigned id1, unsigned id2): id1_(id1), id2_(id2) {}
  int id1_,id2_;
};

struct PlaneConstraint{
  template<typename G>
  void operator()(G& g, double t){
    (void) t;
    for(auto it=g.node_begin(); it!=g.node_end(); ++it){
      if ((*it).position().z<-0.75){
        (*it).position().z=-0.75;
        (*it).value().vel.z=0;
      }
    }
  }
};

struct SphereConstraint{
  template<typename G>
  void operator()(G& g,double t){
    (void) t;
    Point center=Point(0.5,0.5,-0.5);
    double radius=0.15;
    for(auto it=g.node_begin(); it!=g.node_end(); ++it){
      double n=norm((*it).position()-center);
      if(n<radius){
        Point R=((*it).position()-center)/n;
        (*it).position()=R*radius+center;
        (*it).value().vel-=inner_prod((*it).value().vel,R)*R;
      }
    }
  }
};

struct RemoveConstraint{
  template<typename G>
  void operator()(G& g,double t){
    (void) t;
    Point center=Point(0.5,0.5,-0.5);
    double radius=0.15;
    for(auto it=g.node_begin(); it!=g.node_end(); ++it){
      double n=norm((*it).position()-center);
      if(n<radius){
        g.remove_node(*it);
      }
    }
  }
};


template<typename C1, typename C2>
struct Combine2C{
  C1 c1;
  C2 c2;
  template<typename G>
  void operator()(G& g, double t){
    c1(g,t);
    c2(g,t);
  }
  Combine2C(C1 c1, C2 c2) : c1(c1), c2(c2) {}
};

template<typename C1, typename C2>
Combine2C<C1,C2> make_combined_constraint(C1 c1, C2 c2){
  return Combine2C<C1,C2>(c1,c2);
}

template<typename C1, typename C2, typename C3>
Combine2C<Combine2C<C1,C2>,C3> make_combined_constraint(C1 c1, C2 c2, C3 c3){
  return Combine2C<Combine2C<C1,C2>,C3>(Combine2C<C1,C2>(c1,c2),c3);
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
 
  
  //Finding the uid's of the nodes at the points to be held stationary and updating the mass of nodes
  int uid1,uid2;
  for (auto it=graph.node_begin(); it != graph.node_end(); ++it){
    auto node=*it;
    node.value().mass=1/(double)graph.num_nodes();
    if (node.position() == Point(0,0,0)){
      uid1=node.uid();
    }
    if (node.position() == Point(1,0,0)){
      uid2=node.uid();
    }
    //Setting the value for each edge
    for(auto it_edge=node.edge_begin(); it_edge !=node.edge_end(); ++it_edge){
      (*it_edge).value().L=(*it_edge).length();
      (*it_edge).value().K=100;
    }
  }



  //L for Problem 1
  //double L=(*graph.edge_begin()).length();

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

        //Problem1- Make sure L in uncommented above as well to run this 
        //symp_euler_step(graph, t, dt, Problem1Force(100,L),ConstConstraint(uid1,uid2));
       
        //Problem2  
        //symp_euler_step(graph, t, dt, Problem2Force(),ConstConstraint(uid1,uid2));
       
        //Problem3
        //symp_euler_step(graph, t, dt, make_combined_force<GravityForce,MassSpringForce,DampingForce>
        //(GravityForce(), MassSpringForce(),DampingForce(graph.num_nodes())),ConstConstraint(uid1,uid2));
       
        //Problem4
        symp_euler_step(graph, t, dt, make_combined_force<GravityForce,MassSpringForce,DampingForce>
        (GravityForce(),MassSpringForce(),DampingForce(graph.num_nodes())),
        make_combined_constraint<PlaneConstraint,RemoveConstraint,ConstConstraint>
        (PlaneConstraint(),RemoveConstraint(),ConstConstraint(uid1,uid2)));

        viewer.clear();
        node_map.clear();

        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

        // Update viewer with nodes' new positions
        //viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
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
