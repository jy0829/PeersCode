
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
  double K;       //< spring constant
  double L;     //< length
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/* I believe this is the most efficient implementation of the constant constraint
 * there is no reason to design it taking graph in unless we modify the way the
 * nodes are read at the beginning (looking for those 2 specific points)
 * or using global variables. I implemented all 3 options and decided this was the best
 *
 */
struct const_pos_update{
  double dt_;
  const_pos_update(double dt) : dt_(dt){}
  void operator()(Node n){
    if (n.position() == Point(0,0,0) or n.position() == Point(1,0,0)){
      n.value().vel = Point(0,0,0);
    }
    n.position() += n.value().vel * dt_;
  }
};

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

  std::for_each(g.node_begin(), g.node_end(), const_pos_update(dt));
  /*
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    
  }
*/
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position

  std::for_each(g.node_begin(), g.node_end(), const_pos_update(dt));
  /*
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    
  }
  */

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    
  }

  constraint(g);
  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
double K, L;

Problem1Force(const double K, const double L) : K(K), L(L) {};

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t; (void) grav;    // silence compiler warnings

    if(n.position() == Point(0, 0, 0) || n.position() == Point(1,0,0)){
      return Point(0, 0, 0);
    }

    Point xi = n.position();
    Point F_spring = Point(0,0,0);
    /*
    std::for_each(n.edge_begin(), n.edge_end(), [&](Edge& e)
       {
       Point xj = e.node2().position();
       F_spring += -K * (xi-xj) * (norm(xi-xj) -L)/(norm(xi-xj));
       }
      );
    */
    L = *(n.edge_begin()).length();
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i){
      Point xj = (*i).node2().position();
      F_spring +=  -K * (xi-xj) * (norm(xi-xj) -L)/(norm(xi-xj));
    }
    return (n.value().mass * Point(0,0,-grav)) + F_spring;
  }
};

struct Problem2Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
double K, L;

//Problem1Force(const double K, const double L) : K(K), L(L) {};

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) n; (void) t; (void) grav;    // silence compiler warnings

    if(n.position() == Point(0, 0, 0) || n.position() == Point(1,0,0)){
      return Point(0, 0, 0);
    }
    double K, L, edgelen;
    Point xi = n.position();
    Point F_spring = Point(0,0,0);
    /*
    std::for_each(n.edge_begin(), n.edge_end(), [&](Edge& e)
       {
       Point xj = e.node2().position();
       F_spring += -K * (xi-xj) * (norm(xi-xj) -L)/(norm(xi-xj));
       }
      );
    */
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i){
      K = (*i).value().K;
      L = (*i).value().L;
      edgelen = (*i).length();

      Point xj = (*i).node2().position();
      F_spring +=  -K * (xi-xj) * (edgelen -L)/(edgelen);
    }
    return (n.value().mass * Point(0,0,-grav)) + F_spring;
  }
};

/*

template <typename T>
inline T make_combined_force(T n) {
  return n;
}

// Recursive case
template <typename T, typename... Args>
inline T make_combined_force(T n, Args... args) {
  return n + make_combined_force(args...);
}
*/

/*

template <typename...> struct SumTs;
template <typename T1> struct SumTs<T1> { typedef T1 type; };
template <typename T1, typename... Ts>
struct SumTs<T1, Ts...>
{
  typedef typename SumTs<Ts...>::type rhs_t;
  typedef decltype(std::declval<T1>() + std::declval<rhs_t>()) type;
};

//now the sum function
template <typename T>
T sum(const T& v) {
  return v;
}

template <typename T1, typename... Ts>
auto sum(const T1& v1, const Ts&... rest) 
  -> typename SumTs<T1,Ts...>::type //instead of the decltype
{
  return v1 + sum(rest... );
}
*/

/* Dear TA:
 * Pleae excuse all the commented code that I have in my submission.
 * I tried more than one approach for each thing and even though many of them
 * work not all of them satisfied HW requirements.
 *
 * I really tried to implement the force and constraint combiners using variadic
 * functions and it took me a long time and I never really got them to work
 * in this proccess a created a nice function that takes in any number of
 * objects that have the + operator defined and sums them.
 *
 * I have tried really hard to do this assignments to perfection but this one
 * wa really hard to do. I had almost finished all the HW when the new/unoffical
 * sort of required remove specs were given in class. I didn't have enough
 * time to reimplement and I am still unsure of the benefit of doing that
 * unless nodes were always accesible using uid's which is not the case
 * of the design proposed in class.
 *
 *
 */
template<typename F1, typename F2>
struct Combine2{
  F1 f1;
  F2 f2;
  template<typename NODE>
  Point operator()(NODE n, double t){
    return f1(n, t) + f2(n, t);
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


/* This is very similar to the previous force cobminer but combiens constraints
 * we were fored to do this in the HW specs but I believe this is a really bad
 * idea. we are performing O(n) searcher per each constraint.
 * I think it would be better to use a different approach and enforece all
 * constraint over a single search.
 *
 */
template<typename F1, typename F2>
struct Combine2c{
  F1 f1;
  F2 f2;
  template<typename GRAPH>
  void operator()(GRAPH& g){
    f1(g);
    f2(g);
  }
  Combine2c(F1 f1, F2 f2) : f1(f1), f2(f2) {}
};
template<typename F1, typename F2>
Combine2c<F1,F2> combine_constraints(F1 f1, F2 f2){
  return Combine2c<F1,F2>(f1,f2);
}
template<typename F1, typename F2, typename F3>
Combine2<Combine2<F1,F2>,F3> combine_constraints(F1 f1, F2 f2, F3 f3){
  return Combine2c<Combine2c<F1,F2>,F3>(Combine2c<F1,F2>(f1,f2),f3);
}


//UGLY BUT WORKS:

/*
template <typename force1, typename force2, typename force3>
struct make_combined_force{
  force1 f1_;
  force2 f2_;
  force3 f3_;
  bool flag = false;
  make_combined_force(force1 f1, force2 f2, force3 f3) :
    f1_(f1), f2_(f2), f3_(f3){};

  make_combined_force(force1 f1, force2 f2) :
    f1_(f1), f2_(f2), flag(True) {};

  template <typename NODE>
  Point operator()(NODE n, double t){
    if (flag)
      return f1_(n,t) + f2_(n,t);
    return f1_(n,t) + f2_(n,t) + f3_(n,t);
  }
};
*/

struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;    // silence compiler warnings
    /*if(n.position() == Point(0, 0, 0) || n.position() == Point(1,0,0)){
      return Point(0, 0, 0);
    }*/
    return (n.value().mass * Point(0,0,-grav));
  }
};


struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    /*if(n.position() == Point(0, 0, 0) || n.position() == Point(1,0,0)){
      return Point(0, 0, 0);
    }*/
    double K, L, edgelen;
    Point xi = n.position();
    Point F_spring = Point(0,0,0);
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i){
      K = (*i).value().K;
      L = (*i).value().L;
      edgelen = (*i).length();
      Point xj = (*i).node2().position();
      F_spring +=  -K * (xi-xj) * (edgelen -L)/(edgelen);
    }
    return F_spring;
  }
};

struct DampingForce {
  double c;
  DampingForce(const double c) : c(c) {};
  template<typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().vel * -c;
  }
};

struct PlaneConstraint {
  double z=-.75;

  struct planector{
    double val;
    planector(double val) : val(val) {}
    void operator()(Node n){
      if(n.position().z < val){
        n.position().z = val;
        n.value().vel.z = 0;
      }
    }
  };

  template<typename GRAPH>
  void operator()(GRAPH& g){
    std::for_each(g.node_begin(), g.node_end(), planector(z));
  }
};

struct SphereConstraint {
  double r = .15;
  Point c = Point(.5, .5, -.5);


  struct spheretor{
    double r;
    Point c;
    spheretor(double r, Point c) : r(r), c(c) {}
    void operator()(Node n){
      auto x = n.position();
      auto R = (x - c)/norm(x - c);
      if(norm(x-c) < r){
        n.position() = c + r*R;
        n.value().vel -= (n.value().vel * R)*R;
      }
    }
  };

  template<typename GRAPH>
  void operator()(GRAPH& g){
    std::for_each(g.node_begin(), g.node_end(), spheretor(r, c));
  }
};

struct SphereKiller {
  double r = .15;
  Point c = Point(.5, .5, -.5);

  template<typename GRAPH>
  void operator()(GRAPH& g){
    auto start = g.node_begin();
    auto end = g.node_end();

    while(start != end){
      if( norm((*start).position()-c) < r){
        g.remove_node(*start);
      } else {
        ++start;
      }
    }
  }
};

/*
struct ConstantConstraint{
  Node node1;
  Node node2;
  bool found;

  ConstantConstraint() : found(false) {}
  void operator()(GRAPH& g){
    if (!found){
      node1 = std::find_if(g.node_begin, g.node_end)
      node2 = std::find_if(g.node_)
    }
  }

};
*/

int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;
 // size_type constant_aux1, constant_aux2;

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

//mass

struct mass_vel_init {
  GraphType g;

  mass_vel_init(GraphType gra) : g(gra){}
  void operator()(Node n){
    n.value().vel = Point(0,0,0);
    n.value().mass = 1.0/g.num_nodes();
  }
};

mass_vel_init mvi(graph);
std::for_each(graph.node_begin(), graph.node_end(), mvi);

struct l_k_init {
  //GraphType g;
  //l_k_init(GraphType gra) : g(gra){}
  //Edge symm;
  void operator()(Edge e){
    /*
    e.value().L = e.length();
    e.value().K = 100;
    e.symmetric_edge().value().L = e.symmetric_edge().length();
    e.symmetric_edge().value().K = 100;
    */
    EdgeData auxdata;
    auxdata.K = 100;
    auxdata.L = e.length();
    e.initialize(auxdata);
  }
};

l_k_init lki;
std::for_each(graph.edge_begin(), graph.edge_end(), lki);




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
      //double K = 100;
      //double L = (*(graph.edge_begin())).length();

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        /* SORRY for all this mess it is not very clear how this HW has to be
         * submitted, I show a lot of diffrerent combinations.
         * NOTE: the render will be sort of flashing and that is due to
         * removing and readding all nodes and edges at each go which is 
         * required for the last problem.
         */ 
        //symp_euler_step(graph, t, dt, Problem1Force(K,L));
        //symp_euler_step(graph, t, dt, Problem2Force());
        //symp_euler_step(graph, t, dt, sum(GravityForce(), MassSpringForce()));
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes())));
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes())), combine_constraints(SphereConstraint(), PlaneConstraint()));
        // Update viewer with nodes' new positions
        
        //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes())), combine_constraints(SphereKiller(), PlaneConstraint()));
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0/graph.num_nodes())), SphereKiller());

        viewer.clear();
        node_map.clear();
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
