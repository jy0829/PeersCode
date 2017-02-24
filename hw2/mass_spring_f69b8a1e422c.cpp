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
    double K;
    double L;
    EdgeData(double K, double L) : K(K), L(L) {}
    EdgeData() : K(100.0), L(1.0) {}
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
   /*
   if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)){
       n.value().vel += force(n, t) * (dt / n.value().mass);
   }
   */

   // With constraints:
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
    //double K = 100;
    //double L = (*(n.edge_begin())).length();

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        return Point(0,0,0);

    Point fspring = Point(0); 
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
        Point pos = n.position() - (*it).node2().position();
        fspring -=  (*it).value().K * pos/norm(pos) * (norm(pos) - (*it).value().L);
    }

    Point fgrav = n.value().mass * Point(0,0,-grav);

    return fspring + fgrav;
  }
};

struct GravityForce {
   template <typename NODE>
   Point operator()(NODE n, double t) {
     return n.value().mass * Point(0, 0, -grav);
   }
 };

struct MassSpringForce {
    template <typename NODE>
    Point operator()(NODE n, double t) {
      Point fspring = Point(0); 
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
          Point pos = n.position() - (*it).node2().position();
          fspring -=  (*it).value().K * pos/norm(pos) * (norm(pos) - (*it).value().L);
      }
      return fspring;
    }
};

struct DampingForce {
    double c; // damping constant

    DampingForce(double c) : c(c){}

    template <typename NODE>
    Point operator()(NODE n, double t) {
      return -c * n.value().vel;
    }
};

template <typename F1, typename F2>
struct PairForces {
  /** Functor that returns the sum of two functors,
   *  that each represent a force
   */
  F1 f1;
  F2 f2;
  
  PairForces(F1 f1, F2 f2) : f1(f1), f2(f2) {}
  template <typename NODE>
  Point operator()(NODE n, double t){
    return f1(n, t) + f2(n, t);
  }
};

template <typename F1, typename F2, typename F3>
struct TripletForces {
  /** Functor that returns the sum of three functors
   *  that each represent a force
   */
  F1 f1;
  F2 f2;
  F3 f3;
  
  TripletForces(F1 f1, F2 f2, F3 f3) : f1(f1), f2(f2), f3(f3) {}
  template <typename NODE>
  Point operator()(NODE n, double t){
    return f1(n, t) + f2(n, t) + f3(n, t);
  }
};

template <typename F1, typename F2>
PairForces<F1, F2> make_combined_force(F1 f1, F2 f2){
    return PairForces<F1, F2>(f1, f2);
}

template <typename F1, typename F2, typename F3>
TripletForces<F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3){
  return TripletForces<F1, F2, F3>(f1, f2, f3);
}

/* Constraint to keep nodes at a certain fixed position pt
 * There is only one node that possibly violates the constraint
 * at the start.
 * @param[in] pt : the position constrained
 * @a n is the possible node that violates the constraint
 */
struct FixedConstraint {
public:
    FixedConstraint(GraphType& g, Point pt) : graph(g), pt(pt) {
        for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
          if ((*it).position() == pt){
            n = (*it);
            violated = true;
            break;
          }
        }
    };

    void operator()(GraphType& graph, double t) {
        if(violated){
            n.position() = pt;
            n.value().vel = Point(0,0,0);
        }
    }

private:
    bool violated;
    Node n;
    Point pt;
    GraphType graph;

};

struct PlaneConstraint {
public:
  PlaneConstraint(double z): z(z){}
  void operator()(GraphType& g, double t){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      if ((*it).position().z < z){
        (*it).position().z = z;
        (*it).value().vel.z = 0;
      }
    }
  }

private:
  double z;
};

struct SphereConstraint {
public:
  SphereConstraint(Point c, double r): c(c), r(r){}
  void operator()(GraphType& g, double t){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      Point dist = (*it).position() - c;
      auto n = norm(dist);
      if (n < r){
          Point Ri = dist / n;
          (*it).position() = c + r * Ri;
          (*it).value().vel -= dot((*it).value().vel, Ri) * Ri;
      }
    }
  }

private:
  Point c;
  double r;
};

struct RemoveSphereConstraint {
public:
  RemoveSphereConstraint(Point c, double r): c(c), r(r){}
  void operator()(GraphType& g, double t){
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      Point dist = (*it).position() - c;
      auto n = norm(dist);
      if (n < r){
          g.remove_node(*it);
      }
    }
  }

private:
  Point c;
  double r;
};



template <typename C1, typename C2>
struct PairConstraints {
  /** Functor that returns the sum of two functors,
   *  that each represent a constrant
   */
  C1 c1;
  C2 c2;
  PairConstraints(C1 c1, C2 c2) : c1(c1), c2(c2) {}
  void operator()(GraphType& g, double t){
    c1(g, t);
    c2(g, t);
  }
};

template <typename C1, typename C2, typename C3>
struct TripletConstraints {
  /** Functor that returns the sum of two functors,
   *  that each represent a constrant
   */
  C1 c1;
  C2 c2;
  C3 c3;
  TripletConstraints(C1 c1, C2 c2, C3 c3) : c1(c1), c2(c2), c3(c3) {}
  void operator()(GraphType& g, double t){
    c1(g, t);
    c2(g, t);
    c3(g, t);
  }
};

template <typename C1, typename C2>
PairConstraints<C1, C2> make_combined_constraint(C1 c1, C2 c2){
    return PairConstraints<C1, C2>(c1, c2);
}

template <typename C1, typename C2, typename C3>
TripletConstraints<C1, C2, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3){
  return TripletConstraints<C1, C2, C3>(c1, c2, c3);
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for nodes and edges
  int N = graph.size();
  for(auto it = graph.node_begin(); it!=graph.node_end(); ++it){
      // initial velocity are already at 0
      (*it).value().mass = 1.0/N;
      (*it).value().vel = Point(0);
      for(auto it2 = (*it).edge_begin(); it2 != (*it).edge_end(); ++it2){
          (*it2).value() = EdgeData(100.0, (*it2).length());
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
      double t_end = 5.0;

      auto force = make_combined_force(GravityForce(), MassSpringForce(), DampingForce(1.0/N));
      auto constraintPoint1 = FixedConstraint(graph, Point(0,0,0));
      auto constraintPoint2 = FixedConstraint(graph, Point(1,0,0));
      auto planeConstraint = PlaneConstraint(-0.75);
      auto sphereConstraint = RemoveSphereConstraint(Point(0.5, 0.5, -0.5), 0.15);

      auto constraint = make_combined_constraint(constraintPoint1, constraintPoint2, planeConstraint);

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        // Check constraints at start as well to identify the FixedConstraints
        constraint(graph, t);
        sphereConstraint(graph, t);

        symp_euler_step(graph, t, dt, force);
        
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
