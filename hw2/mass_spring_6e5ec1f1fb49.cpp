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
  double c;        //< Node damping coefficient
  NodeData() : vel(0), mass(1), c(1) {}
};

/** Custom structure of data to store with Nodes */
struct EdgeData {
  double K;       //< Edge elastic constant
  double L;       //< Edge initial length
  EdgeData() : K(1), L(1) {}
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
//double symp_euler_step(G& g, double t, double dt, F force) {
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position

  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    if (n.position() != Point(0,0,0) and n.position() != Point(1,0,0))
      n.position() += n.value().vel * dt;
  }

  //Apply constraints before computing the force and velocity
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position() != Point(0,0,0) and n.position() != Point(1,0,0))
      n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  return t + dt;
}

struct NullForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t, (void) n;    // silence compiler warnings
    return Point(0,0,0);
  }
};

template <typename f1, typename f2, typename f3>
struct combined_force{
  f1 force1_;
  f2 force2_;
  f3 force3_;
  combined_force(f1 force1,f2 force2, f3 force3):force1_(force1), force2_(force2), force3_(force3){}
  template <typename NODE>
  Point operator()(NODE n, double t){
    return force1_(n,t)+force2_(n,t)+force3_(n,t);
  }
};

template <typename f1, typename f2>
combined_force<f1,f2,NullForce> make_combined_force(f1 force1, f2 force2){
  return combined_force<f1,f2,NullForce>(force1,force2,NullForce());
}

template <typename f1, typename f2, typename f3>
combined_force<f1,f2,f3> make_combined_force(f1 force1, f2 force2, f3 force3){
  return combined_force<f1,f2,f3>(force1,force2,force3);
}


struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return Point(0,0,-grav*n.value().mass);
  }
};

struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    double c = n.value().mass;
    (void) t;
    return -c*n.value().vel;
  }
};

struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point Fi = Point(0,0,0);
    for (auto ei = n.edge_begin(); ei!= n.edge_end(); ++ei){
      Edge e = *ei;
      Node n2 = e.node2();
      double L = e.length();
      double L_rest = e.value().L;
      Fi -= e.value().K*(n.position()-n2.position())*(L-L_rest)/L;
    }
    (void) t;
    return Fi;
  }
};

struct PlaneConstraint {
  double z_ = -0.75;
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto ni = g.node_begin(); ni!= g.node_end(); ++ni){
      Node n = *ni;
      if (inner_prod(n.position(),Point(0,0,1))<z_){
        n.position() = n.position()*Point(1,1,0)+Point(0,0,-0.75);
        n.value().vel = n.value().vel*Point(1,1,0);
      }
    }
    (void) t;
  }
};

struct SphereConstraint {
  double r_ = 0.15;
  Point c_ = Point(0.5,0.5,-0.5);
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto ni = g.node_begin(); ni!= g.node_end(); ++ni){
      Node n = *ni;
      if (norm_2(n.position()-c_)<r_){
        Point Ri = (n.position()-c_)/norm_2(n.position()-c_);
        n.position() = c_+Ri*r_;
        n.value().vel -= Ri*inner_prod(Ri,n.value().vel);
      }
    }
    (void) t;
  }
};

struct SphereConstraint_remove {
  double r_ = 0.15;
  Point c_ = Point(0.5,0.5,-0.5);
  template <typename GRAPH>
  void operator()(GRAPH& g, double t) {
    for (auto ni = g.node_begin(); ni!= g.node_end(); ++ni){
      Node n = *ni;
      if (norm_2(n.position()-c_)<r_){
        unsigned deleted_node = g.remove_node(n);
        (void) deleted_node;
      }
    }
    (void) t;
  }
};

template <typename c1, typename c2>
struct constraints{
  c1 cons1_;
  c2 cons2_;
  constraints(c1 cons1,c2 cons2):cons1_(cons1), cons2_(cons2){}
  template <typename GRAPH>
  void operator()(GRAPH& g, double t){
    cons1_(g,t);
    cons2_(g,t);
    (void) t;
  }
};

template <typename c1, typename c2>
constraints<c1,c2> combined_constraints(c1 cons1, c2 cons2){
  return constraints<c1,c2>(cons1,cons2);
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
  // Set initial conditions for your nodes, if necessary.
  double m = 1/double(graph.num_nodes());
  for (auto ni = graph.node_begin(); ni!= graph.node_end(); ++ni){
      Node n = *ni;
      n.value().mass = m;
      n.value().c = m;
  }
  // Set initial rest length and K of edges
  for (auto ni = graph.node_begin(); ni!= graph.node_end(); ++ni){
      Node n = *ni;
      for (auto ei = n.edge_begin(); ei!= n.edge_end(); ++ei){
        Edge e = *ei;
        e.value().L = e.length();
        e.value().K = 100;
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

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(),MassSpringForce(),DampingForce()),
                        combined_constraints(PlaneConstraint(),SphereConstraint_remove()));

        // Clear the viewer's nodes and edge
        viewer.clear();
        node_map.clear();
        // Update viewer with nodes' new positions and new edges
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
