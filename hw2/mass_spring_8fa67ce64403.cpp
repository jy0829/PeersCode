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
  double len_;
  double const_;
  EdgeData() : len_(1), const_(100.0) {}
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
   
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply the constraints using functor  
  constraint(g);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if(n.position() != Point(0, 0, 0) && n.position() != Point(1, 0, 0))
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Return the length of an edge as the Euclidean distance between its 
*   incident nodes.
*/
template<>
double Edge::length() const {
  return norm(Node(g_, uid1_).position() - Node(g_, uid2_).position());
};

/** Object that models the gravitational force acting on a node */
struct GravityForce {
  template<typename NODE>
  Point operator()(NODE n) {
    return Point(0, 0, -grav*n.value().mass);
  }
};

/** Object that models the spring force acting on a node */
struct MassSpringForce {
  template<typename NODE> 
  Point operator()(NODE n) {
    Point force = Point(0);
    for(auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      auto e = *ei;
      Point neighbor_pos = e.node2().position();
      assert(e.value().const_ == 100.0);
      force += -e.value().const_ * (n.position() - neighbor_pos) * 
               (1 - e.value().len_ / norm(n.position()- neighbor_pos));
    }
    return force;
  }
};     

/** Object which models the drag force on a node as proportional to the 
*   node's velocity. The object stores the damping coefficient as a parameter
*/
struct DampingForce {
  double c;
  DampingForce(double proportionality_c) : c(proportionality_c) {} 
  template<typename NODE>
  Point operator()(NODE n) {
    return -c*n.value().vel;
  }
};

/** Object which combines foces of types F1 and F2 */
template<typename F1, typename F2>
struct MakeCombinedF {
  F1 f1_;
  F2 f2_;
  MakeCombinedF(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}

  template<typename NODE>
  Point operator()(NODE n) {
    return f1_(n) + f2_(n);
  }
};

/** Function that creates a functor equivalent to the sum of two inputted 
*   functors. 
*/
template<typename F1, typename F2>
MakeCombinedF<F1, F2> make_combined_force(F1 f1, F2 f2) {
  return MakeCombinedF<F1, F2>(f1, f2);
}

/** Function that creates a functor equivalent to the sum of the three
*   inputted functors.
*/
template<typename F1, typename F2, typename F3>
MakeCombinedF<MakeCombinedF<F1, F2>, F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return MakeCombinedF<MakeCombinedF<F1, F2>, F3>(make_combined_force(f1, f2), f3);
}

/** Constraint functor that searches for nodes violating the CHECK function
*   and applies the FIX function.
*/
template<typename CHECK, typename FIX>
struct Constraint {
  CHECK check;
  FIX fix;
  Constraint(CHECK ch, FIX fx) : check(ch), fix(fx) {}  

  template<typename GRAPH>
  void operator()(GRAPH& g) {
    for(auto nodeit = g.node_begin(); nodeit != g.node_end(); ++nodeit) 
       if(check((*nodeit))) fix((*nodeit));
  } 
};    

/** Functor which combines constraints */
template<typename C1, typename C2>
struct MakeCombinedC {
  C1 c1_;
  C2 c2_;  
  MakeCombinedC(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}

  template<typename GRAPH>
  void operator()(GRAPH& g) {
    c1_(g);
    c2_(g);
  }
};

/** Function which applies the two given constraints */
template<typename C1, typename C2>
MakeCombinedC<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  return MakeCombinedC<C1, C2>(c1, c2);
};

/** Function which applies the three given constraints */
template<typename C1, typename C2, typename C3>
MakeCombinedC<MakeCombinedC<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  return MakeCombinedC<MakeCombinedC<C1, C2>, C3>(make_combined_constraint(
            make_combined_constraint(c1, c2), c3));
};


/** Functor which applies the constant node constraint. Complexity is O(1),
 *  but an O(N) search is required for initialization. This implementation is
 *  faster than the more flexible implementation involving a functor for each
 *  constrained point, due to the O(N) initialization cost of each such functor.
 */
template<typename GRAPH>
struct ConstantConstr {
  Point::size_type node0;
  Point::size_type node1;
  //Search for node at the origin and node at Point(1, 0, 0) and store node
  //indices (perform search once when initialized).
  ConstantConstr(GRAPH& g) {
    for(auto nodeit = g.node_begin(); nodeit != g.node_end(); ++nodeit) {
      if((*nodeit).position() == Point(0, 0, 0)) node0 = (*nodeit).index();
      if((*nodeit).position() == Point(1, 0, 0)) node1 = (*nodeit).index();
    }
  }
  
  void operator()(GRAPH& g) {
    g.node(node0).value().vel = Point(0, 0, 0);
    g.node(node1).value().vel = Point(0, 0, 0);
  }
};
 
/** Functor which checks if given node violates the plane constraint. 
 *  The operator returns true if the z-coordinate of the node is below
 *  the z = -.75 plane
 */   
struct PlaneCheck {
  template<typename NODE>
  bool operator()(NODE n) {
    return dot(n.position(), Point(0, 0, 1)) < -.75;
  }
};

/** Functor which fixes the position of a node if it violates the plane 
 *  constraint. The position is set to the point on the plane z = -.75 closest
 *  to the node.
 */
struct PlaneFix {
  template<typename NODE>
  void operator()(NODE n) {
    n.position() = Point(n.position()[0], n.position()[1], -.75);
    n.value().vel -= Point(0, 0, n.value().vel[2]);
  }
};

/** Functor which checks if a given node lies inside the sphere centered at
 *  a given Point with radius specified by a given double.
 */ 
struct SphereCheck {
  Point c;
  double r;   
  SphereCheck(Point center, double radius) : c(center), r(radius) {}

  template<typename NODE>
  bool operator()(NODE n) {
    return norm(n.position() - c) < r;
  }
};

/** Functor which fixes a node relative to the sphere constraint. The functor
 *  sets the position of the node to the point on the sphere centered at c with *  radius r closest to the node's current position.
 */
struct SphereFix {
  Point c;
  double r;
  SphereFix(Point center, double radius) : c(center), r(radius) {}
  
  template<typename NODE>
  void operator()(NODE n) {
    Point unormal = (n.position() - c) / norm(n.position() - c);
    n.position() = r*unormal + c;
    n.value().vel = n.value().vel - dot(n.value().vel, unormal) * unormal;
  }
};

/** Functor which removes nodes */  
struct Rem {
  template<typename NODE>
  void operator()(NODE n) {
    n.g_->remove_node(n);
  }
};

/** Functor which checks if a node is isolated */
struct IsoCheck {
  template<typename NODE>
  bool operator()(NODE n) {
    return n.degree() == 0;
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
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Set initial conditions for your nodes, if necessary.

  // Initialize node mass for constant graph density. 
  double N = graph.num_nodes();
  for(auto it = graph.node_begin(); it != graph.node_end(); ++it)  
    (*it).value().mass = 1/N;
  // Initialize spring-rest length to initial edge length
  for(auto eit = graph.edge_begin(); eit != graph.edge_end(); ++eit) 
    (*eit).value().len_ = (*eit).length();
 
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
      double dt = 0.001 / 2.0;
      double t_start = 0;
      double t_end = 5.0;

      //Create constant, plane, and sphere constraints
      ConstantConstr<GraphType> c1 = ConstantConstr<GraphType>(graph);
      Constraint<SphereCheck, SphereFix> c2 = Constraint<SphereCheck, SphereFix>
          (SphereCheck(Point(.5,.5,-.5),.15), SphereFix(Point(.5,.5,-.5), .15));
      Constraint<SphereCheck, Rem> c22 = Constraint<SphereCheck, Rem>
          (SphereCheck(Point(.5,.5,-.5),.15), Rem());
      Constraint<IsoCheck, Rem> c23 = Constraint<IsoCheck, Rem>(IsoCheck(), Rem());
      Constraint<PlaneCheck, PlaneFix> c3 = Constraint<PlaneCheck, PlaneFix>
          (PlaneCheck(), PlaneFix());

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, //GravityForce(),
       make_combined_force(MassSpringForce(), GravityForce(), 
        DampingForce(1/N)), 
        make_combined_constraint(c1, c22, c23));
    
        //Clear the viewer's nodes and edges
        viewer.clear();
        node_map.clear();
    
        // Update viewer with nodes' new positions and new edges
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
    

        // Update viewer with nodes' new positions
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
