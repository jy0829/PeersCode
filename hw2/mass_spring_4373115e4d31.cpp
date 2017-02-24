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

// Indicies of nodes at positions Point(0,0,0) and Point(0,0,1)
int p0_idx;
int p1_idx;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       //< Edge spring constant
  double L;       //< Edge rest length
  EdgeData() : K(100), L(0) {}
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

    if (n.position() == Point(0,0,0)) {
      p0_idx = n.index();
    }
    if (n.position() == Point(1,0,0)) {
      p1_idx = n.index();
    }
    
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

struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point gravity = Point(0,0,-grav);
    (void) t;
    return n.value().mass * gravity;  
  }
};

struct MassSpringForce {
  /** Return the spring force applying to @a n at time @a t. */ 
  template <typename NODE>
  Point operator()(NODE n, double t) {
    Point f = Point(0,0,0);
    Point p1 = n.position();
    
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      auto e = (*ei);
      Point p2 = e.node2().position();
      double dist = norm(p1-p2);
      f += - e.value().K * ((p1-p2)/dist) * (dist - e.value().L);
    } 
    (void) t;
    return f;
  }
};

struct DampingForce {
  /** Return the damping force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -c_ * n.value().vel;  
  }
  double c_;
  DampingForce(double c) : c_(c) {}
};


template<typename F1, typename F2>
struct CombineForce {
  /** Combine two forces. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n,t) + f2_(n,t);
  }
  CombineForce(F1 f1, F2 f2) : f1_(f1), f2_(f2) {}
  F1 f1_;
  F2 f2_;
};

// Combines two forces using ()operator in struct CombineForce
template <typename F1, typename F2>
CombineForce<F1,F2> make_combined_force(F1 f1, F2 f2) {
  return CombineForce<F1,F2>(f1,f2); 
}

// Combines three forces using ()operator in struct CombineForce
template <typename F1, typename F2, typename F3>
CombineForce<CombineForce<F1,F2>,F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return CombineForce<CombineForce<F1,F2>,F3>(CombineForce<F1,F2>(f1,f2),f3); 
}

struct ConstNodeConstraint {
  template <typename G>
  /** Updates the position of two globally defined nodes */
  void operator()(G& g, double t){
    g.node(p0_idx).position() = Point(0,0,0);
    g.node(p1_idx).position() = Point(1,0,0);
    (void) t;
  }
};

struct PlaneConstraint {
  template <typename G>
  /** Updates the position and velocity of nodes beneath the z = -0.75 plane. */
  void operator()(G& g, double t) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (n.position().z < -0.75) {
        n.position().z = -0.75;
        n.value().vel.z = 0;
      }
    }
    (void) t;
  }
};

struct SphereConstraint {
  template <typename G>
  /** Updates the position and velocity of nodes with the defined sphere. */
  void operator()(G& g, double t) {
    Point c = Point(0.5, 0.5, -0.5);
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point direction = n.position() - c;
      double dist = norm(direction);
      if (dist < 0.15) {
        Point R = direction/dist;
        n.value().vel -= dot(n.value().vel, R) * R;
        n.position() = n.position() + (0.15-dist) * R;
      }
    }
    (void) t;
  }
};

struct RemoveSphereConstraint {
  template <typename G>
  /** Removes nodes that come into contact with the defined sphere. */
  void operator()(G& g, double t) {
    Point c = Point(0.5, 0.5, -0.5);
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point direction = n.position() - c;
      double dist = norm(direction);
      if (dist < 0.15) {
	g.remove_node(n);
      }
    }
    (void) t;
  }
};
 
template<typename C1, typename C2>
struct CombineConstraint {
  /** Applies two constraints. */
  template <typename G>
  void operator()(G& g, double t) {
    c1_(g,t);
    c2_(g,t);
  }
  CombineConstraint(C1 c1, C2 c2) : c1_(c1), c2_(c2) {}
  C1 c1_;
  C2 c2_;
};

// Combines two constraints using the ()operator in struct CombineConstraint. 
template <typename C1, typename C2>
CombineConstraint<C1,C2> make_combined_constraint(C1 c1, C2 c2) {
  return CombineConstraint<C1,C2>(c1,c2); 
}

// Combines three constraints using the ()operator in struct CombineConstraint. 
template <typename C1, typename C2, typename C3>
CombineConstraint<CombineConstraint<C1,C2>,C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  return CombineConstraint<CombineConstraint<C1,C2>,C3>(CombineConstraint<C1,C2>(c1,c2),c3); 
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

  // Sets all the nodes mass initial value
  double mass = (float)1/graph.num_nodes();
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    auto n = (*ni);
    n.value().mass = mass;
  }

  // Sets all the edges length initial value
  for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
    auto n = (*ni);
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      auto e = (*ei);
      e.value().L = norm(e.node1().position() - e.node2().position());
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
      double dt = 0.001/2; // Divided by 2 so grid3 works. Can go back to 0.001.
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt
			, make_combined_force(GravityForce(),MassSpringForce(), DampingForce(mass))
			, make_combined_constraint(ConstNodeConstraint(), PlaneConstraint(), RemoveSphereConstraint()));

	// Clear the viewer's nodes and edges
	viewer.clear();
	node_map.clear();

	// Update viewer with nodes' new positions and new edges
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
