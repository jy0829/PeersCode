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
double N;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  int pin; //Node fixed position (nonzero if pinned)
  NodeData() : vel(0), mass(1), pin(0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double L;       //edge length
  double K;     //edge spring constant
  EdgeData() : L(1.0), K(100.0) {}
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
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  //check constraints before velocity update
  constraint(g,t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

//FORCES

//Gravitational force acting on @a n at time @a t.
struct GravityForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    return (n.value().mass)*Point(0,0,-grav);
  }
};

//Spring force acting on @a n at time @a t.
struct SpringForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    //initialize force 
    Point F_s = Point(0,0,0);
    //iterate over all adjacent nodes
    for (auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      auto n2 = (*i).node2();
      Point n2_dist = n.position()-n2.position();
      //calculate and add spring force
      F_s =F_s - (*i).evalue().K*n2_dist*(norm(n2_dist)-(*i).evalue().L)
             /norm(n2_dist);  
    }
    return F_s;
  }
};

//Damping force acting on @a n at time @a t.
struct DampingForce {

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    return (-1.0/N)*(n.value().vel);
  }
};

//struct to combine forces
template<typename F1, typename F2>
struct CombinedForce2 {
  F1 f1;
  F2 f2;
  template <typename Node>
  Point operator()(Node n, double t) {
    return f1(n,t) + f2(n,t);

  }
};

//function for two types of forces
template<typename F1, typename F2>
CombinedForce2<F1, F2> make_combined_forces(F1 f1, F2 f2) {
  return {f1, f2};
}

//function for three types of forces
template<typename F1, typename F2, typename F3>
CombinedForce2<CombinedForce2<F1, F2>, F3> make_combined_forces(F1 f1, F2 f2, F3 f3) {
  return {{f1, f2}, f3};
}

//CONSTRAINTS

//Fixed poisiton constraint possibly acting on @a n at time @a t+dt.
struct FixedConstraint {
  
  template<typename G>
  void operator()(G& g, double t) {
     (void)t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i) {
      Node n = *i;
      //check if node was constrained upon initialization, make valid
      if (n.value().pin == -1) {
        n.position() = Point(0,0,0);
        n.value().vel = Point(0,0,0);
      }
      else if (n.value().pin == 1) {
        n.position() = Point(1,0,0);
        n.value().vel = Point(0,0,0);
      }
    }
  }
};

//Plane constraint possibly acting on @a n at time @a t+dt.
struct PlaneConstraint {
  //constraint for plane
  template <typename G>
  void operator()(G& g, double t) {
    (void)t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i) {
      Node n = *i;
      if (n.position().z < -0.75) {
        n.position().z = -0.75;
        n.value().vel = Point(0, 0, 0);
      }
    }
  }
};

//First sphere constraint possibly acting on @a n at time @a t+dt.
struct SphereConstraint1 {
  //sphere properties
  Point center = Point(0.5, 0.5, -0.5);
  double radius = 0.15;
  template <typename G>
  void operator()(G& g, double t) {
    (void)t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i) {
      Node n = *i;
      Point diff = n.position()-center;
      Point unit = diff/norm(diff);
      if (norm(diff) < radius) {
        n.position() = radius*unit + center;
        n.value().vel = n.value().vel-(n.value().vel*unit)*unit;
      }
    }
  }
};

//Second sphere constraint possibly acting on @a n at time @a t+dt.
struct SphereConstraint2 {
  //sphere properties
  Point center = Point(0.5, 0.5, -0.5);
  double radius = 0.15;
  template <typename G>
  void operator()(G& g, double t) {
    (void)t;
    auto i = g.node_begin();
    while (i != g.node_end()) {
      Node n = *i;
      Point diff = n.position()-center;
      if (norm(diff) < radius) {
        i = i.remove_node(i);
      }
      else {
        ++i;
      }
    }
  }
};


//struct to combine constraints
template<typename C1, typename C2>
struct CombinedConstraint2 {
  C1 c1;
  C2 c2;
  template <typename G>
  void operator()(G& g, double t) {
    c1(g,t);
    c2(g,t);

  }
};

//function for two types of constraints
template<typename C1, typename C2>
CombinedConstraint2<C1, C2> make_combined_constraints(C1 c1, C2 c2) {
  return {c1, c2};
}

//function for three types of constraints
template<typename C1, typename C2, typename C3>
CombinedConstraint2<CombinedConstraint2<C1, C2>, C3> 
make_combined_constraints(C1 c1, C2 c2, C3 c3) {
  return {{c1, c2}, c3};
}

//function for four types of constraints
template<typename C1, typename C2, typename C3, typename C4>
CombinedConstraint2<CombinedConstraint2<C1, C2>, CombinedConstraint2<C3,C4>> 
make_combined_constraints(C1 c1, C2 c2, C3 c3, C4 c4) {
  return {{c1, c2}, {c3, c4}};
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

  // INITIALIZATION
  //get total number of nodes
  N = graph.num_nodes();
  //iterate over nodes to add mass and check if pinned point
  for (auto i = graph.node_begin(); i != graph.node_end(); ++i) {
    //dereference pointer, add mass info
    (*i).value().mass=1.0/N;
    if ((*i).position() == Point(0,0,0)) {
      (*i).value().pin = -1;
    }
    else if ((*i).position() == Point(1,0,0)) {
      (*i).value().pin = 1;  
    }
    for (auto j = (*i).edge_begin(); j != (*i).edge_end(); ++j) {
      (*j).evalue().L = (*j).length();
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
      double dt = 0.001/2.0;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;

        auto all_forces = make_combined_forces(GravityForce(),SpringForce(), 
             DampingForce());

        auto all_constraints = make_combined_constraints(FixedConstraint(),
             SphereConstraint2(), PlaneConstraint(), SphereConstraint1());
        
        symp_euler_step(graph, t, dt, all_forces, all_constraints);

        //clear viewer and nodes since nodes/edges being updated
        viewer.clear();
        node_map.clear();

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
