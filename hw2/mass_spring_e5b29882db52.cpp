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

/** Custom structure of data to store with Edge */
struct EdgeData {
  double K;       //< Spring Constant
  double L;     //< Initial edge length
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
    //Fixed nodes
    if (n.position() == Point(1, 0, 0))
	return Point(0,0,0);

    if (n.position() == Point(0, 0, 0))
	return Point(0,0,0); 
    
    //Spring Force
    Point force_spring = Point(0,0,0);
    for (auto Iter = n.edge_begin(); Iter != n.edge.end(); ++Iter) {
	Point force_diff = n.position() - (*Iter).node2.position();
	force_spring = force_spring-((*Iter).value().K)*force_diff *(norm(force_diff)-(*Iter).value().L) / norm*(force_diff);
	
    }
    
    //Gravity Force
    auto force_grav = Point(0, 0, -n.value().mass*grav);

    //Combined Force
    auto force_comb = force_spring + force_grav;

    (void) t;
    return force_comb;
  }
};

struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    
    //Gravity Force
    auto force_grav = Point(0, 0, -n.value().mass*grav);
    (void) t;
    return force_grav;
  }
};


struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    
    //Spring Force
    Point force_spring = Point(0,0,0);
    for (auto Iter = n.edge_begin(); Iter != n.edge_end(); ++Iter) {
	Point force_diff = n.position() - (*Iter).node2().position();
	force_spring = force_spring-((*Iter).value().K)*force_diff *(norm(force_diff)-(*Iter).value().L) / norm(force_diff);
	
    }
    (void) t;
    return force_spring;
  }
};

struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    
    //Damping Force
    auto force_damp = -n.value().vel*cons_damp;
    (void) t;
    return force_damp;
  }

  //Damping constant
  double cons_damp;

  //Constructor
  DampingForce() {cons_damp = 0;}
  DampingForce(double new_cons_damp) {cons_damp = new_cons_damp;}
};

template <typename F1, typename F2>
struct CombinedForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    
    //Combined Force
    auto force_comb = f1_(n,t) + f2_(n,t);
    
    return force_comb;
  }

  //Forces
  F1 f1_;
  F2 f2_;

  //Constructor
  CombinedForce() {}
  CombinedForce(F1 new_f1, F2 new_f2) {
  	f1_ = new_f1;
  	f2_ = new_f2;
  }
};

template <typename F1, typename F2>
CombinedForce<F1, F2> make_combined_force(F1 f1, F2 f2) {
	auto force_comb = CombinedForce<F1, F2>(f1,f2);
	return force_comb;
}

template <typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1, F2>,F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
	auto force_comb = CombinedForce<F1, F2>(f1,f2);
	auto force_comb_all = make_combined_force(force_comb,f3);
	return force_comb_all;
}

//Constraints

//PlanetConstraint
struct PlaneConstraint {
	void operator() (GraphType& graph, double t) {
		for (auto Iter = graph.node_begin(); Iter!=graph.node_end();++Iter) {
			if ((*Iter).position().z >= -0.75) {
				(*Iter).position().z = 0;				
			} else {
				(*Iter).position().z = -0.75;			
			}
		}	
	}
};

//PlanetConstraint
struct FixedPoints {
	void operator() (GraphType& graph, double t) {
		for (auto Iter = graph.node_begin(); Iter!=graph.node_end();++Iter) {
			//Fixed nodes
   			 if ((*Iter).position() == Point(1, 0, 0))
				(*Iter).position() = Point(0,0,0);

    			if ((*Iter).position() == Point(0, 0, 0))
				(*Iter).position() = Point(0,0,0);
		}	
	}
};

//Spere Constraint
struct SphereConstraint {
 
    void operator()(GraphType& graph, double t) {
    
    Point c(0.5, 0.5, -0.5);
    double r = 0.15;
    
    for (auto Iter = graph.node_begin(); Iter != graph.node_end(); ++Iter) {
      
      if (norm((*Iter).position() - c) < r) {
        Point R = ((*Iter).position() - c) / norm((*Iter).position() - c);
        
        (*Iter).position() = (R * r) + c;
        
        (*Iter).value().vel -= dot((*Iter).value().vel, R) * R;
        
      }
	}
  }
};

struct RemoveSphere {
  void operator()(GraphType& graph, double t) {
    Point c(0.5, 0.5, -0.5);
    double r = 0.15;

    for (auto Iter = graph.node_begin(); Iter != graph.node_end(); ++Iter) {
      if (norm((*Iter).position() - c) < r) {
        graph.remove_node((*Iter));
      }
    }
  }
};

//Combined Constraint
template <typename C1, typename C2>
struct CombinedConstraint {
  
  /** Calls both constraints applied to @a g at time @a t*/
  void operator()(GraphType& graph, double t) {
    c1_(graph, t);
    c2_(graph, t);
  }
  
  //Constraints
  C1 c1_;
  C2 c2_;
  //Constructor
  CombinedConstraint() {}
  CombinedConstraint(C1 new_c1, C2 new_c2) {
    c1_ = new_c1;
    c2_ = new_c2;
  }
};


template <typename C1, typename C2>
CombinedConstraint<C1, C2> make_combined_constraint(C1 c1, C2 c2) {
  return CombinedConstraint<C1, C2>(c1, c2);
}


template <typename C1, typename C2, typename C3>
CombinedConstraint<CombinedConstraint<C1, C2>, C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  
  CombinedConstraint<C1, C2> cons_comb(c1, c2); 
  auto cons_comb_all = make_combined_constraint(cons_comb, c3);
  return cons_comb_all;
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
  for(auto Iter = graph.node_begin(); Iter != graph.node_end(); ++Iter){
    (*Iter).value().vel = Point(0,0,0);
    (*Iter).value().mass = 1. / graph.num_nodes();
  }

    for(auto Iter = graph.edge_begin(); Iter != graph.edge_end(); ++Iter){
    (*Iter).value().K = 100;
    (*Iter).value().L = (*Iter).length();
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
        //For Problem 1 
        //symp_euler_step(graph, t, dt, Problem1Force());

        //For later problems
        //Forces
        DampingForce force_damp(1 / graph.num_nodes());
		GravityForce force_grav;
        MassSpringForce force_ms;
        
        auto force_all = make_combined_force(force_damp, force_grav, force_ms);
    
    	//Constraints
    	FixedPoints cons_fp;
    	PlaneConstraint cons_plane;
    	SphereConstraint cons_sphere;
    	RemoveSphere cons_rmsphere;
    	auto cons_all = make_combined_constraint(cons_fp, cons_plane, cons_rmsphere);


		//Euler step
        symp_euler_step(graph, t, dt, force_all, cons_all);

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
