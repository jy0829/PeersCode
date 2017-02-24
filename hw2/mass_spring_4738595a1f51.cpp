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
// Spring Constant K
static constexpr double K = 100;

static double c;
//static double L;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,double>;
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
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    // HW2 #1: YOUR CODE HERE
    if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) {
	n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}


/** @struct PlaneConstraint
 * @brief a constraint that acts as a plane that forces nodes and edges to lay around it
 * 			instead of being allowed to enter inside of it.
 * @param Graph @a n, double @a t - input graph and the current time
 * @return void
 */
struct PlaneConstraint {
  /** Return the force applying to @a n at time @a t. */
  template <typename G>
  void operator()(G& g, double t) {
   
    (void)t; // silence the compiler about @a t

    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    	auto n = *it;
    	if (inner_prod(n.position(), Point(0,0,1)) < -0.75) {
		n.position().z = -0.75;
		n.value().vel.z = 0;
	}
    }
  }
};

/** @struct Sphere1Constraint
 * @brief a constraint that acts as a sphere and forces nodes and edges to lay around it
 * 			instead of being allowed to enter inside of it.
 * @param Graph @a n, double @a t - input graph and the current time
 * @return void
 */
struct Sphere1Constraint {
  /** Return the force applying to @a n at time @a t. */
  template <typename G>
  void operator()(G& g, double t) {
    (void)t; // silence the compiler about @a t

    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    	auto n = *it;
	
    	if (norm(n.position() - c) < r) {
		n.position() = ((n.position() - c)/norm(n.position() - c))*r + c;
		Point Ri = ((n.position() - c)/norm(n.position() - c));
		n.value().vel = n.value().vel-(inner_prod(n.value().vel, Ri)*Ri);
	}
    }
  }
};

/** @struct Sphere2Constraint
 * @brief a constraint that acts as a sphere that removes the nodes and edges that hit it 
 * @param Graph @a n, double @a t - input graph and the current time
 * @return void
 */
struct Sphere2Constraint {
  /** Return the force applying to @a n at time @a t. */
  template <typename G>
  void operator()(G& g, double t) {
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    	auto n = *it;
	
	double distFromCenter = norm(n.position() - c);
    	if (distFromCenter < r) {
		g.remove_node(n);
	}
    }
  }
};



/** @struct nullConstraint
 * @brief a constraint that merely takes up an argument space in the
 * 		make_combined_constraints function, serving as an insignificant "constraint".
 * 		It does not do anything.
 * @return void
 */
struct nullConstraint {
  /** Return the force applying to @a n at time @a t. */
  template <typename G>
  void operator()(G& g, double t) {
        (void)g; (void)t; // silence the compiler about @a g and @a t
	return;
  }
};

/** @struct totalC
 * @brief calls all of the individual constraint functions
 * @param c1, c2, c3 - the input constraints
 * @param Graph @a g, double @a t - input graph and the current time
 *
 * tparam T type of the constraints
 *
 * @return void
 */
template <typename C1, typename C2, typename C3>
struct totalC {
	C1 c1_; C2 c2_; C3 c3_;
	totalC(C1 c1, C2 c2, C3 c3) : c1_(c1), c2_(c2), c3_(c3){}

	template <typename G>
	void operator()(G& g, double t) {
		c1_(g, t);
		c2_(g, t);
		c3_(g, t);
		return;
	}
};

/** @brief Returns a combination of all the constraints in a form that is usable by symp_euler_step
 * @param[in] c1, c2, c3 the constraints to combine (c2 and/or c3 are null if not specified)
 *
 * tparam T type of the constraints
 *
 * @return a functor whos () operator returns the combination of all the inputted constraints
 */
template <typename C1, typename C2 = nullConstraint, typename C3 = nullConstraint>
totalC<C1, C2, C3> make_combined_constraints(C1 c1, C2 c2 = nullConstraint(), C3 c3 = nullConstraint()) {
	return totalC<C1, C2, C3>(c1, c2, c3);
}


/** @struct GravityForce
 * @brief a force that acts as a gravity force on the nodes and edges
 * @param Node @a n, double @a t - input node and the current time
 * @return a Point object mimicking a gravity behavior
 */
struct GravityForce {
  /** Return the force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t; // silence the compiler about @a t
    return Point(0,0,-grav*n.value().mass);

  }
};

/** @struct MassSpringForce
 * @brief a force that acts as a mass spring force on the nodes and edges
 * @param Node @a n, double @a t - input node and the current time
 * @return a Point object mimicking a mass spring behavior
 */
struct MassSpringForce {
  /** Return the force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {

    Point forceS = Point(0,0,0);

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
	if (t == 0) {
		// set the initial values for all the edges
		(*it).value() = (*it).length();
	}
	Point p1 = n.position();
	Point p2 = (*it).node2().position();
	double thisL = (*it).value();

	forceS += ((-K)*(p1-p2)/norm(p1-p2))*(norm(p1-p2)-thisL);
    }
    return forceS;

  }
};

/** @struct DampingForce
 * @brief a force that acts as a damping force on the nodes and edges
 * @param Node @a n, double @a t - input node and the current time
 * @return a Point object mimicking a damping behavior
 */
struct DampingForce {
  /** Return the force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t; // silence the compiler about @a t
    return -c*n.value().vel;

  }
};

/** @struct nullForce
 * @brief a force that merely takes up an argument space in the make_combined_forces function,
 *        serving as an insignificant "force"
 * @return a Point object with all entries 0
 */
struct nullForce {
  /** Return the force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)n; (void)t; // silence the compiler about @a n and @a t
    return Point(0,0,0);

  }
};
/** @struct totalF
 * @brief returns a point which is the sum of all of the inputted forces
 * @param F1, F2, F3 - the input forces
 * @param Node @a n, double @a t - input node and the current time
 *
 * tparam T type of the forces
 *
 * @return object of type Point that is the sum of all the inputted forces
 */
template <typename F1, typename F2, typename F3>
struct totalF {
	F1 f1_; F2 f2_; F3 f3_;
	totalF(F1 f1, F2 f2, F3 f3) : f1_(f1), f2_(f2), f3_(f3){}

	template <typename NODE>
	Point operator()(NODE n, double t) {
		return f1_(n, t) + f2_(n, t) + f3_(n, t);
	}
};

/** @brief Returns a sum of all the inputted forces in a form that is usable by symp_euler_step
 * @param[in] f1, f2, f3 the forces to combine (f3 is null if not specified)
 *
 * tparam T type of the forces
 *
 * @return a functor whos () operator returns the sum of all the inputted forces
 */
template <typename F1, typename F2, typename F3 = nullForce>
totalF<F1, F2, F3> make_combined_force(F1 f1, F2 f2, F3 f3 = nullForce()) {
	return totalF<F1, F2, F3>(f1, f2, f3);
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

  NodeData nd;
  nd.vel = Point(0,0,0);
  nd.mass = (double)1/graph.num_nodes();

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
	(*it).value() = nd;
  }

  // the damping constant
  c = (double)1/graph.num_nodes();

  //L = norm(graph.node(0).position()-graph.node(1).position());

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

	GravityForce gf;
	MassSpringForce msf;
	DampingForce df;

	PlaneConstraint pc;
	Sphere1Constraint s1c;
	Sphere2Constraint s2c;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
	
	//symp_euler_step(graph, t, dt, make_combined_force(gf, msf), s1c);
        symp_euler_step(graph, t, dt, make_combined_force(gf, msf, df), make_combined_constraints(pc, s2c));

	// Clear the viewer's nodes and edges
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
