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
#include <math.h>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
//damping constant
double c;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};
/** Custom structure of data to store with Edges */
struct EdgeData{
  double L;	//Edge spring rest length
  double K;	//Edge Spring constant
  EdgeData() : L(0), K(0) {}
};
// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/*Plane constraint*/
struct plane_constraint{
	void operator()(GraphType& g, double t){
		(void) t;
		for(auto it=g.node_begin(); it != g.node_end(); ++it){
			if((*it).position().z < -0.75)	//if this node violates the constraint 
 			{
				(*it).position().z = -0.75;
				(*it).value().vel.z = 0;
			}
		}
	}
};
/*Sphere constraint*/
struct sphere_constraint{
	template<typename G>
	void operator()(G& g, double t){
		(void) t;
		Point c = Point(0.5,0.5,-0.5);
		double r = 0.15;
		for(auto it=g.node_begin(); it != g.node_end(); ++it){
			double dist = norm((*it).position()-c);
			if(dist < r)	//if this node violates the constraint 
 			{
				Point ri = ((*it).position()-c)/dist;
				(*it).position() = c+r*ri;
				(*it).value().vel = (*it).value().vel - ((*it).value().vel*ri)*ri;
			}
		}
	}
};
/*Sphere constraint to remove nodes*/
struct sphere_constraint2{
	template<typename G>
	void operator()(G& g, double t){
		(void) t;
		Point c = Point(0.5,0.5,-0.5);
		double r = 0.15;
		auto it=g.node_begin();
		while(it!=g.node_end()){
			double dist = norm((*it).position()-c);
			if(dist < r)	//if this node violates the constraint 
 			{
				it=g.remove_node(it);
			}
			else
				++it;
		}
	}
};
/* To combine two constraints */
template<typename Cons1, typename Cons2>
struct combined_constraint{
	Cons1 c1_;
	Cons2 c2_;
	combined_constraint(Cons1 c1=Cons1(), Cons2 c2=Cons2()):c1_(c1),c2_(c2){}
	void operator()(GraphType& g, double t){
		c1_(g,t);
		c2_(g,t);
	}
};
/* To combine two constraints, used as a helper function to combined_constraints*/
template<typename Cons1, typename Cons2>
combined_constraint<Cons1,Cons2> make_combined_constraint(Cons1 C1, Cons2 C2){
	return combined_constraint<Cons1, Cons2>(C1,C2);
}
/* To combine three constraints, used as a helper function to combined_constraints*/
template<typename Cons1, typename Cons2, typename Cons3>
combined_constraint<combined_constraint<Cons1,Cons2>, Cons3> make_combined_constraint(Cons1 C1, Cons2 C2, Cons3 C3){
	return combined_constraint<combined_constraint<Cons1, Cons2>, Cons3>(combined_constraint<Cons1,Cons2>(C1,C2),C3);
}

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
	
    //Check constraint and then update
    if(n.position()!=Point(0,0,0) && n.position()!=Point(1,0,0))
    {
	//apply constraints
	auto c = make_combined_constraint(sphere_constraint2(),plane_constraint());
	c(g,t);
	//update velocties
	n.value().vel += force(n, t) * (dt / n.value().mass);
    }
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
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    if(n.position()==Point(0,0,0) ||n.position()==Point(1,0,0))
	return Point(0,0,0);    
 
    double mi = n.value().mass;
    
    Point spring = Point(0,0,0);
    for (auto it = n.edge_begin(); it!=n.edge_end(); ++it)
    {
	Point p1 = n.position();
        Point p2 = (*it).node2().position();
	if(p1.x==p2.x && p1.y==p2.y && p1.z==p2.z)
		p2=(*it).node1().position();
	
        double L = (*it).value().L;
	double K = (*it).value().K;

 	Point xi_xj = p1-p2;
        
	double ed = norm(xi_xj);
        spring+=(-1.0)*K*(xi_xj)*(ed-L)/(double)ed;
    }
    
    Point f_grav = Point(0,0,-grav*mi);
    
    return spring+f_grav;
  }
};
struct GravityForce{
  /* Returns the force of gravity applied to @a n at time @a t*/
  template <typename NODE>
  Point operator()(NODE n, double t){
	(void) t;
       	double mi = n.value().mass;
	return Point(0,0,-grav*mi);
  }
};

struct MassSpringForce{
  /* Returns the spring forces applied to @a n at time @a t*/
  
  Point operator()(Node n, double t){
	(void) t;
 	Point spring = Point(0,0,0);
    	for (auto it = n.edge_begin(); it!=n.edge_end(); ++it)
    	{
		Point p1 = n.position();
        	Point p2 = (*it).node2().position();
		if(p1.x==p2.x && p1.y==p2.y && p1.z==p2.z)
			p2=(*it).node1().position();
	
	        double L = (*it).value().L;
		double K = (*it).value().K;
	
	 	Point xi_xj = p1-p2;
        
		double ed = norm(xi_xj);
	        spring+=(-1.0)*K*(xi_xj)*(ed-L)/(double)ed;
	}
	return spring;
    }
};
struct DampingForce{
	Point operator()(Node n, double t){
		(void) t;
		return (-1)*c*n.value().vel;
	}
};
/* Returns a Point, which is the sum of values returned from two functions */
template<typename Force1, typename Force2>
struct combined_force{
	Force1 f1_;
	Force2 f2_;
	combined_force(Force1 f1=Force1(), Force2 f2=Force2()):f1_(f1),f2_(f2){}
	Point operator()(Node n, double t){
		return f1_(n,t)+f2_(n,t);
	}
};
/* Combines two functions and returns the sum of their values*/
template<typename Force1, typename Force2>
combined_force<Force1,Force2> make_combined_force(Force1 F1, Force2 F2){
	return combined_force<Force1, Force2>(F1,F2);
}
/* Combines three functions and returns the sum of their values*/
template<typename Force1, typename Force2, typename Force3>
combined_force<combined_force<Force1,Force2>, Force3> make_combined_force(Force1 F1, Force2 F2, Force3 F3){
	return combined_force<combined_force<Force1, Force2>, Force3>(combined_force<Force1,Force2>(F1,F2),F3);
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

   //Initial conditions for node
  for(auto it = graph.node_begin(); it!=graph.node_end();++it)
  {
  	(*it).value().mass = (double)1/graph.num_nodes(); //mass
	(*it).value().vel = Point(0,0,0);		  //velocity
	
  }
  //Initial conditions for edge
  for(auto it=graph.edge_begin();it!=graph.edge_end();++it)
  {
	(*it).value().L = (*it).length();		//spring rest length
	(*it).value().K = 100;				//spring constant
  }
  c = (double)1/graph.num_nodes();	//damping constant

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
        auto f = make_combined_force(GravityForce(),MassSpringForce(), DampingForce());
	
	symp_euler_step(graph, t, dt, f);
        
	//Clear the viewer's nodes and edges
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
