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


static constexpr double grav = 9.81; // Gravity in meters/sec^2
static double damping_constant; // Damping constant
unsigned fixed_node_1; // Fixed node constraint #1
unsigned fixed_node_2; // Fixed node constrain #2

/** Sphere struct used in constraint */
struct Sphere {
  Point center;
  float radius;
  Sphere(Point p, float r) : center(p),
                             radius(r) {}
                             
} sphere = {Point(0.5,0.5,-0.5),0.15}; // initialize sphere

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double spring_constant;
  double rest_length;
  EdgeData() : spring_constant(0), rest_length(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using NodeIter = typename GraphType::node_iterator;
using EdgeIter = typename GraphType::edge_iterator;
using size_type = unsigned;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type suppors NodeData type
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a constraints functor called as @a constraints(&g, t). Void return. 
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if(n.position() == Point(0,0,0)) fixed_node_1 = std::distance(g.node_begin(), it); 
    if(n.position() == Point(1,0,0)) fixed_node_2 = std::distance(g.node_begin(), it); 

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  
  // Apply constraints
  constraints(&g, t);
  
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Gravity force functor */
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {    
    (void)t;
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** MassSpring functor */
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {   
    (void)t;
    Point force_spring(0,0,0);
    Point delta_x(0,0,0);
    for(auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      auto e_info = e.value();
      Node adj = e.node2();
      delta_x = n.position() - adj.position();
      force_spring += -e_info.spring_constant * (delta_x / norm_2(delta_x)) * (norm_2(delta_x) - e_info.rest_length);
    }
    return force_spring;
  }
};

/** Damping force functor */
struct DampingForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    return -damping_constant * n.value().vel;
  }
};

/** Functor that combines two (force) types at a time */
template<class f1, class f2>
struct combine {
  /** Return the force applying to @a n at time @a t.
   * */
  f1 f1_;
  f2 f2_;
  
  combine(f1 x, f2 y) : f1_(x), 
                        f2_(y) {
  }

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    return f1_(n,t) + f2_(n,t); // Combine forces
  }
};

/** Function that combines all forces */
template<class f1, class f2, class f3>
combine<combine<f1,f2>,f3> make_combined_force(f1 x, f2 y, f3 z) {
  return  combine<combine<f1,f2>,f3>(combine<f1,f2>(x,y),z);
}

/** Fixed point constraints */
struct constant_constraint {
  void operator()(GraphType *g, double t) {
    (void)t;
    g->node(fixed_node_1).position() = Point(0,0,0);
    g->node(fixed_node_1).value().vel = Point(0,0,0);
    g->node(fixed_node_2).position() = Point(1,0,0);
    g->node(fixed_node_2).value().vel = Point(0,0,0);
  }
};

/** Plane constraint that prevents animation from falling past Z-plane = -0.75 */
struct quarter_constraint {
  void operator()(GraphType *g, double t) {
   (void)t;
   for(auto it = g->node_begin(); it != g->node_end(); ++it) {
      Point *curr = &(*it).position();
      Point *vel = &(*it).value().vel;
      if(curr->z < -0.75) {
        curr->z = -0.75;
        vel->z = 0;
      }
    } 
   }
};

/** Sphere constraint. If nodes from graph violating this constraint */
struct sphere_constraint {
  void operator()(GraphType *g, double t) {
    (void)t;
    for(auto it = g->node_begin(); it != g->node_end(); ++it) {
      Point *curr = &(*it).position();
      //Point *vel = &(*it).value().vel;
      float dist = norm(*curr - sphere.center); 
      if(dist < sphere.radius) {
        //Point surface_normal = (*curr - sphere.center) / dist;
        //*curr = *curr + surface_normal * (sphere.radius - dist); 
        //*vel = *vel - dot(*vel, surface_normal)*surface_normal;  
        g->remove_node(it);
      }    
    }
   }
};

/** Struct used to combine constraints */
template<class c1, class c2, class c3> 
struct constraints {
  c1 const_1;
  c2 const_2;
  c3 const_3;
  constraints(c1 x, c2 y, c3 z) : const_1(x),
                                  const_2(y),
                                  const_3(z) {
  }
  void operator()(GraphType *g, double t) {
    const_1(g, t);
    const_2(g, t);
    const_3(g, t);
  }
};

/** Function that combines constraints. Passed as a functor */
template<class c1, class c2, class c3>
constraints<c1,c2,c3> make_combined_constraints(c1 constraint_1, c2 constraint_2, c3 constraint_3) {
  return  constraints<c1,c2,c3>(constraint_1, constraint_2, constraint_3);
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

  // Set initial conditions for your nodes, if necessary. 
  // Set spring rest length
  damping_constant = 1.0/graph.num_nodes();

  for(NodeIter it = graph.node_begin(); it != graph.node_end(); ++it) {
    Node it_ = *it;
    it_.value().vel = Point(0,0,0);
    it_.value().mass = 1.0/graph.num_nodes();
  }

  for(NodeIter it = graph.node_begin(); it != graph.node_end(); ++it) {
    Node node_ = *it;
    for(auto et = node_.edge_begin(); et != node_.edge_end(); ++et) {
      Edge rev_edge = *et;
      rev_edge.value().rest_length = rev_edge.length();
      rev_edge.value().spring_constant = 100;
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
      //double dt = 0.001 / 2;
      double dt = 0.001;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        // Combine Forces acting on graph
        auto mcf =  make_combined_force(GravityForce(), MassSpringForce(), DampingForce()); 
        
        // Combine constraints acting on graph
        auto consts = make_combined_constraints(quarter_constraint(), constant_constraint(), sphere_constraint());
        
        // Symplectic Euler
        symp_euler_step(graph, t, dt, mcf, consts);
          
        // Clear / update viewer 
        viewer.clear();
        node_map.clear();

        // Add a splash of colour
        // Lambda function that determines color scheme
        auto color_function = [&](const Node& a) {
          double ratio = cos(a.position().z + a.position().y);
          return CME212::Color::make_heat( ratio );        
        };
  
        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), color_function, node_map);
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
