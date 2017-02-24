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

#define DEBUX(x) \
  do { \
    std::cerr << __FILE__ << ": " << __LINE__ << ": "; \
    std::cerr << #x << "->" << (x) << std::endl; \
  } while (0)
  
   


// Gravity in meters/sec^2
static constexpr double grav = 9.81, c = 0.0001;




/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(Point(0,0,0)), mass(1) {}
};

/** Custom structure of data to store with Edges. */
struct EdgeData {
  double K;
  double L;
  EdgeData(double K_temp = 0, double L_temp = 0): K(K_temp),L(L_temp) {}
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
 * @tparam G::node_value_type supports a concept with data:
                 Point vel;
                 double mass;
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
    (void) t;
    //if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
     // return Point(0,0,0);
    //}
    Point force(0,0,0), gravity(0,0,-grav);
    force += gravity*n.value().mass;
    double l =0; 
    for (auto it = n.edge_begin();it != n.edge_end();++it) {
      auto e = *it;
      auto p = e.arrow();
      l = norm(p);
      force -= e.value().K*p*(1-e.value().L/l);
    }
    return force;
  }
};

/** Struct for computing gravitational force. */
struct GravityForce {
  /** Return the gravitational force applying to @a n at time @a t.
   * @param n: Node object for which we compute force.
   * @param t: a double, representing time.
   * @returns a Point which represents the force vector.
   */
  Point operator()(Node n, double t) {
    (void) t;
    return Point(0,0,-grav)*n.value().mass;
  }
};

/** Struct for computing total spring force. */
struct MassSpringForce {
  /** Return the spring force applying to @a n at time @a t due to its neighbors.
   * @param n: Node object for which we compute force.
   * @param t: a double, representing time.
   * @returns a Point which represents the force vector.
   */
  Point operator()(Node n, double t) {
    (void) t;
    Point force(0,0,0);
    double l =0; 
    for (auto it = n.edge_begin();it != n.edge_end();++it) {
      auto e = *it;
      auto p = e.arrow();
      l = norm(p);
      force -= e.value().K*p*(1-e.value().L/l);
    }
    return force;
  }
};

/** Struct for computing damping force. */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t.
   *  c is a global static constant.
   * @param n: Node object for which we compute force.
   * @param t: a double, representing time.
   * @returns a Point which represents the force vector.
   */
  Point operator()(Node n, double t) {
    (void) t;
    return n.value().vel*(-c);
  }
};

/** Creates a combination force struct for two forces.
* @tparam F1, @tparam F2: a force class with the following method:
    Point operator()(NODE n, double t);
*/
template <typename F1, typename F2>
struct CombinedForce {
  F1 f1;
  F2 f2;
  CombinedForce(F1 &f1_temp,F2 &f2_temp): f1(f1_temp), f2(f2_temp) {}
  /** Returns the combination force value for two forces.
   * @tparam NODE: NODE class that works with the existing force classes.
   * @param n: Node for which we compute forces.
   * @param t: time at which we compute force.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n,t) + f2(n,t);
  }  
};

template <typename F1, typename F2>
CombinedForce<F1,F2> make_combined_force(F1 f1, F2 f2) {
  return CombinedForce<F1,F2>(f1,f2);
}

template <typename F1, typename F2, typename F3>
CombinedForce<F1,CombinedForce<F2,F3>> make_combined_force(F1 f1, F2 f2, F3 f3) {
  return make_combined_force(f1,CombinedForce<F2,F3>(f2,f3));
}


/** Creates the constraint of setting a node to a fixed point.
 * @a node: Node object which we intend to keep fixed.
 * @a point: Point object which specificies the fixed position of @a node.
*/
struct ConstraintFixedNode {
  Node node;
  Point point;
  /**Sets up the values of node and the point to which it will be fixed.
  * @param n: Node for which we specify a fixed point.
  * @param p: Point at which we fix position of @a n.
  */
  ConstraintFixedNode(Node &n, Point p): node(n), point(p) {}
  

  /**Enforces fixed node constraint.
  * @param n: node on which we check constraint.
  * @param t: time of simulation.
  * @post if node == n
           node.position() == point;
           node.value().vel == Point(0,0,0);
  */
  void operator()(Node n, double t) {
    (void) t;
    if (node == n) {
      n.position() = point;
      n.value().vel = Point(0,0,0);
    }
  }
};

/** Creates the constraint of being in a half plane.
 * @a normal: The normal to the plane.
 * @a point: The value that defines the plane, @a v . normal < val;
*/
struct ConstraintPlane {
  Point normal;
  double val;
  ConstraintPlane(Point n, double v): normal(n), val(v) {}

  /**Enforces half-plane constraint.
  * @param n: node on which we check constraint.
  * @param t: time of simulation.
  * @post if dot(old node.position(),normal) > val
               dot(node.position(),normal) == val
               dot(node.value().vel,normal) == Point(0,0,0);
  */
  void operator()(Node &n, double t) {
    (void) t;
    double x = dot(n.position(),normal);
    if (x < val) {
      n.position() = n.position() + normal*(val-x)/dot(normal,normal);
      n.value().vel = n.value().vel - normal*dot(n.value().vel,normal)/dot(normal,normal);
    }
  }
};

/** Creates the constraint of being outside a sphere.
 * @a center: a Point, the center to the sphere.
 * @a radius: a double, the radius of the sphere.
*/
struct ConstraintSphere {
  Point center;
  double radius;
  ConstraintSphere(Point c, double r): center(c), radius(r) {}

  /**Enforces sphere constraint.
  * @param n: node on which we check constraint.
  * @param t: time of simulation.
  * @post if norm(old node.position()-center) < radius)
              norm(node.position()-center) == radius;
              dot(node.value().vel,node.position()-center) == Point(0,0,0);
  */
  void operator()(Node &n, double t) {
    (void) t;
    Point x = n.position()-center;
    if (norm(x) < radius) {
      n.position() = center + x*radius/norm(x);
      n.value().vel = n.value().vel - x*dot(n.value().vel,x)/dot(x,x);
    }
  }
};

/** Creates the constraint of being outside a sphere.
 * @a center: a Point, the center to the sphere.
 * @a radius: a double, the radius of the sphere.
*/
struct ConstraintSphereRemoval {
  Point center;
  double radius;
  ConstraintSphereRemoval(Point c, double r): center(c), radius(r) {}

  /**Enforces sphere constraint.
  * @param n: node on which we check constraint.
  * @param t: time of simulation.
  * @post if norm(old node.position()-center) < radius)
              node is removed using node_remove (see post-conditions of remove_node in Graph.hpp).
  */
  void operator()(Node &n, double t) {
    (void) t;
    Point x = n.position()-center;
    if (norm(x) < radius) {
      n.remove();
    }
  }
};


/** Combines constraints to create a new constraint object.
* @tparam C1,C2: constraint concepts with
              void operator()(Node, double);
*/
template<typename C1, typename C2>
struct CombinedConstraint {
  C1 c1;
  C2 c2;
  CombinedConstraint(C1 d1, C2 d2): c1(d1), c2(d2) {}
  /*Executes individual constraint functors.*/
  void operator()(Node &n, double t) {
    c1(n,t);
    c2(n,t);
  }
};

/** mcc: make combined constraint.
 * @tparam C1, C2: constraint structs.
      They have void operator()(Node n, double t);
 * @param c1: constraint object of struct C1.
 * @param c2: constraint object of struct C2.
*/
template<typename C1, typename C2>
CombinedConstraint<C1,C2> mcc(C1 c1, C2 c2) {
  return CombinedConstraint<C1,C2>(c1,c2);
}

/*Creates a constrain enforcer functor which iterates through all nodes to enforce constraint.
 * @tparam C: constraint struct.
      Has void operator()(Node n, double t);
 */
template <typename C>
struct ConstraintEnforcer {
  C c;
  ConstraintEnforcer(C d): c(d) {}
  /* Functor which enforces constraints on all nodes.
  * @param g: GraphType object for which we enforce constraints on nodes.
  * @param t: double representing time.
  * @post: post-conditions of constraint c are satisfied for all nodes.
  */
  void operator()(GraphType &g, double t) {
    for (auto it = g.node_begin(); it!= g.node_end(); ++it) {
      auto n=*it;
      c(n,t);
    }
  }
};

/*Creates a constrain enforcer functor which iterates through all nodes to enforce constraint.
 * @tparam C: constraint struct.
      Has void operator()(Node n, double t);
 */
template <typename C>
ConstraintEnforcer<C> make_constraint_enforcer(C c) {
  return ConstraintEnforcer<C>(c);
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

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  double K = 100;
  double M = 1.0/(double)(graph.num_nodes());
  Node n0, n1;
  for (auto it = graph.node_begin(); it!= graph.node_end();++it) {
    auto n = *it;
    n.value().vel = Point(0,0,0);
    n.value().mass = M;

    //Get nodes for which we put a fixed position constraint.
    if (n.position()==Point(0,0,0)) n0 = n;
    if (n.position()==Point(1,0,0)) n1 = n;
  }
  for (auto it = graph.edge_begin(); it!= graph.edge_end();++it) {
    auto e = *it;
    e.value().L = norm(e.arrow());
    e.value().K = K;
  }

  ConstraintFixedNode fixed_0(n0,Point(0,0,0));
  ConstraintFixedNode fixed_1(n1,Point(1,0,0));
  ConstraintPlane plane(Point(0,0,1),-0.75);
  ConstraintSphere sphere(Point(0.5,0.5,-0.5),0.15);
  ConstraintSphereRemoval sphere2(Point(0.5,0.5,-0.5),0.15);
  auto constraint = mcc(fixed_0,mcc(fixed_1,mcc(plane,sphere2)));
  auto enforcer = make_constraint_enforcer(constraint);
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
      GravityForce G;
      MassSpringForce MS;
      DampingForce D;
      auto F = make_combined_force(G,MS,D);
      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        symp_euler_step(graph, t, dt, F);
        viewer.clear();
        node_map.clear();
        enforcer(graph,t);
        // Update viewer with nodes' new positions
        
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 100)
          std::this_thread::sleep_for(std::chrono::microseconds(10));
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
