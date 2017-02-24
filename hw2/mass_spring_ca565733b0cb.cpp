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
    double K; // spring constant
    double L; // edge length
    EdgeData() {}
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

//added the constraint
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

    //constraint(g,t);

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

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

      //prob3 TEMP
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      { n.value().vel = Point(0, 0, 0);}

  }

    //constraint(g,t);

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
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    Point spring_force = Point(0,0,0);
    Point gravity_force = Point(0,0,-grav);

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
        auto adjacent_edge = *it;
        Point diff_points = (n.position()-adjacent_edge.node2().position());
    spring_force += -adjacent_edge.value().K*diff_points/norm_2(diff_points)*(norm_2(diff_points)-adjacent_edge.value().L);
    }
    return spring_force + n.value().mass * gravity_force;
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
   // return Point(0);
  }
};


struct GravityForce {
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        return n.value().mass*Point(0,0,-grav);
    }
};

struct MassSpringForce {
     template <typename NODE>
     Point operator()(NODE n, double t) {
        (void) t;
        Point spring_force = Point(0,0,0);
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
            {
             auto adjacent_edge = *it;
             Point diff_points = (n.position()-adjacent_edge.node2().position());
             spring_force += -adjacent_edge.value().K*diff_points/norm(diff_points)*(norm(diff_points)-adjacent_edge.value().L);
            }
        return spring_force ;
        }
};


struct DampingForce{
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void) t; //to avoid the warning
        return n.value().vel*const_dump;
    }
    double const_dump;
};

template <typename F1, typename F2>
struct Combined_force{
    F1 f1;
    F2 f2;
    template<typename NODE>
    Point operator()(NODE n, double t)
    {
        (void) t; //to avoid warning
        return f1(n,t)+f2(n,t); //call c1(g,t)
    }
};

template <typename F1, typename F2>
Combined_force<F1,F2> make_combined_force(F1 f1, F2 f2){
   return {f1,f2};
}

template <typename F1, typename F2, typename F3>
Combined_force<Combined_force<F1,F2>, F3> make_combined_force(F1 f1,F2 f2,F3 f3){
    auto miniforce = make_combined_force(f1,f2);
    return make_combined_force(miniforce,f3);
}

struct Constraint_Plane {
    //template<typename GraphType>
    void operator()(GraphType& g, double t) {
        (void) t;
        Point p(0,0,1);
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            auto n = *it;
            if (dot(n.position(),p) < -0.75) {
                n.position().z = -0.75; //nearest point
                n.value().vel.z = 0;
            }
        }
    }
};

struct Constraint_Sphere1 {

    //template<typename GraphType>
    void operator()(GraphType& g, double t)
    {
        (void) t;
        double r = 0.15;
        Point center(0.5,0.5,-0.5);
        for (auto it = g.node_begin(); it != g.node_end(); ++it)
        {
            auto n = *it;
            if (norm_2(n.position()-center)<r) {
                n.position() = r * (n.position()-center)/norm_2(n.position()-center);
                auto Ri = n.position() / r;
                n.value().vel = n.value().vel - dot(n.value().vel, Ri) * Ri;
                }
        }
     }
};


struct Constraint_Sphere2 {

    //template<typename GraphType>
    void operator()(GraphType& g, double t)
    {
        (void) t;
        double r = 0.15;
        Point center(0.5,0.5,-0.5);
        for (auto it = g.node_begin(); it != g.node_end(); ++it)
        {
            auto n = *it;
            if ((norm_2(n.position()-center)<r) || (n.degree()==0))
                // added that extra condition  for n.degree() ==0 to eliminate empty points
            {
                 g.remove_node(n);
            }
        }
    }
};


template <typename C1, typename C2, typename C3>
struct Combined_constraints{
    C1 c1;
    C2 c2;
    C3 c3;
    //template<typename NODE>
    void operator()(GraphType& g, double t)
    {
        c1(g,t);
        c2(g,t);
        c3(g,t);
    }
};

//this is useful because I can call everything via this function
//without worrying about templates
//returns an object
template <typename C1, typename C2, typename C3>
Combined_constraints<C1,C2, C3> make_combined_constraints(C1 c1, C2 c2, C3 c3){
    return {c1,c2,c3};
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
      //Hw 2 Problem 2 - ADD THE DIAGONALS
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

  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    Point p(0,0,0);
    n.value().vel = p;
    n.value().mass = 1.0/((double)graph.num_nodes());
  }

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it){
      auto e = *it;
      e.value().K = 100.0;
      e.value().L = e.length();
  }

    double const_dump = -1.0/double(graph.num_nodes());

    DampingForce damping_force{const_dump};

    GravityForce gravity_force;
    MassSpringForce mass_spring_force;

    auto CForce = make_combined_force(mass_spring_force, gravity_force, damping_force);

    //Add Problem 4 Graph Constraints
    Constraint_Plane constraint1;
    Constraint_Sphere1 constraint2;
    Constraint_Sphere2 constraint3;

    //Activate Constraint 1, 2 and 3
    auto CConditions = make_combined_constraints(constraint1, constraint3, constraint2);

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
      double dt = 0.001; //was 0.001 - will change to 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;
        //symp_euler_step(graph, t, dt, Problem1Force(), CConditions);//; HERE I SWITHCED FOR 3

        //  symp_euler_step(graph, t, dt, CForce); // Part 4
        symp_euler_step(graph, t, dt, CForce, CConditions); //Part 4 and Part 5

        //Hw 2 Problme 5
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(),graph.edge_end(), node_map);

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
