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

struct EdgeData
{
    double k;
    double length;

    EdgeData() : k(1), length(.1){
       // std::cout << "hello " << std::endl;
    }
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

//WITH Constraints problem 3
template<typename G, typename F>
double symp_euler_step_p3(G& g, double t, double dt, F force) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        //Skip nodes (0,0,0) and (1,0,0)
        if ((*it).position()==Point(0,0,0) or (*it).position()==Point(1,0,0))
        {
            continue;
        }
        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        else
        n.position() += n.value().vel * dt;
    }

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        if ((*it).position()==Point(0,0,0) or (*it).position()==Point(1,0,0))
        {
            continue;
        }
        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }

    return t + dt;
}

template <typename G, typename F, typename C1>
double symp_euler_step_p4(G& g, double t, double dt, F force, C1 cons) {
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
   // constraint(-0.75);
    cons(g,t);
    //CombinedConstraint(constraint1,constraint2);
    return t + dt;
}

template <typename G, typename F, typename C>
double symp_euler_step_p5(G& g, double t, double dt, F force, C constraint) {
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
    }

    return t + dt;
}

struct Problem1Force {
    /** Return the force applying to @a n at time @a t.
     *
     * For HW2 #1, this is a combination of mass-spring force and gravity,
     * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
     * model that by returning a zero-valued force. */

    Problem1Force(double kcons, double l):k(kcons),length(l){};

    template <typename NODE>
    Point operator()(NODE n, double t) {
        // HW2 #1: YOUR CODE HERE
        (void) t;    // silence compiler warnings

        Point SpringF(0,0,0);
        Point xi = n.position();

        // Returning zero force in corners

        if(xi == Point(0,0,0) or xi == Point(1,0,0))
        {
            return Point(0,0,0);
        }


        for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit )
        {
            Point xj = (*eit).node2().position();
            SpringF += (k) * (xi - xj) / (norm(xi - xj)) * (norm(xi - xj) - length);
        }

        Point Ftotal = SpringF + Point(0,0,-grav)*n.value().mass;
        return Ftotal;
    }
    double k,length;
};
/** Force function object for HW2 #1. */
struct Problem2Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */


    template <typename NODE>
    Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;    // silence compiler warnings

    Point SpringF(0,0,0);
    Point xi = n.position();

    // Returning zero force in corners

    if(xi == Point(0,0,0) or xi == Point(1,0,0))
    {
      return Point(0,0,0);
    }


    for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit )
    {
      Point xj = (*eit).node2().position();
      SpringF += (-(*eit).value().k) * (xi-xj)/(norm(xi-xj))*(norm(xi-xj)-(*eit).value().length);
        /*
        std::cout << "k: " << (-(*eit).value().k) << std::endl;
        std::cout << "delta: " << (norm(xi-xj)-(*eit).value().length) << std::endl;
        std::cout << "norm: " << norm(xi-xj) << std::endl;
         */
    }
        //std::cout<< SpringF.x << " " << SpringF.y << " " << SpringF.z << std::endl;

     Point Ftotal = SpringF + Point(0,0,-grav)*n.value().mass;
     return Ftotal;
  }
  double k,length;
};

struct GravityForce
{
    template<typename Node>
    Point operator()(Node n, double t)
    {
        (void) t;
        Point gravityF = Point(0,0,-grav)*n.value().mass;
        return gravityF;
    }
};

struct SpringForce
{
    template<typename Node>


    Point operator()(Node n, double t)
    {
        (void) t;

        Point springF(0,0,0);
        Point xi = n.position();

        for (auto eit = n.edge_begin(); eit != n.edge_end(); ++eit )
        {
            Point xj = (*eit).node2().position();
            springF += (-(*eit).value().k) * (xi-xj)/(norm(xi-xj))*(norm(xi-xj)-(*eit).value().length);

        }

        return springF;
    }
};

struct DampingForce
{
    DampingForce(double damping):c(damping){};

    template<typename Node>
    Point operator()(Node n, double t)
    {
        (void) t;

        Point DampF = -c*n.value().vel;
        return DampF;
    }
    double c;
};

template<typename F1, typename F2>

struct CombineForce
{
    F1 f1;
    F2 f2;

    CombineForce(F1 force1, F2 force2):f1(force1),f2(force2){};

    template<typename Node>

    Point operator()(Node n, double t)
    {
        (void) t;
        return f1(n,t) + f2(n,t);
    }
};

struct ConstraintPlane
{
    ConstraintPlane(double inputz) : zplane(inputz) {};

    template<typename Graph>
    void operator()(Graph &graph, double t)
    {
        (void) t;

        for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit)
        {
            auto currnode = (*nit);

            if (currnode.position().z < zplane)
            {
                currnode.position().z = zplane;
                currnode.value().vel.z = 0;
            }
        }
    }

    double zplane;
};

struct ConstraintSphere
{
    ConstraintSphere(Point cen, double rad):center(cen),radius(rad){};

    template<typename Graph>
    void operator()(Graph& graph, double t)
    {
        (void) t;

        for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit)
        {
            auto currnode = (*nit);
            Point x = currnode.position();

            Point R = (x-center)/norm(x-center);

            if ( norm(currnode.position() - center) < radius)
            {
               // std::cout << "center " << center << "radius" << radius << std::endl;
                currnode.position()=radius*R+center;
                currnode.value().vel -= (currnode.value().vel*R)*R;
            }
        }

    }

    Point center;
    double radius;
};

struct ConstraintSphereRemove
{
    ConstraintSphereRemove(Point cen, double rad):center(cen),radius(rad){};

    template<typename Graph>
    void operator()(Graph& graph, double t)
    {
        (void) t;

        for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit)
        {
            auto currnode = (*nit);
            //auto targetid = currnode.index();
            Point x = currnode.position();

            //Point R = (x-center)/norm(x-center);

            if ( norm(currnode.position() - center) <= radius)
            {
               auto it= graph.remove_node((*nit));
                (void) it;
               continue;
            }
        }
    }
    Point center;
    double radius;
};

template<typename F1, typename F2>
CombineForce<F1,F2> make_combined_force(F1 f1, F2 f2)
{
    return CombineForce<F1,F2>(f1,f2);
}

template<typename F1, typename F2, typename F3>
CombineForce<CombineForce<F1,F2>,F3> make_combined_force(F1 f1,F2 f2,F3 f3)
{
    auto Force1and2 = make_combined_force(f1,f2);
    return CombineForce<CombineForce<F1,F2>,F3>(Force1and2,f3);
};

template<typename C1, typename C2>
struct CombinedConstraints
{
    C1 c1;
    C2 c2;

    CombinedConstraints(C1 constr1, C2 constr2):c1(constr1),c2(constr2){};

    template<typename Graph>
    void operator()(Graph& graph, double t)
    {
        c1(graph, t);
        c2(graph, t);
    }
};

template<typename C1, typename C2, typename C3>
CombinedConstraints<CombinedConstraints<C1,C2>,C3> make_combined_constraints(C1 c1,C2 c2,C3 c3)
{
    auto c1and2 = make_combined_constraints(c1,c2);
    return CombinedConstraints<CombinedConstraints<C1,C2>,C3>(c1and2,c3);
};

template<typename C1,typename C2>
CombinedConstraints<C1,C2> make_combined_constraints(C1 c1, C2 c2)
{
    return CombinedConstraints<C1,C2>(c1,c2);
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
#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.

  for (auto nit = graph.node_begin(); nit!=graph.node_end();++nit)
  {
    (*nit).value().mass=1.0/(double)graph.num_nodes();
  }

  double len = norm((*(graph.edge_begin())).node2().position()-(*(graph.edge_begin())).node1().position());

  for (auto eit = graph.edge_begin(); eit!=graph.edge_end();++eit)
  {
     // std::cout << "starting setting value" << std::endl;
      EdgeData& Edgeinfo = (*eit).value();
      Edgeinfo.length = len;
      Edgeinfo.k = 100.0;
     // std::cout << "k = " << (*eit).value().k << "; len = " << (*eit).value().length << "; " << std::endl;
  }


/*
  for (auto nit = graph.node_begin(); nit!=graph.node_end();++nit)
  {
      for (auto aeit=(*nit).edge_begin();aeit!=(*nit).edge_end();++aeit)
      {
          //EdgeData& Edgeinfo = (*aeit).value();
          (*aeit).value().length = len;
          (*aeit).value().k = 100.0;
          std::cout << "k = " << (*aeit).value().k << "; len = " << (*aeit).value().length << "; " << std::endl;

      }
  }
*/
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
      double dt = 0.0005;
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        //std::cout << "t = " << t << std::endl;

        //symp_euler_step(graph, t, dt, Problem1Force(100,.010101));
        //symp_euler_step(graph, t, dt, Problem2Force());
        //symp_euler_step_p3(graph, t, dt, make_combined_force(GravityForce(),SpringForce()));
        //symp_euler_step_p3(graph, t, dt, make_combined_force(GravityForce(),SpringForce(),DampingForce(1/graph.num_nodes())));
        //symp_euler_step_p4(graph, t, dt, Problem2Force(),ConstraintPlane(-.75));
        symp_euler_step_p4(graph, t, dt, Problem2Force(),make_combined_constraints(ConstraintSphere(Point(0.5,0.5,-0.5),.15),ConstraintPlane(-.75)));
        // symp_euler_step_p5(graph, t, dt, Problem2Force(),ConstraintSphereRemove(Point(0.5,0.5,-0.5),.15));

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
