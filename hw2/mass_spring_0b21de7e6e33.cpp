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

#include <vector>


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// Spring constant
static constexpr double K = 100;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edge */
struct EdgeData {
  double K;
  double L;     //< Node mass
  EdgeData(double K, double L) : K(K), L(L) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

using NodeIter  = typename GraphType::node_iterator;
using EdgeIter  = typename GraphType::edge_iterator;
using IncIter  = typename GraphType::incident_iterator;


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
  int ii = 0;
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Constraint
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      continue;
    }

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    ii ++;
  }

  return t + dt;
}

class Force {
  public:
    Point operator()(Node n, double t) {
      return force(n, t);
    };
  private:
    // user vertual method so that (Force*) p pointing to descendant objects will have
    // corresponding force method
    virtual Point force(Node n, double t) {
      return Point(0, 0, 0);
    }
};

class MassSpringForce : public Force {
  private:
    Point force(Node n, double t) {
      Point f = Point(0, 0, 0);
      double dist = 0;
      Point n1, n2;
      int ii = 0;
      for (IncIter it = n.edge_begin(); it != n.edge_end(); ++it) {
        n1 = (*it).node1().position(), n2 = (*it).node2().position();
        dist = (*it).length();
        f += - (*it).value().K * (n1 - n2) * (dist - (*it).value().L) / dist;
        ii += 1;
      }
      return f;
    }
};

class GravityForce : public Force {
  private:
    Point force(Node n, double t) {
      return n.value().mass * Point(0, 0, -grav);
    }
};

class DampingForce : public Force {
  public:
    DampingForce(double C) : C(C) {};
  private:
    double C;
    Point force(Node n, double t) {
      return -C * n.value().vel;
    }
};

struct make_combined_force {
  Force* f1;
  Force* f2;
  Force* f3;
  unsigned int num_forces;
  make_combined_force(Force* f1, Force* f2) : f1(f1), f2(f2), num_forces(2) {
  };
  make_combined_force(Force* f1, Force* f2, Force* f3) : f1(f1), f2(f2), f3(f3), num_forces(3) {};
  Point operator()(Node n, double t) {
    Point f;
    if (num_forces == 2) {
      f = (*f1)(n, t) + (*f2)(n, t);
    } else {
      f = (*f1)(n, t) + (*f2)(n, t) + (*f3)(n, t);
    }
    return f;
  }
};

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
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }
    Point f = n.value().mass * Point(0, 0, -grav);
    double dist = 0;
    for (IncIter it = n.edge_begin(); it != n.edge_end(); ++it) {
      Point n1 = (*it).node1().position(), n2 = (*it).node2().position();
      dist = norm(n1 - n2);
      f += - K * (n1 - n2) / dist * (dist - (*it).length());
    }
    return f;
  }
};

/** Constraints for HW2 #4 **/
class Constraint {
  private:
    virtual bool violate(NodeIter n) = 0;
    virtual void reset(NodeIter n) = 0;
  public:
    void operator()(GraphType& g) {
      for (NodeIter it = g.node_begin(); it != g.node_end(); ++it) {
        if (violate(it)) {
          reset(it);
        }
      }
    };
};

class PlaneConstraint : public Constraint {
  public:
    PlaneConstraint(double z) : z(z) {};
  private:
    double z;
    bool violate(NodeIter n) {
      return (*n).position().z < z;
    };
    void reset(NodeIter n) {
      (*n).position().z = z;
      (*n).value().vel.z = 0;
    }
};

class SphereConstraint : public Constraint {
  public:
    SphereConstraint(Point c, double r) : c(c), r(r) {};
  private:
    Point c;
    double r;
    bool violate(NodeIter n) {
      return norm((*n).position() - c) < r;
    };
    void reset(NodeIter n) {
      Point rd = (*n).position() - c;
      Point R = rd  / norm(rd);
      (*n).position() = c + R * r;
      (*n).value().vel = (*n).value().vel - dot((*n).value().vel, R) * R ;
    }
};

class SphereRemoveConstraint : public Constraint {
  public:
    SphereRemoveConstraint(Point c, double r) : c(c), r(r) {};
    void operator()(GraphType& g) {
      int i = 0;
      int num_nodes = g.num_nodes();
      for (NodeIter it = g.node_begin(); i < num_nodes; i ++ ) {
        if (violate(it)) {
          it = g.remove_node(it);
        } else {
          ++it;
        }
      }
    };
  private:
    Point c;
    double r;
    bool violate(NodeIter n) {
      return norm((*n).position() - c) < r;
    };
    virtual void reset(NodeIter n) {
      return;
    };
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
    graph.add_edge(nodes[t[0]], nodes[t[1]], 
        EdgeData(K, norm(nodes[t[0]].position() - nodes[t[1]].position())));
    graph.add_edge(nodes[t[0]], nodes[t[2]], 
        EdgeData(K, norm(nodes[t[0]].position() - nodes[t[2]].position())));
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]], 
        EdgeData(K, norm(nodes[t[0]].position() - nodes[t[3]].position())));
    graph.add_edge(nodes[t[1]], nodes[t[2]], 
        EdgeData(K, norm(nodes[t[1]].position() - nodes[t[2]].position())));
    graph.add_edge(nodes[t[1]], nodes[t[3]], 
        EdgeData(K, norm(nodes[t[1]].position() - nodes[t[3]].position())));
    graph.add_edge(nodes[t[2]], nodes[t[3]], 
        EdgeData(K, norm(nodes[t[2]].position() - nodes[t[3]].position())));
  }

struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};
  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  double m = 1.0 / graph.size();
  double C = 1.0 / graph.size();
  for (NodeIter it = graph.node_begin(); it != graph.node_end(); ++it) {
    NodeData nd;
    (*it).value().vel = Point(0, 0, 0);
    (*it).value().mass = m;
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

  PlaneConstraint pc(-0.75);
  SphereConstraint sc(Point(0.5, 0.5, -0.5), 0.15);
  SphereRemoveConstraint src(Point(0.5, 0.5, -0.5), 0.15);

  auto sim_thread = std::thread([&](){

      // Begin the mass-spring simulation
      double dt = 0.001;
      double t_start = 0;
      double t_end = 50.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {
        GravityForce gf;
        MassSpringForce msf;
        DampingForce df(C);
        Force* gfp = &gf;
        Force* msfp = &msf;
        Force* dfp = &df;
        symp_euler_step(graph, t, dt, 
          make_combined_force(gfp, msfp, dfp));
        //sc(graph);
        //pc(graph);
        src(graph);

        // Clear the viewer’s nodes and edges
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
