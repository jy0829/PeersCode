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
  double c;
  NodeData() : vel(0), mass(1), c(0) {}
};

// HW2
struct EdgeData {
    double K;
    double L;
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

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraint(g, t);
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
    (void) t;
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        return Point(0,0,0);
    }
    // start with gravity
    Point F_total = Point(0,0,-n.value().mass*grav);
    // iterate over incident edges and add spring force
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        Edge e = *it;
        F_total += -e.value().K*(e.node1().position()-e.node2().position())/e.length()*(e.length()-e.value().L);
    }
    return F_total;
  }
};

struct GravityForce {
    Point operator()(Node n, double t) {
        (void) t;
        return Point(0,0,-n.value().mass*grav);
    }
};

struct MassSpringForce {
    Point operator()(Node n, double t) {
        (void) t;
        Point F_total;
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
            Edge e = *it;
            F_total += -e.value().K*(e.node1().position()-e.node2().position())/e.length()*(e.length()-e.value().L);
        }
    return F_total;
    }
};

struct DampingForce {
    Point operator()(Node n, double t) {
        (void) t;
        return -n.value().c*n.value().vel;
    }
};

struct ZeroForce {
    Point operator()(Node n, double t) {
        (void) t;
        (void) n;
        return Point(0,0,0);
    }
};

template <typename F1, typename F2, typename F3>
struct TripleForce {
    F1 f1_;
    F2 f2_;
    F3 f3_;
    TripleForce(F1 f1, F2 f2, F3 f3): f1_(f1), f2_(f2), f3_(f3) {}
    Point operator()(Node n, double t) {
        return f1_(n,t) + f2_(n,t) + f3_(n,t);
    }
};

template <typename F1, typename F2>
TripleForce<F1,F2,ZeroForce> make_combined_force(F1 f1, F2 f2) {
    return TripleForce<F1,F2,ZeroForce>(f1, f2, (ZeroForce()));
}
template <typename F1, typename F2, typename F3>
TripleForce<F1,F2,F3> make_combined_force(F1 f1, F2 f2, F3 f3) {
    return TripleForce<F1,F2,F3>(f1, f2, f3);
}

struct FixedNodeConstraint {
    unsigned node_uid;
    Point fixed_pos;
    FixedNodeConstraint(GraphType& g, Point fixed) : node_uid(0), fixed_pos(fixed) {
        // find first node with position equal to fixed position
        for (auto it=g.node_begin(); it!=g.node_end(); ++it) {
            if ((*it).position() == fixed_pos) {
                node_uid = (*it).index();
                break;
            }
        }
    }
    void operator()(GraphType& g, double t) {
        (void) t;
        g.node(node_uid).position() = fixed_pos;
        g.node(node_uid).value().vel = Point(0,0,0);
    }
};

struct NullConstraint {
    void operator()(GraphType& g, double t) {
        (void) g;
        (void) t;
    }
};

struct PlaneConstraint {
    Point normal_;
    double v_;
    // @pre normal is a unit vector
    PlaneConstraint(Point normal, double v): normal_(normal), v_(v) {}
    void operator()(GraphType& g, double t) {
        (void) t;
        for (auto it=g.node_begin(); it!=g.node_end(); ++it) {
            Node n = *it;
            if (dot(n.position(),normal_) < v_) {
                // constraint is violated, so set pos to be nearest point on plane, vel to zero
                Point& pos = n.position();
                Point pointOnPlane = v_*normal_;
                pos = (pos - normal_*dot(pos-pointOnPlane,normal_));
                Point& vel = n.value().vel;
                vel = vel - normal_*dot(vel,normal_);
            }
        }
    }
};

struct SphereConstraint {
    Point center_;
    double r_;
    SphereConstraint(Point center, double r): center_(center), r_(r) {}
    void operator()(GraphType& g, double t) {
        (void) t;
        for (auto it=g.node_begin(); it!=g.node_end(); ++it) {
            Node n = *it;
            if (norm(n.position()-center_) < r_) {
                /* HW 4 */
                /*
                // constraint violated, set pos to nearest point on sphere
                Point& pos = n.position();
                pos = pos - center_;        // shift relative to sphere center
                pos = pos * r_/norm(pos);    // scale to length=r
                pos = pos + center_;        // shift back
                // velocity
                Point& vel = n.value().vel;
                Point normal = (pos-center_)/norm(pos-center_); // unit normal for sphere is position
                vel = vel - dot(vel,normal)*normal;
                */
                /* HW 5 */
                g.remove_node(n);
            }
        }
    }

};

template <typename C1, typename C2, typename C3, typename C4>
struct QuadConstraint {
    C1 c1_;
    C2 c2_;
    C3 c3_;
    C4 c4_;
    QuadConstraint(C1 c1, C2 c2, C3 c3, C4 c4): c1_(c1), c2_(c2), c3_(c3), c4_(c4) {}
    void operator()(GraphType& g, double t) {
        c1_(g, t);
        c2_(g, t);
        c3_(g, t);
        c4_(g, t);
    }
};

template <typename C1, typename C2, typename C3, typename C4>
QuadConstraint<C1,C2,C3,C4> make_combined_constraint(C1 c1, C2 c2, C3 c3, C4 c4) {
    return QuadConstraint<C1,C2,C3,C4>(c1, c2, c3, c4);
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
  for (auto it = graph.node_begin(); it!=graph.node_end(); ++it) {
      Node n = *it;
      n.value().mass = 1.0/graph.size();
      n.value().c = 1.0/graph.size();
  }
  for (auto it = graph.edge_begin(); it!=graph.edge_end(); ++it) {
      Edge e = *it;
      e.value().L = e.length();
      e.value().K = 100;
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
      double dt = 0.0005; // smaller timestep needed for grid3
      double t_start = 0;
      double t_end = 5.0;

      for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {

        symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(), DampingForce()),
                make_combined_constraint(FixedNodeConstraint(graph, Point(0,0,0)),
                                         FixedNodeConstraint(graph, Point(1,0,0)),
                                         PlaneConstraint(Point(0,0,1), -0.75),
                                         SphereConstraint(Point(0.5,0.5,-0.5), 0.15)) );
        // clear viewer
        viewer.clear();
        node_map.clear();

        // Update viewer with nodes' new positions
        viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
        viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
        viewer.set_label(t);

        // These lines slow down the animation for small graphs, like grid0_*.
        // Feel free to remove them or tweak the constants.
        if (graph.size() < 1000)
          std::this_thread::sleep_for(std::chrono::milliseconds(5));
      }

    });  // simulation thread

  viewer.event_loop();

  // If we return from the event loop, we've killed the window.
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
