// #include "progetto.hpp"

#include <algorithm>
#include <vector>

namespace pf {
constexpr double G{6.67e-11};
constexpr double TIME_STEP{0.001};
constexpr double SOFTENING{1e-12};

struct Position {
  double x;
  double y;
};

struct Velocity {
  double v_x;
  double v_y;
};

struct Acceleration {
  double a_x;
  double a_y;
};

class Particle {
  Position p_;
  Velocity v_;
  Acceleration a_;
  double m_;

 public:
  Particle(Position p, Velocity v, double m)
      : p_{p}, v_{v}, a_{0., 0.}, m_{m} {}
};

std::vector<Particle> system{};

void step() {}

}  // namespace pf