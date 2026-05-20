// #include "progetto.hpp"

#include <algorithm>
#include <cmath>
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
  double vx;
  double vy;
};

struct Acceleration {
  double ax;
  double ay;
};

class Particle {
  Position p_;
  Velocity v_;
  double m_;

  Acceleration a_;

 public:
  Particle(Position p, Velocity v, double m) : p_{p}, v_{v}, m_{m} {}

  void setAcceleration(Acceleration a, int i) {}

  void stepPosition(Particle& a) {
    a.p_.x =
        a.p_.x + a.v_.vx * TIME_STEP + a.a_.ax * TIME_STEP * TIME_STEP * 0.5;
    a.p_.y =
        a.p_.y + a.v_.vy * TIME_STEP + a.a_.ay * TIME_STEP * TIME_STEP * 0.5;
  }

  Acceleration stepAcceleration(Particle a, Particle b) {
    Acceleration result;
    double distx = b.p_.x - a.p_.x;
    double disty = b.p_.y - a.p_.y;

    result.ax = (G * b.m_) * (-distx) /
                pow(((distx * distx) + SOFTENING * SOFTENING), 1.5);
    result.ay = (G * b.m_) * (-disty) /
                pow(((disty * disty) + SOFTENING * SOFTENING), 1.5);

    return result;
  }
};

auto operator!=(Particle a, Particle b) {
  if (a.p_ != b.p_ || a.v_ != b.v_ || a.m_ != b.m_ || a.a_ != b.a_) return true;
}

void step(std::vector<Particle>& system) {
  auto N{system.size()};

  for (auto& s : system) {
    Particle::stepPosition(s);
  }

  for (int i{1}; i != N; ++i) {
    Acceleration sum{0., 0.};
    for (int j{1}; j != N; ++j) {
      if (j != i) {
        sum += Particle::stepAcceleration;
      }
    }
    Particle::setAcceleration(sum, i);
  }
}

}  // namespace pf