// #include "progetto.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace pf {
constexpr double G{6.67e-11};
constexpr double TIME_STEP{0.001};
constexpr double SOFTENING{1e-12};

struct vec {
  double x;
  double y;
};

auto operator+(const vec &v1, const vec &v2) {
  return vec{v1.x + v2.x, v1.y + v2.y};
}

auto operator*(const vec &v, double c) { return vec{v.x * c, v.y * c}; }

auto operatorpow(const vec &v, double c) {
  return vec{pow(v.x, c), pow(v.y, c)};
}

double dt = TIME_STEP;

class Particle {
  vec position_;
  vec acceleration_;
  vec velocity_;
  double mass_;

public:
  Particle(vec p, vec v, double m) : position_{p}, velocity_{v}, mass_{m} {}

  void setAcceleration(vec a) { acceleration_ = a; }

  void computePosition() {
    position_ = position_ + (velocity_ * dt) + (acceleration_ * dt * dt * 0.5);
  }

  vec computeAcceleration(const Particle &a, const Particle &b) {
    vec r{b.position_.x - a.position_.x, b.position_.y - a.position_.y};

    double r2 = r.x * r.x + r.y * r.y + SOFTENING * SOFTENING;
    double denom = std::pow(r2, 1.5);

    vec result{r * (G * b.mass_ / denom)};
    return result;
  }

  void computeVelocity() {} // ,,,
};

void step(std::vector<Particle> &system) {
  auto N{system.size()};

  for (auto s : system) {
    s.computePosition();
  }

  for (auto i : system) {
    vec sum{0., 0.};
    for (auto j : system) {
      if (&j != &i) {
        sum = sum + i.computeAcceleration(i, j);
      }
    }
    i.setAcceleration(sum);
  }

  // ciclo che aggiorna le velocità  ...
}

} // namespace pf