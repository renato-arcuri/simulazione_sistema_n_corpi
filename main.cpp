// #include "progetto.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace pf {

// costanti
struct Config {
  constexpr double G{6.67e-11};
  constexpr double DT{0.001};
  constexpr double SOFTENING{1e-12};
}

struct Vec {
  double x;
  double y;
};

auto operator+(const Vec &v1, const Vec &v2) {
  return Vec{v1.x + v2.x, v1.y + v2.y};
}
auto operator*(const Vec &v, double c) { return Vec{v.x * c, v.y * c}; }
auto norm2(Vec &v) { return v.x * v.x + v.y * v.y; }

class Particle {
  Vec position_;
  Vec acceleration_;
  Vec velocity_;
  double mass_;

 public:
  Particle(Vec p, Vec v, double m) : position_{p}, velocity_{v}, mass_{m} {}

  double getMass() { return mass_; }
  auto getPosition() { return position_; }
  auto getVelocity() { return velocity_; }
  auto getAcceleration() { return acceleration_; }

  void setAcceleration(Vec a) { acceleration_ = a; }
  void setVelocity(Vec v) { velocity_ = v; }

  void computePosition() {
    getPosition() = getPosition() + (getVelocity() * DT) +
                    (getAcceleration() * DT * DT * 0.5);
  }

  Vec computeAcceleration(const Particle &a, const Particle &b, Config cost) {
    Vec r{b.getPosition().x - a.getPosition().x,
          b.getPosition().y - a.getPosition().y};

    double r2 = norm2(r) + cost.SOFTENING * cost.SOFTENING;
    double denom = pow(r2, 1.5);

    Vec result{r * (G * b.getMass() / denom)};
    return result;
  }

  void computeVelocity(Vec a_new, Config cost) {
    getVelocity() = getVelocity() + 0.5 * cost.DT * (getAcceleration() + a_new);
  }
};

void step(std::vector<Particle> &system) {
  // ricalcolo posizioni
  for (auto s : system) {
    s.computePosition();
  }

  // ricalcolo accelerazioni, ricalcolo velocità
  for (auto i : system) {
    vec sum{0., 0.};
    for (auto j : system) {
      if (&j != &i) {
        sum += i.computeAcceleration(i, j);
      }
    }
    i.computeVelocity(sum, Config);
    i.setAcceleration(sum);
  }
}

double KineticEnergy(std::vector<Particle> &system) {
  double K{0.};
  for (auto p : system) {
    K += 0.5 * p.getMass() * norm2(p);
  }

  return K;
}

double PotentialEnergy(std::vector<Particle> &system, Config cost) {
  double U{0.};
  for (auto i : system) {
    for (auto j : system) {
      if (&j != &i) {
        r = Vec{j.getPosition().x - i.getPosition().x,
                j.getPosition().y - i.getPosition().y};
        double d{sqrt(norm2(r))};
        U += (i.getMass() * j.getMass() * cost.G) / d;
        return -U;
      }
    }
  }

  double MechanicalEnergy(std::vector<Particle> & system) {
    return KineticEnergy(system) + PotentialEnergy(system);
  }

}  // namespace pf

int main() {
  /*
  std::vector<Particle> system;
  int nStep{};

  Vec p1 = {0.1, 0.2};
  Vec p2 = {2.0, 3.0};
  Vec p3 = {4.0, 9.2};

  Vec v1{2.3, 0.1};
  Vec v2{8.4, 3.4} Vec v3{3.4, 2.0}

  double m1{1.};
  double m2{5.};
  double m3{2.}

  Particle c1= {p1, m1, v1}

      Particle c2 Particle c3 {}
  system.push_back()
  */
}