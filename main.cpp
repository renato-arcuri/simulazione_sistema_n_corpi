#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace pf {

struct Config {
  double G{6.67e-11};
  double dt{0.001};
  double softening{1e-12};
};

struct Vec {
  double x;
  double y;
};

Vec operator+(const Vec &v1, const Vec &v2) {
  return Vec{v1.x + v2.x, v1.y + v2.y};
}
Vec operator-(const Vec &v1, const Vec &v2) {
  return Vec{v1.x - v2.x, v1.y - v2.y};
}
Vec operator*(const Vec &v, double c) { return Vec{v.x * c, v.y * c}; }
Vec operator*(double c, const Vec &v) { return v * c; }

double norm2(const Vec &v) { return v.x * v.x + v.y * v.y; }
double norm(const Vec &v) { return std::sqrt(norm2(v)); }

class Particle {
  Vec position_;
  Vec velocity_;
  Vec acceleration_;
  double mass_;

 public:
  Particle(Vec p, Vec v, double m)
      : position_{p}, velocity_{v}, acceleration_{0., 0.}, mass_{m} {
    if (mass_ <= 0.) {
      throw std::invalid_argument{"mass must be positive"};
    }
  }

  const Vec &getPosition() const { return position_; }
  const Vec &getVelocity() const { return velocity_; }
  const Vec &getAcceleration() const { return acceleration_; }
  double getMass() const { return mass_; }

  void setAcceleration(const Vec &a) { acceleration_ = a; }

  void computePosition(const Config &cfg) {
    position_ = position_ + getVelocity() * cfg.dt +
                getAcceleration() * (0.5 * cfg.dt * cfg.dt);
  }

  void computeVelocity(const Vec &a_old, const Config &cfg) {
    velocity_ = velocity_ + 0.5 * cfg.dt * (a_old + acceleration_);
  }
};

Vec accelerationContribution(const Particle &a, const Particle &b,
                             const Config &cfg) {
  Vec r{b.getPosition() - a.getPosition()};
  double r2{norm2(r) + cfg.softening * cfg.softening};
  double denom{r2 * std::sqrt(r2)};
  return (cfg.G * b.getMass() / denom) * r;
}

Vec totalAcceleration(const std::vector<Particle> &system, const Particle &p,
                      const Config &cfg) {
  Vec sum{0., 0.};
  for (const auto &other : system) {
    if (&other != &p) {
      sum = sum + accelerationContribution(p, other, cfg);
    }
  }
  return sum;
}

void step(std::vector<Particle> &system, const Config &cfg) {
  for (auto &p : system) {
    p.computePosition(cfg);
  }

  for (auto &p : system) {
    Vec a_old{p.getAcceleration()};
    p.setAcceleration(totalAcceleration(system, p, cfg));
    p.computeVelocity(a_old, cfg);
  }
}

double kineticEnergy(const std::vector<Particle> &system) {
  double K{0.};
  for (const auto &p : system) {
    K += 0.5 * p.getMass() * norm2(p.getVelocity());
  }
  return K;
}

double potentialEnergy(const std::vector<Particle> &system, const Config &cfg) {
  double U{0.};
  auto const N{system.size()};
  for (std::size_t i{0}; i < N; ++i) {
    for (std::size_t j{i + 1}; j < N; ++j) {
      Vec r{system[j].getPosition() - system[i].getPosition()};
      double d{norm(r)};
      U += cfg.G * system[i].getMass() * system[j].getMass() / d;
    }
  }
  return -U;
}

double mechanicalEnergy(const std::vector<Particle> &system,
                        const Config &cfg) {
  return kineticEnergy(system) + potentialEnergy(system, cfg);
}

}  // namespace pf

int main() {
  pf::Config cfg{};

  constexpr double m{1.0};
  std::vector<pf::Particle> system{
      pf::Particle{{-0.97000436, 0.24308753}, {0.4662036850, 0.4323657300}, m},
      pf::Particle{{0.97000436, -0.24308753}, {0.4662036850, 0.4323657300}, m},
      pf::Particle{{0., 0.}, {-0.93240737, -0.86473146}, m},
  };

  for (auto &p : system) {
    p.setAcceleration(pf::totalAcceleration(system, p, cfg));
  }

  constexpr int N_STEPS{10000};
  constexpr int PRINT_EVERY{500};

  std::cout << "step\tE_tot\t\tK\t\tU\n";
  for (int s{0}; s < N_STEPS; ++s) {
    pf::step(system, cfg);
    if (s % PRINT_EVERY == 0) {
      double K{pf::kineticEnergy(system)};
      double U{pf::potentialEnergy(system, cfg)};
      std::cout << s << "\t" << K + U << "\t" << K << "\t" << U << "\n";
    }
  }
}