#include "system.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace pf {

Vec operator+(const Vec &v1, const Vec &v2) {
  return Vec{v1.x + v2.x, v1.y + v2.y};
}
Vec operator-(const Vec &v1, const Vec &v2) {
  return Vec{v1.x - v2.x, v1.y - v2.y};
}
Vec operator*(double c, const Vec &v) { return Vec{v.x * c, v.y * c}; }
Vec operator*(const Vec &v, double c) { return Vec{v.x * c, v.y * c}; }

double norm2(const Vec &v) { return v.x * v.x + v.y * v.y; }
double norm(const Vec &v) { return std::sqrt(norm2(v)); }

Particle::Particle(Vec p, Vec v, double m)
    : position_{p}, velocity_{v}, acceleration_{0., 0.}, mass_{m} {
  if (mass_ <= 0.) {
    throw std::invalid_argument{"mass must be positive"};
  }
}

const Vec &Particle::getPosition() const { return position_; }
const Vec &Particle::getVelocity() const { return velocity_; }
const Vec &Particle::getAcceleration() const { return acceleration_; }
double Particle::getMass() const { return mass_; }

void Particle::setAcceleration(const Vec &a) { acceleration_ = a; }

void Particle::computePosition(const Config &cfg) {
  position_ =
      position_ + velocity_ * cfg.dt + acceleration_ * 0.5 * cfg.dt * cfg.dt;
}

void Particle::computeVelocity(const Vec &a_old, const Config &cfg) {
  velocity_ = velocity_ + 0.5 * cfg.dt * (a_old + acceleration_);
}

Vec Particle::momentum() const { return velocity_ * mass_; }

double Particle::angularMomentum() const {
  return mass_ * (position_.x * velocity_.y - position_.y * velocity_.x);
}

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

Vec totalMomentum(const std::vector<Particle> &system) {
  Vec M{0., 0.};
  for (const auto &s : system) {
    M = M + s.momentum();
  }
  return M;
}

double totalAngularMomentum(const std::vector<Particle> &system) {
  double L{0.};
  for (const auto &s : system) {
    L = L + s.angularMomentum();
  }
  return L;
}

auto read_from(std::filesystem::path const &p) {
  std::ifstream is(p);
  if (!is) {
    throw std::filesystem::filesystem_error{
        "read_from", p, std::make_error_code(std::errc::invalid_argument)};
  }

  std::vector<Particle> particles;
  double mass;
  double rx;
  double ry;
  double vx;
  double vy;

  while (is >> mass >> rx >> ry >> vx >> vy) {
    Vec position{rx, ry};
    Vec velocity{vx, vy};

    particles.push_back(position, velocity, mass);
  }

  if (!is.eof() && is.fail()) {
    throw std::runtime_error("Reading from file error");
  }

  return particles;
}

}  // namespace pf