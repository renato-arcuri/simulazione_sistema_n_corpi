#include "system.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
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
    : position_{p}, velocity_{v}, acceleration_{0., 0.}, mass_{m}, k_{0.f} {
  if (mass_ <= 0.) {
    throw std::invalid_argument{"mass must be positive"};
  }
}

const Vec &Particle::getPosition() const { return position_; }
const Vec &Particle::getVelocity() const { return velocity_; }
const Vec &Particle::getAcceleration() const { return acceleration_; }
double Particle::getMass() const { return mass_; }
float Particle::getRadius() const { return k_; }

void Particle::setAcceleration(const Vec &a) { acceleration_ = a; }
void Particle::setRadius(const float &k) { k_ = k; }

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
  double denom{(std::pow(norm(r), 3.0))};
  return (cfg.G * b.getMass() / denom) * r;
}

Vec totalAcceleration(const std::vector<Particle> &system, const Particle &p,
                      const Config &cfg) {
  return std::accumulate(system.begin(), system.end(), Vec{0., 0.},
                         [&p, &cfg](const Vec &sum, const Particle &other) {
                           if (&other == &p) return sum;
                           return sum + accelerationContribution(p, other, cfg);
                         });
}

void computeRadius(std::vector<Particle> &system) {
  if (system.empty()) return;

  double min_dist{std::numeric_limits<double>::max()};
  for (std::size_t i{0}; i < system.size(); ++i) {
    for (std::size_t j{i + 1}; j < system.size(); ++j) {
      double d{norm(system[i].getPosition() - system[j].getPosition())};
      min_dist = std::min(min_dist, d);
    }
  }

  double max_radius{min_dist * 0.1};

  auto max_m_it = std::max_element(system.begin(), system.end(),
                                   [](const Particle &a, const Particle &b) {
                                     return a.getMass() < b.getMass();
                                   });
  double m_max = max_m_it->getMass();

  for (auto &p : system) {
    double ratio{std::pow(p.getMass() / m_max, 1. / 3.)};
    float rad{static_cast<float>(max_radius * ratio)};
    p.setRadius(rad);
  }
}

void manageCollisions(std::vector<Particle> &system) {
  for (std::size_t i{0}; i < system.size(); ++i) {
    for (std::size_t j{i + 1}; j < system.size(); ++j) {
      double r{norm(system[i].getPosition() - system[j].getPosition())};
      double r_min =
          static_cast<double>(system[i].getRadius() + system[j].getRadius());

      if (r <= r_min) {
        double m_new{system[i].getMass() + system[j].getMass()};
        Vec r_new{(system[i].getMass() * system[i].getPosition() +
                   system[j].getMass() * system[j].getPosition()) *
                  (1.0 / m_new)};
        Vec v_new{(system[i].getMass() * system[i].getVelocity() +
                   system[j].getMass() * system[j].getVelocity()) *
                  (1.0 / m_new)};

        float rad_new =
            std::sqrt(system[i].getRadius() * system[i].getRadius() +
                      system[j].getRadius() * system[j].getRadius());

        system[i] = Particle{r_new, v_new, m_new};
        system[i].setRadius(rad_new);

        system.erase(system.begin() + static_cast<long>(j));
        --j;
      }
    }
  }
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

  manageCollisions(system);
}

double kineticEnergy(const std::vector<Particle> &system) {
  return std::accumulate(
      system.begin(), system.end(), 0.0, [](double sum, const Particle &p) {
        return sum + 0.5 * p.getMass() * norm2(p.getVelocity());
      });
}

double potentialEnergy(const std::vector<Particle> &system, const Config &cfg) {
  double U{0.};
  auto const N{system.size()};
  for (std::size_t i{0}; i < N; ++i) {
    for (std::size_t j{i + 1}; j < N; ++j) {
      Vec r{system[j].getPosition() - system[i].getPosition()};
      double inv_d{1.0 / norm(r)};
      U += cfg.G * system[i].getMass() * system[j].getMass() * inv_d;
    }
  }
  return -U;
}

double mechanicalEnergy(const std::vector<Particle> &system,
                        const Config &cfg) {
  return kineticEnergy(system) + potentialEnergy(system, cfg);
}

Vec totalMomentum(const std::vector<Particle> &system) {
  return std::accumulate(
      system.begin(), system.end(), Vec{0., 0.},
      [](const Vec &sum, const Particle &p) { return sum + p.momentum(); });
}

double totalAngularMomentum(const std::vector<Particle> &system) {
  return std::accumulate(
      system.begin(), system.end(), 0.0,
      [](double sum, const Particle &p) { return sum + p.angularMomentum(); });
}

std::vector<Particle> read_from(std::filesystem::path const &p, Config &cfg) {
  std::ifstream is(p);
  if (!is) {
    throw std::filesystem::filesystem_error{
        "read_from", p, std::make_error_code(std::errc::invalid_argument)};
  }

  std::vector<Particle> particles;
  std::string line;

  while (std::getline(is, line)) {
    if (line.empty()) {
      continue;
    }
    if (line[0] == '#') {
      continue;
    }

    if (line[0] == '@') {
      std::istringstream ss(line.substr(1));
      ss >> cfg.G >> cfg.zoom >> cfg.dt >> cfg.speed;
      if (cfg.zoom < 0.) {
        throw std::invalid_argument{" zoom must be at least 0 "};
      }
    } else {
      // riga di particella: mass rx ry vx vy
      std::istringstream ss(line);
      double mass, rx, ry, vx, vy;
      if (ss >> mass >> rx >> ry >> vx >> vy) {
        particles.push_back(Particle{{rx, ry}, {vx, vy}, mass});
      }
    }
  }

  if (particles.empty()) {
    throw std::runtime_error{"no particles found in file"};
  }

  return particles;
}
}  // namespace pf