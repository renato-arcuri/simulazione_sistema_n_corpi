#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

namespace pf {

struct Config {
  double G{6.67e-11};
  double dt{0.001};
  double softening{1.e-12};
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
Vec operator*(double c, const Vec &v) { return Vec{v.x * c, v.y * c}; }
Vec operator*(const Vec &v, double c) { return Vec{v.x * c, v.y * c}; }

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
    position_ =
        position_ + velocity_ * cfg.dt + acceleration_ * 0.5 * cfg.dt * cfg.dt;
  }

  void computeVelocity(const Vec &a_old, const Config &cfg) {
    velocity_ = velocity_ + 0.5 * cfg.dt * (a_old + acceleration_);
  }

  Vec momentum() const { return velocity_ * mass_; }

  double angularMomentum() const {
    return mass_ * ((position_.x * velocity_.y) - (position_.y * velocity_.x));
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

}  // namespace pf

int main() {
  pf::Config cfg{};

  ////CONFIG INIZIALE (con questa non va...)
  // const double m1{2.0e24};
  // const double m2{1.0e17};
  // const double m3{1.0e13};

  // const pf::Vec r1{0.0, 0.0};
  // const pf::Vec r2{10.0, 0.0};
  // const pf::Vec r3{10.7, 0.0};

  // const pf::Vec v1{0.0, 0.0};
  // const pf::Vec v2{0.0, 115.498};
  // const pf::Vec v3{0.0, 115.498 + 2.58};

  // std::vector<pf::Particle> system{{r1, v1, m1}, {r2, v2, m2}, {r3, v3, m3}};
  /*
    pf::Particle{
        {-0.97000436, 0.24308753e10}, {0.4662036850, 0.4323657300}, m},
    pf::Particle{{0.97000436, -0.24308753}, {0.4662036850,
    0.4323657300}, m}, pf::Particle{{0., 0.}, {-0.93240737,
    -0.86473146}, m},};
  */

  // SIMULAZIONE CON SISTEMA SOLARE RIDOTTO

  std::vector<pf::Particle> system;
  // --- SOLE ---
  const double m1{6.56e11};
  const pf::Vec r1{0.0, 0.0};
  const pf::Vec v1{0.0, 0.0};

  // --- MERCURIO ---
  const double m2{1.08e5};
  const pf::Vec r2{2.7, 0.0};
  const pf::Vec v2{0.0, 4.02};

  // --- VENERE ---
  const double m3{1.61e6};
  const pf::Vec r3{5.1, 0.0};
  const pf::Vec v3{0.0, 2.93};

  // --- TERRA ---
  const double m4{1.97e6};
  const pf::Vec r4{7.0, 0.0};
  const pf::Vec v4{0.0, 2.50};

  // --- MARTE ---
  const double m5{2.11e5};
  const pf::Vec r5{10.7, 0.0};
  const pf::Vec v5{0.0, 2.02};

  // --- GIOVE ---
  const double m6{6.26e8};
  const pf::Vec r6{-15.0, 0.0};
  const pf::Vec v6{0.0, -1.71};

  // --- SATURNO ---
  const double m7{1.87e8};
  const pf::Vec r7{27.5, 0.0};
  const pf::Vec v7{0.0, 1.26};

  // --- URANO ---
  const double m8{2.86e7};
  const pf::Vec r8{-55.3, 0.0};
  const pf::Vec v8{0.0, -0.89};

  // --- NETTUNO ---
  const double m9{3.37e7};
  const pf::Vec r9{86.7, 0.0};
  const pf::Vec v9{0.0, 0.71};

  system.push_back(pf::Particle{r1, v1, m1});
  system.push_back(pf::Particle{r2, v2, m2});
  system.push_back(pf::Particle{r3, v3, m3});
  system.push_back(pf::Particle{r4, v4, m4});
  system.push_back(pf::Particle{r5, v5, m5});
  system.push_back(pf::Particle{r6, v6, m6});
  system.push_back(pf::Particle{r7, v7, m7});
  system.push_back(pf::Particle{r8, v8, m8});
  system.push_back(pf::Particle{r9, v9, m9});

  //--------------------

  for (auto &p : system) {
    p.setAcceleration(pf::totalAcceleration(system, p, cfg));
  };

  constexpr int N_STEPS{1000000};
  constexpr int PRINT_EVERY{50000};

  std::cout << "\n";
  std::cout << "-------------------------------------------------------------"
               "------------------------------------\n";

  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout.precision(4);

  std::cout.width(8);
  std::cout << "Step";
  std::cout.width(15);
  std::cout << "E_tot";
  std::cout.width(15);
  std::cout << "Kinetic";
  std::cout.width(15);
  std::cout << "Potential";
  std::cout.width(15);
  std::cout << "P";
  std::cout.width(15);
  std::cout << "L_z" << "\n";

  std::cout << "-------------------------------------------------------------"
               "------------------------------------\n";

  for (int s{0}; s < N_STEPS; ++s) {
    pf::step(system, cfg);

    if (s % PRINT_EVERY == 0) {
      double E_tot{pf::mechanicalEnergy(system, cfg)};
      double K{pf::kineticEnergy(system)};
      double U{pf::potentialEnergy(system, cfg)};
      pf::Vec P{pf::totalMomentum(system)};
      double L{pf::totalAngularMomentum(system)};

      std::cout.width(8);
      std::cout << s;
      std::cout.width(15);
      std::cout << E_tot;
      std::cout.width(15);
      std::cout << K;
      std::cout.width(15);
      std::cout << U;
      std::cout.width(15);
      std::cout << pf::norm(P);
      std::cout.width(15);
      std::cout << L << "\n";
    }
  }

  std::cout << "-------------------------------------------------------------"
               "------------------------------------\n";

  //--------------------

  // PARTE GRAFICA

  sf::RenderWindow window(sf::VideoMode(800, 600), "N-bodies Simulation");
  window.setFramerateLimit(60);

  bool Timeout{false};

  std::vector<sf::Color> palette = {sf::Color::Red, sf::Color::Yellow,
                                    sf::Color::Green, sf::Color::Cyan};

  while (window.isOpen()) {
    sf::Event event;

    // Events
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed) {
        window.close();
      }

      if (event.type == sf::Event::KeyPressed &&
          event.key.code == sf::Keyboard::Space) {
        Timeout = !Timeout;
      }
    }

    // Step
    if (!Timeout) {
      for (int i{0}; i != 15; ++i) {
        pf::step(system, cfg);
      }  // serve farlo 15 volte perché il dt=1 ms mentre il frame
    }

    // Draw
    window.clear(sf::Color::Black);

    double sum_r{0.0};
    for (const auto &p : system) {
      sum_r += pf::norm(p.getPosition());
    }
    double mean_r{sum_r / system.size()};

    float percentage{0.60f};
    float target_radius{400.f * percentage};
    float zoom{target_radius / static_cast<float>(mean_r)};

    for (auto i = 0; i < system.size(); ++i) {
      const auto &p{system[i]};

      float radius{static_cast<float>(std::pow(p.getMass(), 1.0 / 5.0)) *
                   0.05f};
      if (radius < 3.f) radius = 3.f;
      if (radius > 20.f) radius = 20.f;

      sf::CircleShape particle(radius);
      particle.setOrigin(radius, radius);

      float X{400.f + static_cast<float>(p.getPosition().x) * zoom};
      float Y{300.f - static_cast<float>(p.getPosition().y) * zoom};
      particle.setPosition(X, Y);

      particle.setFillColor(palette[i % palette.size()]);

      window.draw(particle);
    }

    window.display();
  }
  return 0;
};  // main