#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <filesystem>
#include <vector>

namespace pf {

struct Config {
  // valori di default se non configurati diversamente nel file
  double G{6.67e-11};
  double zoom{0.};
  double dt{0.001};
  int speed{17};
};

struct Vec {
  double x;
  double y;
};

Vec operator+(const Vec &v1, const Vec &v2);
Vec operator-(const Vec &v1, const Vec &v2);
Vec operator*(double c, const Vec &v);
Vec operator*(const Vec &v, double c);

double norm2(const Vec &v);
double norm(const Vec &v);

class Particle {
  Vec position_;
  Vec velocity_;
  Vec acceleration_;
  double mass_;
  float k_;

 public:
  Particle(Vec p, Vec v, double m);
  const Vec &getPosition() const;
  const Vec &getVelocity() const;
  const Vec &getAcceleration() const;
  double getMass() const;
  float getRadius() const;

  void setAcceleration(const Vec &a);
  void setRadius(const float &k);

  void computePosition(const Config &cfg);

  void computeVelocity(const Vec &a_old, const Config &cfg);

  Vec momentum() const;

  double angularMomentum() const;
};

Vec accelerationContribution(const Particle &a, const Particle &b,
                             const Config &cfg);
Vec totalAcceleration(const std::vector<Particle> &system, const Particle &p,
                      const Config &cfg);
void computeRadius(std::vector<Particle> &system);
void manageCollisions(std::vector<Particle> &system);
void step(std::vector<Particle> &system, const Config &cfg);
double kineticEnergy(const std::vector<Particle> &system);
double potentialEnergy(const std::vector<Particle> &system, const Config &cfg);
double mechanicalEnergy(const std::vector<Particle> &system, const Config &cfg);
Vec totalMomentum(const std::vector<Particle> &system);
double totalAngularMomentum(const std::vector<Particle> &system);
std::vector<Particle> read_from(std::filesystem::path const &p, Config &cfg);

}  // namespace pf

#endif