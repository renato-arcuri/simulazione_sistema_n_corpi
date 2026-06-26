#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "system.hpp"

#include "doctest.h"

/*
TEST_CASE("Vector operations") {
  pf::Vec v1{3., 4.};
  pf::Vec v2{1., 2.};

  auto s = v1 + v2;
  CHECK(s.x == doctest::Approx(4.));
  CHECK(s.y == doctest::Approx(6.));

  auto d = v1 - v2;
  CHECK(d.x == doctest::Approx(2.));
  CHECK(d.y == doctest::Approx(2.));

  CHECK(pf::norm2(v1) == doctest::Approx(25.));
  CHECK(pf::norm(v1) == doctest::Approx(5.));
}
TEST_CASE("Particle mass must be positive") {
  CHECK_THROWS_AS(pf::Particle({0., 0.}, {0., 0.}, 0.), std::invalid_argument);

  CHECK_THROWS_AS(pf::Particle({0., 0.}, {0., 0.}, -1.), std::invalid_argument);
}

TEST_CASE("Particle momentum") {
  pf::Particle p({0., 0.}, {2., 3.}, 4.);

  auto mom = p.momentum();

  CHECK(mom.x == doctest::Approx(8.));
  CHECK(mom.y == doctest::Approx(12.));
}
TEST_CASE("Angular momentum") {
  pf::Particle p({1., 0.}, {0., 2.}, 3.);

  CHECK(p.angularMomentum() == doctest::Approx(6.));
}
TEST_CASE("Acceleration contribution") {
  pf::Config cfg;

  pf::Particle a({0., 0.}, {3., 0.}, 100000.);
  pf::Particle b({1., 0.}, {2., 1.}, 2000000.);

  auto acc = pf::accelerationContribution(a, b, cfg);

  CHECK(acc.x == doctest::Approx(13.34));
  CHECK(acc.y == doctest::Approx(0.));
}
TEST_CASE("Total acceleration") {
  pf::Config cfg;

  std::vector<pf::Particle> sys{{{0., 0.}, {0., 0.}, 1.},
                                {{1., 0.}, {0., 0.}, 2e11},
                                {{0., 1.}, {2., 3.}, 2e11}};

  auto a = pf::totalAcceleration(sys, sys[0], cfg);

  CHECK(a.x == doctest::Approx(13.34));
  CHECK(a.y == doctest::Approx(13.34));
}
TEST_CASE("Kinetic energy") {
  std::vector<pf::Particle> sys{{{-2., 3.}, {3., 4.}, 2.}};

  CHECK(pf::kineticEnergy(sys) == doctest::Approx(25.));
}
TEST_CASE("Potential energy") {
  pf::Config cfg;

  std::vector<pf::Particle> sys{{{0., 3.}, {2., 0.}, 2e11},
                                {{4., 0.}, {3., 1.}, 5.}};

  CHECK(pf::potentialEnergy(sys, cfg) == doctest::Approx(-13.34));
}
TEST_CASE("Total momentum") {
  std::vector<pf::Particle> sys{{{0., 0.}, {1., 0.}, 2.},
                                {{0., 0.}, {-1., 0.}, 2.}};

  auto P = pf::totalMomentum(sys);

  CHECK(P.x == doctest::Approx(0.));
  CHECK(P.y == doctest::Approx(0.));

  std::vector<pf::Particle> sys2{{{2., 0.}, {1., 0.}, 2.},
                                 {{0., 1.}, {-1., 0.}, 2.}};

  auto P2 = pf::totalMomentum(sys2);
  CHECK(P2.x == doctest::Approx(2.));
  CHECK(P2.y == doctest::Approx(0.));
}
TEST_CASE("Free particle motion") {
  pf::Config cfg;
  pf::Particle p({0., 0.}, {2000., 3000.}, 1.);
  p.setAcceleration({0., 0.});

  p.computePosition(cfg);

  CHECK(p.getPosition().x == doctest::Approx(2.));
  CHECK(p.getPosition().y == doctest::Approx(3.));
}
*/