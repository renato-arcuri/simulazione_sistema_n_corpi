#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "system.hpp"

#include "doctest.h"

TEST_CASE("Vector operations") {
  pf::Vec v1{3., 4.};
  pf::Vec v2{1., 2.};
  pf::Vec v3{0., 0.};
  double a{3.};
  double b{0.};

  auto s = v1 + v2 + v3;
  CHECK(s.x == 4.);
  CHECK(s.y == 6.);

  auto d = v1 - v2 - v3;
  CHECK(d.x == 2.);
  CHECK(d.y == 2.);

  CHECK(pf::norm2(v1) == 25.);
  CHECK(pf::norm(v1) == 5.);
  CHECK(pf::norm(v3) == 0.);

  auto m1 = a * v1;
  auto m2 = v1 * a;
  auto m3 = b * v1;

  CHECK(m1.x == 9.);
  CHECK(m1.y == 12.);
  CHECK(m2.x == 9.);
  CHECK(m2.y == 12.);
  CHECK(m3.x == 0.);
  CHECK(m3.y == 0.);
}

TEST_CASE("Particle members") {
  CHECK_THROWS_AS(pf::Particle({0., 0.}, {0., 0.}, 0.), std::invalid_argument);
  CHECK_THROWS_AS(pf::Particle({0., 0.}, {0., 0.}, -1.), std::invalid_argument);

  pf::Particle p{{1., 2.}, {3., 4.}, 5.};
  pf::Vec a{1., 2.};

  CHECK(p.getPosition().x == 1.);
  CHECK(p.getPosition().y == 2.);
  CHECK(p.getVelocity().x == 3.);
  CHECK(p.getVelocity().y == 4.);
  p.setAcceleration(a);
  CHECK(p.getAcceleration().x == 1.);
  CHECK(p.getAcceleration().y == 2.);
}

TEST_CASE("Particle momentum") {
  pf::Particle p1({0., 0.}, {2., 3.}, 4.);
  auto mom = p1.momentum();

  CHECK(mom.x == 8);
  CHECK(mom.y == 12.);

  pf::Particle p2({-34., 2.}, {0., -2.}, 2.);
  auto mom = p2.momentum();

  CHECK(mom.x == 0.);
  CHECK(mom.y == -4);
}
TEST_CASE("Angular momentum") {
  pf::Particle p1({1., 0.}, {0., 2.}, 3.);
  CHECK(p1.angularMomentum() == doctest::Approx(6.));
  
  pf::Particle p2({1., 0.}, {-1., 0.}, 2.);
  CHECK(p2.angularMomentum() == doctest::approx(0.));
}

/*
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