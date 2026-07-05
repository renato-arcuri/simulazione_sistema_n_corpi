#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "system.hpp"

#include <fstream>

#include "doctest.h"

TEST_CASE("Vector operations 1") {
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

TEST_CASE("computeRadius") {
  pf::Config cfg;

  pf::Particle a({0., 0.}, {3., 0.}, 100000.);
  pf::Particle b({1000., 0.}, {2., 1.}, 2000000.);
  std::vector<pf::Particle> system{a, b};

  pf::computeRadius(system);

  CHECK(system[0].getRadius() == doctest::Approx(36.8403));
  CHECK(system[1].getRadius() == doctest::Approx(100.));
}

TEST_CASE("Acceleration contribution") {
  pf::Config cfg;

  pf::Particle a({0., 0.}, {3., 0.}, 100000.);
  pf::Particle b({1., 0.}, {2., 1.}, 2000000.);

  pf::Vec acc = pf::accelerationContribution(a, b, cfg);

  CHECK(acc.x == doctest::Approx(0.0001334));
  CHECK(acc.y == doctest::Approx(0.));

  pf::Particle a2({0., -1.}, {-2., -2.}, 2345.);
  pf::Particle b2({0., 1.}, {44., 2.}, 3343.);

  pf::Vec acc2a = pf::accelerationContribution(a2, b2, cfg);
  pf::Vec acc2b = pf::accelerationContribution(b2, a2, cfg);

  CHECK(-acc2a.x == doctest::Approx(acc2b.x));
  CHECK(-acc2a.y == doctest::Approx(acc2b.y));
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

TEST_CASE("Free particle motion: computePosition(),setAcceleration()") {
  pf::Config cfg;
  pf::Particle p({0., 0.}, {2000., 3000.}, 1.);
  p.setAcceleration({0., 0.});

  p.computePosition(cfg);

  CHECK(p.getPosition().x == doctest::Approx(2.));
  CHECK(p.getPosition().y == doctest::Approx(3.));
  CHECK(p.getAcceleration().x == 0.0);
  CHECK(p.getAcceleration().y == 0.0);
}

TEST_CASE("Velocity Verlet: step()") {
  pf::Config cfg;
  std::vector<pf::Particle> system;

  pf::Particle p1({10., 5.}, {-1., 3.}, 1.);
  pf::Particle p2({-10., -5.}, {1., -3.}, 1.);
  system.push_back(p1);
  system.push_back(p2);

  pf::step(system, cfg);

  pf::step(system, cfg);

  CHECK(system[0].getAcceleration().x == doctest::Approx(8.44726e-08));
  CHECK(system[0].getAcceleration().y == doctest::Approx(-2.63612e-07));
  CHECK(system[0].getVelocity().x == doctest::Approx(-0.9999873));
  CHECK(system[0].getVelocity().y == doctest::Approx(2.9999604));

  CHECK(system[1].getAcceleration().x == doctest::Approx(-8.44726e-08));
  CHECK(system[1].getAcceleration().y == doctest::Approx(2.63612e-07));
  CHECK(system[1].getVelocity().x == doctest::Approx(0.9999873));
  CHECK(system[1].getVelocity().y == doctest::Approx(-2.9999604));
  // CHECK(p1.getAcceleration().x == doctest::Approx(8.44726e-08));
  // CHECK(p1.getAcceleration().y == doctest::Approx(-2.63612e-07));

  // CHECK(p2.getAcceleration().x == doctest::Approx(-8.44726e-08));
  // CHECK(p2.getAcceleration().y == doctest::Approx(2.63612e-07));

  // CHECK(p1.getVelocity().x == doctest::Approx(-0.9999873));
  // CHECK(p1.getVelocity().y == doctest::Approx(2.9999604));

  // CHECK(p2.getVelocity().x == doctest::Approx(0.9999873));
  // CHECK(p2.getVelocity().y == doctest::Approx(-2.9999604));

  CHECK(system.size() == static_cast<size_t>(2));
}

TEST_CASE("Manage Collisions") {
  SUBCASE("STD") {
    pf::Config cfg;
    cfg.dt = 0.1;
    cfg.G = 1.0;
    std::vector<pf::Particle> system;
    pf::Particle p1({2.0, 0.0}, {-2.0, 0.0}, 1.0);
    pf::Particle p2({-1.0, 0.0}, {1.0, 0.0}, 1.0);

    system.push_back(p1);
    system.push_back(p2);

    system[0].setAcceleration({0., 0.});
    system[1].setAcceleration({0., 0.});
    pf::computeRadius(system);
    pf::Vec mom{pf::totalMomentum(system)};

    for (int i{0}; i != 10; ++i) {
      pf::step(system, cfg);  // manageCollisions is in step()
    }
    // CHECK(system[0].getPosition().x == doctest::Approx(0.0));
    // CHECK(system[0].getVelocity().x == doctest::Approx(-0.5));
    // CHECK(pf::totalMomentum(system).x == doctest::Approx(mom.x));

    // CHECK(system[0].getMass() == doctest::Approx(2.0));

    // CHECK(system.size() == static_cast<size_t>(1));

    CHECK(pf::totalMomentum(system).x == doctest::Approx(mom.x));
    CHECK(system[0].getMass() == doctest::Approx(2.0));
    CHECK(system.size() == static_cast<size_t>(1));
  }

  SUBCASE("OVERLAP") {
    pf::Config cfg;
    std::vector<pf::Particle> system;

    system.push_back(pf::Particle({5., 3.}, {10., 0.}, 2.));
    system.push_back(pf::Particle({5., 3.}, {-10., 0.}, 2.));

    pf::manageCollisions(system);

    CHECK(system[0].getPosition().x == doctest::Approx(5.0));
    CHECK(system[0].getPosition().y == doctest::Approx(3.0));

    CHECK(system[0].getVelocity().x == doctest::Approx(0.0));
    CHECK(system[0].getVelocity().y == doctest::Approx(0.0));

    CHECK(system.size() == static_cast<size_t>(1));
  }
}

TEST_CASE("Energies, Momentum and Angular Momentum") {
  SUBCASE("Kinetic energy") {
    pf::Config cfg;
    std::vector<pf::Particle> sys{{{-2., 3.}, {3., 4.}, 2.}};

    CHECK(pf::kineticEnergy(sys) == doctest::Approx(25.));
  }

  SUBCASE("Potential energy") {
    pf::Config cfg;
    std::vector<pf::Particle> sys{{{0., 3.}, {2., 0.}, 2e11},
                                  {{4., 0.}, {3., 1.}, 5.}};

    CHECK(pf::potentialEnergy(sys, cfg) == doctest::Approx(-13.34));
  }

  SUBCASE("Momentum") {
    pf::Particle p1({0., 0.}, {2., 3.}, 4.);
    pf::Vec mom = p1.momentum();

    CHECK(mom.x == doctest::Approx(8.));
    CHECK(mom.y == doctest::Approx(12.));

    pf::Particle p2({-34., 2.}, {0., -2.}, 2.);
    pf::Vec mom2 = p2.momentum();

    CHECK(mom2.x == doctest::Approx(0.));
    CHECK(mom2.y == doctest::Approx(-4.));
  }

  SUBCASE("Angular momentum") {
    pf::Particle p1({1., 0.}, {0., 2.}, 3.);
    CHECK(p1.angularMomentum() == doctest::Approx(6.));

    pf::Particle p2({1., 0.}, {-1., 0.}, 2.);
    CHECK(p2.angularMomentum() == doctest::Approx(0.));
  }

  SUBCASE("Total Momentum") {
    std::vector<pf::Particle> sys{{{0., 0.}, {1., 0.}, 2.},
                                  {{0., 0.}, {-1., 0.}, 2.}};

    auto P = pf::totalMomentum(sys);

    CHECK(P.x == doctest::Approx(0.));
    CHECK(P.y == doctest::Approx(0.));

    // std::vector<pf::Particle> sys2{{{2., 0.}, {1., 0.}, 2.},
    //                                {{0., 1.}, {-1., 0.}, 2.}};

    // auto P2 = pf::totalMomentum(sys2);
    // CHECK(P2.x == doctest::Approx(0.));
    // CHECK(P2.y == doctest::Approx(0.));
    std::vector<pf::Particle> sys2{{{2., 0.}, {1., 0.}, 2.},
                                   {{0., 1.}, {0., 1.}, 2.}};
    auto P2 = pf::totalMomentum(sys2);
    CHECK(P2.x == doctest::Approx(2.));
    CHECK(P2.y == doctest::Approx(2.));
  }

  SUBCASE("Total Angular momentum") {
    pf::Config cfg;
    cfg.G = 1.;
    cfg.dt = 0.01;
    std::vector<pf::Particle> system;

    system.push_back(pf::Particle{{1.0, 0.0}, {0.0, 1.0}, 2.0});
    system.push_back(pf::Particle{{-1.0, 0.0}, {0.0, -1.0}, 2.0});

    CHECK(pf::totalAngularMomentum(system) == doctest::Approx(4.0));
  }
}

TEST_CASE("Coservation of quantities over time") {
  pf::Config cfg;
  std::vector<pf::Particle> sys{{{-2., 3.}, {3., 4.}, 2.},
                                {{0., 3.}, {2., 0.}, 2},
                                {{4., 0.}, {3., 1.}, 5.}};
  for (auto &p : sys) {
    p.setAcceleration(pf::totalAcceleration(sys, p, cfg));
  }
  pf::step(sys, cfg);

  double E_tot{pf::mechanicalEnergy(sys, cfg)};
  double L{pf::totalAngularMomentum(sys)};
  pf::Vec P{pf::totalMomentum(sys)};

  for (int i{0}; i < 20; ++i) {
    pf::step(sys, cfg);
  }

  SUBCASE("Meccanical energy") {
    CHECK(pf::mechanicalEnergy(sys, cfg) == doctest::Approx(E_tot));
  }

  SUBCASE("Total Momentum") {
    CHECK(pf::totalMomentum(sys).x == doctest::Approx(P.x));
    CHECK(pf::totalMomentum(sys).y == doctest::Approx(P.y));
  }

  SUBCASE("Total Angular momentum") {
    CHECK(pf::totalAngularMomentum(sys) == doctest::Approx(L));
  }
}

TEST_CASE(" Input: read_from() ") {
  SUBCASE(" non-existent file ") {
    pf::Config cfg;
    CHECK_THROWS_AS(pf::read_from("non-existent_file.txt", cfg),
                    std::filesystem::filesystem_error);
  }

  SUBCASE(" empty file ") {
    pf::Config cfg;
    std::string const path{"empty_file.txt"};
    { std::ofstream os(path); }
    CHECK_THROWS_AS(pf::read_from("empty_file.txt", cfg), std::runtime_error);
    std::filesystem::remove(path);
  }

  SUBCASE(" # and @") {
    pf::Config cfg;
    std::string const path{"file.txt"};
    {
      std::ofstream os(path);
      os << "# comment\n";
      os << "1.0 10.0 10.0 -1.0 -1.0\n";
      os << "@ 1.0 1.0 1.0 1\n";
    }
    std::vector<pf::Particle> system = pf::read_from(path, cfg);

    CHECK(system.size() == 1);
    CHECK(system[0].getMass() == doctest::Approx(1.0));

    CHECK(system[0].getPosition().x == doctest::Approx(10.0));
    CHECK(system[0].getPosition().y == doctest::Approx(10.0));

    CHECK(system[0].getVelocity().x == doctest::Approx(-1.0));
    CHECK(system[0].getVelocity().y == doctest::Approx(-1.0));

    CHECK(cfg.G == 1.0);
    CHECK(cfg.zoom == 1.0);
    CHECK(cfg.dt == 1.0);
    CHECK(cfg.speed == 1);

    std::filesystem::remove(path);
  }

  SUBCASE(" without @ ") {
    pf::Config cfg;
    pf::Config cfg_default = cfg;
    std::string const path{"file.txt"};
    {
      std::ofstream os(path);
      os << "# comment\n";
      os << "1.0 10.0 10.0 -1.0 -1.0\n";
    }
    std::vector<pf::Particle> system = pf::read_from(path, cfg);

    CHECK(cfg_default.G == 6.67e-11);
    CHECK(cfg_default.zoom == 0.0);
    CHECK(cfg_default.dt == 0.001);
    CHECK(cfg_default.speed == 17);
    std::filesystem::remove(path);
  }

  SUBCASE("Negative mass") {
    pf::Config cfg;
    std::string const path{"file.txt"};
    {
      std::ofstream os(path);
      os << "-1. 10. 5. -1. 3.\n";
    }

    CHECK_THROWS_AS(pf::read_from("file.txt", cfg), std::invalid_argument);
    std::filesystem::remove(path);
  }

  SUBCASE("Negative zoom") {
    pf::Config cfg;
    std::string const path{"file.txt"};
    {
      std::ofstream os(path);
      os << "@1. -1. 1. 1. \n";
      os << "1. 1. 1. 1. 1.\n";
    }

    CHECK_THROWS_AS(pf::read_from("file.txt", cfg), std::invalid_argument);
    std::filesystem::remove(path);
  }

  SUBCASE(" STD ") {
    pf::Config cfg;
    cfg.dt = 300.0;
    cfg.G = 1.0;
    std::string const path{"file.txt"};
    {
      std::ofstream os(path);
      os << "1. 10. 5. -1. 3.\n";
      os << "1. -10. -5. 1. -3.\n";
    }
    std::vector<pf::Particle> system = pf::read_from(path, cfg);

    pf::step(system, cfg);

    CHECK(system[0].getAcceleration().x == doctest::Approx(8.44726e-08));
    CHECK(system[0].getAcceleration().y == doctest::Approx(-2.63612e-07));
    CHECK(system[0].getVelocity().x == doctest::Approx(-0.9999873));
    CHECK(system[0].getVelocity().y == doctest::Approx(2.9999604));

    CHECK(system[1].getAcceleration().x == doctest::Approx(-8.44726e-08));
    CHECK(system[1].getAcceleration().y == doctest::Approx(2.63612e-07));
    CHECK(system[1].getVelocity().x == doctest::Approx(0.9999873));
    CHECK(system[1].getVelocity().y == doctest::Approx(-2.9999604));

    CHECK(system.size() == static_cast<size_t>(2));

    std::filesystem::remove(path);
  }
}
