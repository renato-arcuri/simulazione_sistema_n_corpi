#include <SFML/Graphics.hpp>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>

#include "system.hpp"

int main() {
  try {
    pf::Config cfg{};

    // una riga del file deve essere fatto in questo modo(ogni riga è una
    // particella): mass rx ry vx vy
    const std::filesystem::path f{"./input/solar_system.txt"};
    std::vector<pf::Particle> system{pf::read_from(f)};

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

    // std::vector<pf::Particle> system{{r1, v1, m1}, {r2, v2, m2}, {r3, v3,
    // m3}};

    for (auto &p : system) {
      p.setAcceleration(pf::totalAcceleration(system, p, cfg));
    };

    int N_step{0};
    double E_tot{pf::mechanicalEnergy(system, cfg)};
    double K{pf::kineticEnergy(system)};
    double U{pf::potentialEnergy(system, cfg)};
    pf::Vec P{pf::totalMomentum(system)};
    double L{pf::totalAngularMomentum(system)};

    //--------------------

    // PARTE GRAFICA

    sf::RenderWindow window(sf::VideoMode(800, 600), "N-bodies Simulation");
    window.setFramerateLimit(60);  // 60 fps

    bool Timeout{false};

    std::vector<sf::Color> palette = {sf::Color::Red, sf::Color::Yellow,
                                      sf::Color::Green, sf::Color::Cyan};
    std::vector<float> radius;

    double p_min = system[0].getMass();
    for (auto &p : system) {
      if (p.getMass() < p_min) {
        p_min = p.getMass();
      }
    }

    for (const auto &p : system) {
      float rad = static_cast<float>(std::pow(p.getMass() / p_min, 1.0 / 6.0));
      if (rad < 5.f) {
        rad = 5.f;
      }
      radius.push_back(rad);
    }

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
      if (Timeout == false) {
        for (int i{0}; i != 60; ++i) {
          pf::step(system, cfg);
        }  // più volte si itera più velocemente avverrà la simulazione
      }

      // Draw
      window.clear(sf::Color::Black);

      double sum_r{0.0};
      for (const auto &p : system) {
        sum_r += pf::norm(p.getPosition());
      }
      double mean_r{sum_r / static_cast<double>(system.size())};

      float percentage{0.60f};
      float target_rad{400.f * percentage};
      float zoom{target_rad / static_cast<float>(mean_r)};

      // constexpr float zoom{200.f};

      for (std::size_t i = 0; i < system.size(); ++i) {
        const auto &p{system[i]};

        sf::CircleShape particle(radius[i]);
        particle.setOrigin(radius[i], radius[i]);

        float X{400.f + static_cast<float>(p.getPosition().x) * zoom};
        float Y{300.f - static_cast<float>(p.getPosition().y) * zoom};
        particle.setPosition(X, Y);

        particle.setFillColor(palette[i % palette.size()]);

        window.draw(particle);
      }

      ++N_step;

      if (N_step % 50 == 0) {
        std::cout << "\n";

        std::cout.width(8);
        std::cout << N_step;
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

        // std::cout
        //     <<
        //     "-------------------------------------------------------------"
        //        "------------------------------------\n";
      }
      window.display();
    }  // game loop
  } catch (std::filesystem::filesystem_error const &e) {
    std::cerr << e.path1();
    return EXIT_FAILURE;
  } catch (std::exception const &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "unknown exception\n";
    return EXIT_FAILURE;
  }

};  // main