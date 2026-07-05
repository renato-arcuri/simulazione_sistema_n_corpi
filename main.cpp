#include <SFML/Graphics.hpp>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>

#include "system.hpp"

int main(int argc, char *argv[]) {
  try {
    // Input
    pf::Config cfg{};
    if (argc < 2) {
      std::cerr << "Use: " << argv[0] << " <path_file_input>\n";
      return EXIT_FAILURE;
    }
    const std::filesystem::path f{argv[1]};
    auto system{pf::read_from(f, cfg)};

    for (auto &p : system) {
      p.setAcceleration(pf::totalAcceleration(system, p, cfg));
    };

    int N_step{0};
    double E_tot{0.};
    double K{0.};
    double U{0.};
    pf::Vec P{0., 0.};
    double L{0.};

    // Graphics
    sf::RenderWindow window(sf::VideoMode(800, 600), "N-bodies Simulation");
    window.setFramerateLimit(60);  // 60 FPS

    bool Timeout{false};
    std::vector<sf::Color> palette = {sf::Color::Red,   sf::Color::Yellow,
                                      sf::Color::Green, sf::Color::Cyan,
                                      sf::Color::Blue,  sf::Color::Magenta};

    pf::computeRadius(system);

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

    float zoom{0.};

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
        for (int i{0}; i != cfg.speed; ++i) {
          pf::step(system, cfg);
          ++N_step;
        }
      }

      // Draw
      window.clear(sf::Color::Black);

      if (cfg.zoom > 0.) {
        zoom = static_cast<float>(cfg.zoom);
      } else {
        double max_d{1.e-6};
        if (!system.empty()) {
          auto it = std::max_element(
              system.begin(), system.end(),
              [](const pf::Particle &a, const pf::Particle &b) {
                return pf::norm(a.getPosition()) < pf::norm(b.getPosition());
              });
          max_d = std::max(max_d, pf::norm(it->getPosition()));
        }
        if (max_d < 1e-6) {
          max_d = 1e-6;
        }

        float target_zoom{350.f / static_cast<float>(max_d)};
        zoom += 0.05f * (target_zoom - zoom);
      }

      for (std::size_t i = 0; i < system.size(); ++i) {
        const auto &p{system[i]};

        float X{400.f + static_cast<float>(p.getPosition().x) * zoom};
        float Y{300.f - static_cast<float>(p.getPosition().y) * zoom};

        float current_radius{p.getRadius() * zoom};
        if (current_radius < 3.f) {
          current_radius = 3.f;
        }

        sf::CircleShape particle(current_radius);
        particle.setOrigin(current_radius, current_radius);
        particle.setPosition(X, Y);

        particle.setFillColor(palette[i % palette.size()]);
        window.draw(particle);
      }

      // Terminal
      if (N_step % (cfg.speed * 60) == 0) {  // print every second
        std::cout << "\n";

        E_tot = pf::mechanicalEnergy(system, cfg);
        K = pf::kineticEnergy(system);
        U = pf::potentialEnergy(system, cfg);
        P = pf::totalMomentum(system);
        L = pf::totalAngularMomentum(system);

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
      }
      window.display();
    }  // game loop

  } catch (std::filesystem::filesystem_error const &e) {
    std::cerr << e.what();
    return EXIT_FAILURE;
  } catch (std::exception const &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "unknown exception\n";
    return EXIT_FAILURE;
  }
}  // main