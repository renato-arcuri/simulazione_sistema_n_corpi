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
    // la prima riga inizia con @ e contiene i dati i configurazione, se non c'è
    // vengono usati quelli di default
    const std::filesystem::path f{"./input/figure8.txt"};
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

    // PARTE GRAFICA

    sf::RenderWindow window(sf::VideoMode(800, 600), "N-bodies Simulation");
    window.setFramerateLimit(60);  // 60 fps

    bool Timeout{false};

    std::vector<sf::Color> palette = {sf::Color::Red,   sf::Color::Yellow,
                                      sf::Color::Green, sf::Color::Cyan,
                                      sf::Color::Blue,  sf::Color::Magenta};
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
        for (int i{0}; i != 50; ++i) {
          pf::step(system, cfg);
        }  // più volte si itera più velocemente avverrà la simulazione
      }

      // Draw
      window.clear(sf::Color::Black);

      float zoom;
      if (cfg.zoom > 0.) {
        zoom = static_cast<float>(cfg.zoom);
      } else {
        double sum_r{0.0};
        for (const auto &p : system) {
          sum_r += pf::norm(p.getPosition());
        }
        double mean_r{sum_r / static_cast<double>(system.size())};
        if (mean_r < 1e-6) {
          mean_r = 1e-6;
        }
        constexpr float percentage{0.60f};
        constexpr float target_rad{400.f * percentage};
        zoom = target_rad / static_cast<float>(mean_r);
      }

      for (std::size_t i = 0; i < system.size(); ++i) {
        const auto &p{system[i]};

        float X{400.f + static_cast<float>(p.getPosition().x) * zoom};
        float Y{300.f - static_cast<float>(p.getPosition().y) * zoom};

        sf::CircleShape particle(radius[i]);
        particle.setOrigin(radius[i], radius[i]);
        particle.setPosition(X, Y);
        particle.setFillColor(palette[i % palette.size()]);
        window.draw(particle);
      }

      ++N_step;

      if (N_step % 50 == 0) {
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