#include <SFML/Graphics.hpp>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>

#include "system.hpp"

int main(int argc, char *argv[]) {
  try {
    pf::Config cfg{};

    // una riga del file deve essere fatto in questo modo(ogni riga è una
    // particella): mass rx ry vx vy
    // la prima riga inizia con @ e contiene i dati i configurazione, se non c'è
    // vengono usati quelli di default
    if (argc < 2) {
      std::cerr << "Uso: " << argv[0] << " <percorso_file_input>\n";
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

    // PARTE GRAFICA

    sf::RenderWindow window(sf::VideoMode(800, 600), "N-bodies Simulation");
    window.setFramerateLimit(60);  // 60 fps

    bool Timeout{false};

    std::vector<sf::Color> palette = {sf::Color::Red,   sf::Color::Yellow,
                                      sf::Color::Green, sf::Color::Cyan,
                                      sf::Color::Blue,  sf::Color::Magenta};

    // std::vector<float> radius;

    double p_min = system[0].getMass();
    for (auto &p : system) {
      if (p.getMass() < p_min) {  // !!!!!!!!!!!!!!!!!!!!! da cambiare
        p_min = p.getMass();
      }
    }

    for (auto &p : system) {
      float rad = static_cast<float>(2 + std::log10(10 * p.getMass() / p_min));
      p.setRadius(rad);
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
        for (int i{0}; i != cfg.speed; ++i) {
          pf::step(system, cfg);
          ++N_step;
        }  // più volte si itera più velocemente avverrà la simulazione
      }

      /*
      if (!system.empty()) {
        double current_p_min = system[0].getMass();
        for (const auto &p : system) {
          if (p.getMass() < current_p_min) {
            current_p_min = p.getMass();
          }
        }
        for (auto &p : system) {
          float rad = static_cast<float>(
              (p.getMass() / current_p_min, 10.0));
          if (rad < 5.f) {
            rad = 5.f;
          }
          p.setRadius(rad);
        }
      }
      */

      // Draw
      window.clear(sf::Color::Black);

      float zoom;
      if (cfg.zoom != 0) {
        zoom = static_cast<float>(cfg.zoom);
      } else {
        double sum_d{0.0};
        for (const auto &p : system) {
          sum_d += pf::norm(p.getPosition());
        }
        double mean_d{sum_d / static_cast<double>(system.size())};
        if (mean_d < 1e-6) {
          mean_d = 1e-6;
        }

        zoom = 50.f / static_cast<float>(mean_d);
      }

      for (std::size_t i = 0; i < system.size(); ++i) {
        const auto &p{system[i]};

        float X{400.f + static_cast<float>(p.getPosition().x) * zoom};
        float Y{300.f - static_cast<float>(p.getPosition().y) * zoom};

        float current_radius = p.getRadius() * zoom;
        if (current_radius < 3.f) {
          current_radius = 3.f;
        }
        sf::CircleShape particle(current_radius);
        particle.setOrigin(current_radius, current_radius);
        particle.setPosition(X, Y);

        particle.setFillColor(palette[i % palette.size()]);
        window.draw(particle);
      }

      if (N_step % cfg.speed == 0) {
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