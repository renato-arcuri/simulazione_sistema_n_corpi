#include <SFML/Graphics.hpp>
#include <cmath>
#include <exception>
#include <iostream>

#include "system.hpp"

int main() {
  try {
    pf::Config cfg{};

    // una riga del file deve essere fatto in questo modo(ogni riga è una
    // particella): mass rx ry vx vy
    const std::filesystem::path f{"..."};
    std::vector<pf::Particle> system{read_from(f)};

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
    /*
      pf::Particle{
          {-0.97000436, 0.24308753e10}, {0.4662036850, 0.4323657300}, m},
      pf::Particle{{0.97000436, -0.24308753}, {0.4662036850,
      0.4323657300}, m}, pf::Particle{{0., 0.}, {-0.93240737,
      -0.86473146}, m},};
    */

    // SIMULAZIONE CON SISTEMA SOLARE RIDOTTO

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
      if (rad < 2.f) {
        rad = 2.f;
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
        for (int i{0}; i != 15; ++i) {
          pf::step(system, cfg);
        }  // serve farlo 15 volte perché il dt=1 ms mentre il frame 15ms
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

      if (N_step % 500 == 0) {
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