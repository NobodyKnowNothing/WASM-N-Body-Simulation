#include "newtonian-dynamics.hpp"

void solve(std::vector<particle*> particles, int method) {
    switch (method) {
        case 0:
            reset_forces(particles);

            for (size_t i = 0; i < particles.size(); i++) {
                for (size_t j = i + 1; j < particles.size(); j++) {
                    compute_gravity(particles[i], particles[j]);
                }
                std::cout << "X: " << particles[i]->x << ", Y: " << particles[i]->y << std::endl;
            }

            update_particles(particles);
            
            break;
        case 1:
            reset_forces(particles);

            for (size_t i = 0; i < particles.size(); i++) {
                for (size_t j = i + 1; j < particles.size(); j++) {
                    second_order_rk(particles[i], particles[j]);
                }
                std::cout << "X: " << particles[i]->x << ", Y: " << particles[i]->y << std::endl;
            }

            update_particles(particles);
            
            break;
        case 2:
            update_velocities(particles, 0.5);

            update_positions(particles);

            reset_forces(particles);
            for (size_t i = 0; i < particles.size(); i++) {
                for (size_t j = i + 1; j < particles.size(); j++) {
                    compute_gravity(particles[i], particles[j]);
                }
                std::cout << "X: " << particles[i]->x << ", Y: " << particles[i]->y << std::endl;
            }

            update_velocities(particles, 0.5);
            break;
    }
}

int main() {
    std::cout << "START" << std::endl;

    int method = 0;

    particle* p1 = new particle{};
    particle* p2 = new particle{};
    
    p1->x = -10.0;
    p1->Vy = -2.5;

    p2->x = 10.0;
    p2->Vy = 2.5;

    std::vector<particle*> particles = {p1, p2};

    if (method == 2) {
        for (size_t i = 0; i < std::size(particles); i++) {
                for (size_t j = i + 1; j < std::size(particles); j++) {
                    compute_gravity(particles[i], particles[j]);
                }
                std::cout << "X: " << particles[i]->x << ", Y: " << particles[i]->y << std::endl;
        }
    }

    for (int i = 0; i < 15; i++) {
        solve(particles, method);
        std::cout << "---" << std::endl;
    }
    std::cout << "DONE" << std::endl;
}