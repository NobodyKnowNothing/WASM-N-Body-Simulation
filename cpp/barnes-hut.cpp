#include "barnes-hut.hpp"


int main() {
    std::vector<particle*> particles = {};
    add_particle(particles, 10, 0, 0, 0);
    add_particle(particles, -10, 0, 0, 0);

    double dt = 0.1;
    
    setup_verlet(particles);

    for (int i = 0; i < 200; i++) {
        verlet(particles, dt);
    }

}