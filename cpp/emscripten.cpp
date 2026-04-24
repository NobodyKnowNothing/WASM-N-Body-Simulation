#include "barnes-hut.hpp"
#include <emscripten.h>

std::vector<particle*> particles;



extern "C" {
    EMSCRIPTEN_KEEPALIVE
    void add_particle_(double x, double y, double Vx, double Vy, double mass, double radius) {
        add_particle(particles, x, y, Vx, Vy, 0, 0, mass, radius);
    }

    EMSCRIPTEN_KEEPALIVE
    double particle_get_x_(int index) {
        return particles[index]->x;
    }

    EMSCRIPTEN_KEEPALIVE
    double particle_get_y_(int index) {
        return particles[index]->y;
    }

    EMSCRIPTEN_KEEPALIVE
    void setup_verlet_(int i) {
        if (particles.size() == 0) return;
        qtnode* qtroot = init_qtroot(particles);
        traverse_tree(qtroot, particles[i], compute_gravity);
        delete qtroot;
    }

    EMSCRIPTEN_KEEPALIVE
    void verlet_(double dt, bool lyap = false) {
        if (particles.size() == 0) return;
        if (lyap) {
            const double epsilon = 1e-6;
            std::vector<particle> particlez = particles;
            for (const& particle part: particlez) {
                part->x += epsilon;
                part->y += epsilon;
            }
            verlet(particlez, dt);
        }
        verlet(particles, dt);
        if (lyap) {
            int size = particles.size();
            std::vector<std::vector<double>> deltaF(2, std::vector<double>(size, 0));
            for (int i = 0; i < size; i++) {
                deltaF.push_back({particlez->Fx - particles->Fx, particlez->Fy - particles->Fy});
            }

        }
    }

    EMSCRIPTEN_KEEPALIVE
    double mean_vel_() {
        return mean_vel(particles);
    }

    EMSCRIPTEN_KEEPALIVE
    double variance_vel_() {
        return variance_vel(particles);
    }

    EMSCRIPTEN_KEEPALIVE
    double std_dev_vel_() {
        return std::sqrt(variance_vel(particles));
    }
}