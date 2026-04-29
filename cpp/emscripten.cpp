#include "barnes-hut.hpp"
#include <emscripten.h>
#include <random>

std::vector<particle*> particles;
std::vector<particle*> particlez;
std::vector<double> deltaF;
int frame = 0;
int lastFrame = 0;
double lyap_sum = 0;


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
    void verlet_(double dt, bool lyap = false, int currFrame = 0) {
        if (particles.size() == 0) return;
        const double epsilon = 1e-6;
        if (lyap) {
            if (lastFrame == 0) {
                double shift = epsilon/particles.size();

                std::random_device rd;
                std::mt19937 gen(rd());
                for (auto *part: particlez) {
                    std::uniform_real_distribution<double> dist(0.0, shift);
                    double shiftDiff = dist(gen);
                    part->x += shiftDiff;
                    part->y += shift - shiftDiff;
                }
                lastFrame = currFrame;
            }
            verlet(particlez, dt);
        }
        verlet(particles, dt);
        if (lyap) {
            int size = particles.size();
            double dist_sq = 0;
            if (currFrame % 15 == 0) {
                for (int i = 0; i < size; i++) {
                    double deltaFy = particlez[i]->Fy - particles[i]->Fy;
                    double deltaFx = particlez[i]->Fx - particles[i]->Fx;
                    dist_sq += deltaFx*deltaFx + deltaFy*deltaFy;
                }
                double current_dist = std::sqrt(dist_sq);
                
                lyap_sum += std::log(current_dist / epsilon);

                particlez = particles;
                double scale = epsilon / current_dist; 
                for (int i = 0; i < size; i++) {
                    particlez[i]->x = particles[i]->x + (particles[i]->x - particlez[i]->x)*scale;
                    particlez[i]->y = particles[i]->y + (particles[i]->y - particlez[i]->y)*scale;
                }
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

    EMSCRIPTEN_KEEPALIVE
    double get_lyap_sum_() {
        return lyap_sum;
    }

    EMSCRIPTEN_KEEPALIVE
    void reset_() {
        for (particle* part: particles) {
            delete part;
        }
        for (particle* part: particlez) {
            delete part;
        }
        particles = {};
        particlez = {};
        deltaF = {};
        frame = 0;
        lastFrame = 0;
        lyap_sum = 0;
    }
}