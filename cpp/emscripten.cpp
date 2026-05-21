#include "barnes-hut.hpp"
#include <emscripten.h>
#include <random>

std::vector<particle*> particles;
std::vector<particle*> particlez;

std::vector<double> particle_pos; // pos buffer

double stats_buffer[5];

std::vector<double> deltaF;
int step = 0;
double lyap_sum = 0;
double DT = 0;

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    void add_particle_(double x, double y, double Vx, double Vy, double mass, double radius) {
        add_particle(particles, x, y, Vx, Vy, 0, 0, mass, radius);
    }

    EMSCRIPTEN_KEEPALIVE
    double* get_particle_positions_() {
        particle_pos.clear();
        particle_pos.reserve(particles.size() * 2);

        for (auto *index: particles) {
            particle_pos.push_back(index->x);
            particle_pos.push_back(index->y);
        }

        return particle_pos.data();
    }

    EMSCRIPTEN_KEEPALIVE
    int get_particle_count_() {
        return particle_pos.size();
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
        DT = dt; 
        const double epsilon = 1e-6;
        if (lyap) {
            for (int i = particlez.size(); i < particles.size(); i++) {
                particlez.push_back(new particle(*particles[i]));
            }
            if (step == 0) {
                double shift = epsilon*epsilon/particles.size();

                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<double> dist(0.0, shift);
                std::uniform_int_distribution<int> coin(0, 1);
                for (auto *part: particlez) {
                    double shiftDiff = dist(gen);
                    part->x += coin(gen) ? std::sqrt(shiftDiff) : -std::sqrt(shiftDiff);
                    part->y += coin(gen) ? std::sqrt(shift - shiftDiff) : -std::sqrt(shift - shiftDiff);
                }
            }
            verlet(particlez, dt);
        }
        verlet(particles, dt);
        if (lyap) {
            int size = particles.size();
            double dist_sq = 0;
            if (step % 20 == 0) {
                for (int i = 0; i < size; i++) {
                    double deltay = particles[i]->y - particlez[i]->y;
                    double deltax = particles[i]->x - particlez[i]->x;
                    
                    double deltaVy = particles[i]->Vy - particlez[i]->Vy;
                    double deltaVx = particles[i]->Vx - particlez[i]->Vx;
                    dist_sq += deltax*deltax + deltay*deltay + deltaVx*deltaVx + deltaVy*deltaVy;
                }
                double current_dist = std::sqrt(dist_sq);
                
                lyap_sum += std::log(current_dist / epsilon);
                double scale = epsilon / current_dist; 
                for (int i = 0; i < size; i++) {
                    particlez[i]->x = particles[i]->x + (particlez[i]->x - particles[i]->x)*scale;
                    particlez[i]->y = particles[i]->y + (particlez[i]->y - particles[i]->y)*scale;
                    
                    particlez[i]->Vx = particles[i]->Vx + (particlez[i]->Vx - particles[i]->Vx)*scale;
                    particlez[i]->Vy = particles[i]->Vy + (particlez[i]->Vy - particles[i]->Vy)*scale;
                }
            }
            step += 1;
        }
    }

    EMSCRIPTEN_KEEPALIVE
    double* get_simulation_stats_() {
        stats_buffer[0] = mean_vel(particles);
        stats_buffer[1] = variance_vel(particles);
        stats_buffer[2] = std::sqrt(stats_buffer[1]);
        stats_buffer[3] = (step > 0 && DT > 0) ? (lyap_sum / (step * DT)) : 0;
        stats_buffer[4] = ken_en_sum(particles) + 0.5 * pot_en_sum;

        return stats_buffer;
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
    double get_lyap_expo_() {
        return lyap_sum / (step * DT);
    }
    
    EMSCRIPTEN_KEEPALIVE
    double get_ham_sum_() {
        return ken_en_sum(particles) + 0.5*pot_en_sum;
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
        lyap_sum = 0;
        step = 0;
    }
}