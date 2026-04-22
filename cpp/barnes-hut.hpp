#ifndef BARNES_HUT_H
#define BARNES_HUT_H


#include <algorithm>
#include "newtonian-dynamics.hpp"

struct qtnode {
    double CoMx = 0;
    double CoMy = 0;
    double cx = 0;
    double cy = 0;
    double totalMass = 0;
    double length = 0;
    std::vector<particle*> particles;
    qtnode* NW = nullptr;
    qtnode* NE = nullptr;
    qtnode* SW = nullptr;
    qtnode* SE = nullptr;
    ~qtnode() {
        delete NE;
        delete NW;
        delete SE;
        delete SW;
    }
};

const double theta = 1;

inline void traverse_tree(qtnode* current, particle* i, void (*method)(particle*, particle*)) {
    if (current == nullptr) return;
    bool part_size = (current->particles.size() == 1);
    particle* curr_part;
    if (part_size) {
        curr_part = current->particles[0];
        if (curr_part == i) return;
    }
    double dist = std::max(particle_distance(i->x, i->y, current->CoMx, current->CoMy), 0.0001);
    double len = current->length;
    if (!part_size && len/dist >= theta) {
        traverse_tree(current->NE, i, method);
        traverse_tree(current->NW, i, method);
        traverse_tree(current->SW, i, method);
        traverse_tree(current->SE, i, method);
        return;
    }
    particle j{.x = current->CoMx, .y = current->CoMy, .mass = current->totalMass};
    method(i, &j);
    if (part_size && dist < (i->radius + curr_part->radius)) {
        compute_collision(i, curr_part);
    }
}


inline void qt_aux(qtnode* header, int maxSize) {
    std::vector<particle*> subsets[4];
    for (const auto particle: header->particles) {
        if (particle->x > header->cx) {
            if (particle->y > header->cy) {
                subsets[0].push_back(particle);
            } else {
                subsets[3].push_back(particle);
            }
        } else {
            if (particle->y > header->cy) {
                subsets[1].push_back(particle);
            } else {
                subsets[2].push_back(particle);
            }
        }
    }
    int s0 = subsets[0].size(); 
    if (s0 > 0) {
        header->NE = new qtnode{};
        header->NE->particles = subsets[0];
        header->NE->length = header->length/2;
        header->NE->cy = header->cy + header->NE->length/2;
        header->NE->cx = header->cx + header->NE->length/2;
        if (s0 > maxSize) {
            qt_aux(header->NE, maxSize);
            header->CoMx += header->NE->CoMx * header->NE->totalMass;
            header->CoMy += header->NE->CoMy * header->NE->totalMass;
            header->totalMass += header->NE->totalMass;
        } else {
            for (const auto *particle : subsets[0]) {
                header->NE->CoMx += particle->x * particle->mass;
                header->NE->CoMy += particle->y * particle->mass;
                header->NE->totalMass += particle->mass;
            }
            
            header->CoMx += header->NE->CoMx;
            header->CoMy += header->NE->CoMy;

            header->NE->CoMx /= header->NE->totalMass;
            header->NE->CoMy /= header->NE->totalMass;

            header->totalMass += header->NE->totalMass;
        }
    }

    int s1 = subsets[1].size(); 
    if (s1 > 0) {
        header->NW= new qtnode{};
        header->NW->particles = subsets[1];
        header->NW->length = header->length/2;
        header->NW->cy = header->cy + header->NW->length/2;
        header->NW->cx = header->cx - header->NW->length/2;
        if (s1 > maxSize) {
            qt_aux(header->NW, maxSize);
            header->CoMx += header->NW->CoMx * header->NW->totalMass;
            header->CoMy += header->NW->CoMy * header->NW->totalMass;
            header->totalMass += header->NW->totalMass;
        } else {
            for (const auto *particle : subsets[1]) {
                header->NW->CoMx += particle->x * particle->mass;
                header->NW->CoMy += particle->y * particle->mass;
                header->NW->totalMass += particle->mass;
            }
            
            header->CoMx += header->NW->CoMx;
            header->CoMy += header->NW->CoMy;

            header->NW->CoMx /= header->NW->totalMass;
            header->NW->CoMy /= header->NW->totalMass;

            header->totalMass += header->NW->totalMass;
        }
    }
    int s2 = subsets[2].size(); 
    if (s2 > 0) {
        header->SW = new qtnode{};
        header->SW->particles = subsets[2];
        header->SW->length = header->length/2;
        header->SW->cy = header->cy - header->SW->length/2;
        header->SW->cx = header->cx - header->SW->length/2;
        if (s2 > maxSize) {
            qt_aux(header->SW, maxSize);
            header->CoMx += header->SW->CoMx * header->SW->totalMass;
            header->CoMy += header->SW->CoMy * header->SW->totalMass;
            header->totalMass += header->SW->totalMass;
        } else {
            for (const auto *particle : subsets[2]) {
                header->SW->CoMx += particle->x * particle->mass;
                header->SW->CoMy += particle->y * particle->mass;
                header->SW->totalMass += particle->mass;
            }
            
            header->CoMx += header->SW->CoMx;
            header->CoMy += header->SW->CoMy;

            header->SW->CoMx /= header->SW->totalMass;
            header->SW->CoMy /= header->SW->totalMass;

            header->totalMass += header->SW->totalMass;
        }
    }
    int s3 = subsets[3].size(); 
    if (s3 > 0) {
        header->SE = new qtnode{};
        header->SE->particles = subsets[3];
        header->SE->length = header->length/2;
        header->SE->cy = header->cy - header->SE->length/2;
        header->SE->cx = header->cx + header->SE->length/2;
        if (s3 > maxSize) {
            qt_aux(header->SE, maxSize);
            header->CoMx += header->SE->CoMx * header->SE->totalMass;
            header->CoMy += header->SE->CoMy * header->SE->totalMass;
            header->totalMass += header->SE->totalMass;
        } else {
            for (const auto *particle : subsets[3]) {
                header->SE->CoMx += particle->x * particle->mass;
                header->SE->CoMy += particle->y * particle->mass;
                header->SE->totalMass += particle->mass;
            }
            
            header->CoMx += header->SE->CoMx;
            header->CoMy += header->SE->CoMy;

            header->SE->CoMx /= header->SE->totalMass;
            header->SE->CoMy /= header->SE->totalMass;

            header->totalMass += header->SE->totalMass;
        }
    }
    header->CoMx /= header->totalMass;
    header->CoMy /= header->totalMass;
}

inline qtnode* init_qtroot(std::vector<particle*> particles, int maxSize = 1) {
    qtnode* qtroot = new qtnode{};

    auto [min_x, max_x] = std::minmax_element(particles.begin(), particles.end(), [](const particle* a, const particle* b) {return a->x < b->x;});
    
    auto [min_y, max_y] = std::minmax_element(particles.begin(), particles.end(), [](const particle* a, const particle* b) {return a->y < b->y;});

    double x1 = (*min_x)->x;
    double x2 = (*max_x)->x;

    double y1 = (*min_y)->y;
    double y2 = (*max_y)->y;

    double xSpan = (x2 - x1);
    double ySpan = y2 - y1;
    double span = (std::abs(xSpan) > std::abs(ySpan)) ? xSpan : ySpan;
    
    qtroot->length = std::abs(span);
    qtroot->cx = x1 + xSpan/2;
    qtroot->cy = y1 + ySpan/2;
    qtroot->particles = particles;

    qtnode* header = qtroot;
    qt_aux(header, maxSize);
    return qtroot;
}

inline void verlet(std::vector<particle*> particles, double dt = 1.0) {
    update_velocities(particles, 0.5*dt);

    update_positions(particles, dt);

    qtnode* qtroot = init_qtroot(particles);

    reset_forces(particles);
    for (size_t i = 0; i < particles.size(); i++) {
        traverse_tree(qtroot, particles[i], compute_gravity);
    }
    
    update_velocities(particles, 0.5*dt);

    delete qtroot;
}

inline void setup_verlet(std::vector<particle*> particles){
    qtnode* qtroot = init_qtroot(particles);

    for (size_t i = 0; i < particles.size(); i++) {
        traverse_tree(qtroot, particles[i], compute_gravity);
    }

    delete qtroot;
}

inline double mean_vel(std::vector<particle*> particles) {
    double totalVx = 0;
    double totalVy = 0;

    for (const particle *part: particles) {
        totalVx += part->Vx;
        totalVy += part->Vy;
    }

    return std::hypot(totalVx, totalVy)/particles.size();
};

inline double variance_vel(std::vector<particle*> particles) {
    double mean = mean_vel(particles);
    double mean_diff = 0;
    for (const particle *part: particles) {
        mean_diff += std::hypot(part->Vx, part->Vy) - mean;
    }
    return std::pow(mean_diff, 2) / particles.size();
}

#endif