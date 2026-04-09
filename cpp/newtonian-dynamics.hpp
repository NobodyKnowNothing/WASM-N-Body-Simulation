#ifndef NEWTONIAN_DYNAMICS_H
#define NEWTONIAN_DYNAMICS_H

#include <iostream>
#include <cmath>
#include <vector>

// Structure definition
struct particle {
    double x;
    double y;
    double Vx;
    double Vy;
    double Fx;
    double Fy;
    double mass = 1000000000000.0;
    double radius = 1.0;
};

// Global Constants
extern const double G = 0.000000000066743;

inline double particle_distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x1-x2, 2) + std::pow(y1-y2, 2));
}

inline double gravity_force(double mass1, double mass2, double pos1, double pos2, double distance) {
    if (distance == 0.0) return 0.0;
    return -1*G*mass1*mass2*(pos1-pos2)/std::pow(distance, 3); // Teach yourself how this was derived
}

inline void reset_forces(std::vector<particle*> particles) {
    for (particle *farticle: particles) {
        farticle->Fx = 0;
        farticle->Fy = 0;
    }
}

inline void update_particles(std::vector<particle*> particles) {
    for (particle *obj: particles) {
        obj->x += obj->Vx;
        obj->y += obj->Vy;
        obj->Vx += (obj->Fx/obj->mass);
        obj->Vy += (obj->Fy/obj->mass);
    }
}

inline void update_velocities(std::vector<particle*> particles, double step_modifier = 1.0) {
    for (particle *obj: particles) {
        obj->Vx += step_modifier*(obj->Fx/obj->mass);
        obj->Vy += step_modifier*(obj->Fy/obj->mass);
    }
}

inline void update_positions(std::vector<particle*> particles, double step_modifier = 1.0) {
    for (particle *obj: particles) {
        obj->x += step_modifier*obj->Vx;
        obj->y += step_modifier*obj->Vy;
    }
}

inline void compute_gravity(particle *j, particle *k) {
    double r = particle_distance(j->x, j->y, k->x, k->y);

    double Fx = gravity_force(j->mass, k->mass, j->x, k->x, r);
    double Fy = gravity_force(j->mass, k->mass, j->y, k->y, r);

    j->Fx += Fx;
    j->Fy += Fy;

    k->Fx -= Fx;
    k->Fy -= Fy;

}

// May be crude and unusual and not what I want it to be, may have to update entire system by half a timestep
inline void second_order_rk(particle *i, particle *j) {
    double r = particle_distance(i->x, i->y, j->x, j->y); 

    double k1x = gravity_force(i->mass, j->mass, i->x, j->x, r);
    double kV1x = i->Vx + k1x/(i->mass); // Takes force computes velocity
    double jV1x = j->Vx - k1x/(j->mass);

    double k1y = gravity_force(i->mass, j->mass, i->y, j->y, r);
    double kV1y = i->Vy + k1y/(i->mass); // Takes force computes velocity
    double jV1y = j->Vy - k1y/(j->mass);



    double r1 = particle_distance(i->x+kV1x,i->y+kV1y, j->x+jV1x, j->y+jV1y); // Gotta recompute 
    
    double k2x = gravity_force(i->mass, j->mass, (i->x)+kV1x, (j->x)+jV1x, r1);
    double k2y = gravity_force(i->mass, j->mass, (i->y)+kV1y, (j->y)+jV1y, r1);
    
    i->Fx += (k1x+k2x)/2;
    i->Fy += (k1y+k2y)/2;
    j->Fx -= (k1x+k2x)/2;
    j->Fy -= (k1y+k2y)/2;
}

inline double dist_mag(particle* i, particle* j) {
    double x = i->x - j->x;
    double y = i->y - j->y;
    return std::hypot(x, y);
}


inline std::vector<double> normal(particle* i, particle* j) {
    double x = i->x - j->x;
    double y = i->y - j->y;
    double mag = dist_mag(i, j);
    return {x/mag, y/mag};
    
}

inline double impulse_scalar(particle* i, particle* j, std::vector<double> n, double e = 1.0) {
    double j1 = -(1 + e)/(1/i->mass + 1/j->mass);
    double x = i->Vx - j->Vx;
    double y = i->Vy - j->Vy;
    return j1*(x*n[0] + y*n[1]);
}

inline void compute_collision(particle* i, particle* k) {
    std::vector<double> n = normal(i, k);
    if (((i->Vx-k->Vx)*n[0] + (i->Vy-k->Vy)*n[1]) >= 0) return;
    double j = impulse_scalar(i, k, n);
    std::vector<double> n1 = {n[0]*j, n[1]*j};

    i->Vx += n1[0]/i->mass;
    i->Vy += n1[1]/i->mass;

    k->Vx -= n1[0]/k->mass;
    k->Vy -= n1[1]/k->mass;

    double overlap = (i->radius + k->radius) - dist_mag(i, k);

    i->x += n[0]*overlap/2;
    i->y += n[1]*overlap/2;

    k->x -= n[0]*overlap/2;
    k->y -= n[1]*overlap/2;
}


inline void add_particle(std::vector<particle*>& particles, double x, double y, double Vx = 0, double Vy = 0, double Fx = 0, double Fy = 0, double mass = 1000000000000.0, double radius = 2.0) {
    particles.push_back(new particle{.x = x, .y = y, .Vx = Vx, .Vy = Vy, .Fx = Fx, .Fy = Fy, .mass = mass, .radius = radius});
}

#endif