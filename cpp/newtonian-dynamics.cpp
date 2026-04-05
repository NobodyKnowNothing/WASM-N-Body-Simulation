#include <iostream>
#include <cmath>
#include <span>

typedef struct {
    double x;
    double y;
    double Vx;
    double Vy;
    double Fx;
    double Fy;
    double mass = 1000000000000.0;
    double radius = 1.0;
} particle;

const double G = 0.000000000066743;

inline double particle_distance(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x1-x2, 2) + std::pow(y1-y2, 2));
}

inline double gravity_force(double mass1, double mass2, double pos1, double pos2, double distance) {
    if (distance == 0.0) return 0.0;
    return -1*G*mass1*mass2*(pos1-pos2)/std::pow(distance, 3); // Teach yourself how this was derived
}

inline void reset_forces(std::span<particle*> particles) {
    for (particle *farticle: particles) {
        farticle->Fx = 0;
        farticle->Fy = 0;
    }
}

inline void update_particles(std::span<particle*> particles) {
    for (particle *obj: particles) {
        obj->Vx += (obj->Fx/obj->mass);
        obj->Vy += (obj->Fy/obj->mass);
        obj->x += obj->Vx;
        obj->y += obj->Vy;
    }
}

inline void update_velocities(std::span<particle*> particles, double step_modifier = 1) {
    for (particle *obj: particles) {
        obj->Vx += step_modifier*(obj->Fx/obj->mass);
        obj->Vy += step_modifier*(obj->Fy/obj->mass);
    }
}

inline void update_positions(std::span<particle*> particles) {
    for (particle *obj: particles) {
        obj->x += obj->Vx;
        obj->y += obj->Vy;
    }
}

inline void euler_method(particle *j, particle *k) {
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

    double k1y = gravity_force(i->mass, j->mass, i->y, j->y, r);
    double kV1y = i->Vy + k1y/(i->mass); // Takes force computes velocity

    double r1 = particle_distance(i->x+0.5*(kV1x),i->y+0.5*(kV1y), j->x+0.5*(-kV1x), j->y+0.5*(-kV1y)); // Gotta recompute 
    
    double k2x = gravity_force(i->mass, j->mass, (i->x)+0.5*(kV1x), (j->x)+0.5*(-kV1x), r1);
    double k2y = gravity_force(i->mass, j->mass, (i->y)+0.5*(kV1y), (j->y)+0.5*(-kV1y), r1);
    
    i->Fx += k2x;
    i->Fy += k2y;

    j->Fx -= k2x;
    j->Fy -= k2y;
}

inline void solve(std::span<particle*> particles, int method) {
    switch (method) {
        case 0:
            reset_forces(particles);

            for (size_t i = 0; i < particles.size(); i++) {
                for (size_t j = i + 1; j < particles.size(); j++) {
                    euler_method(particles[i], particles[j]);
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

            for (size_t i = 0; i < particles.size(); i++) {
                for (size_t j = i + 1; j < particles.size(); j++) {
                    euler_method(particles[i], particles[j]);
                }
                std::cout << "X: " << particles[i]->x << ", Y: " << particles[i]->y << std::endl;
            }

            update_velocities(particles, 0.5);
            break;
    }
}

int main() {
    std::cout << "START" << std::endl;

    int method = 2;

    particle* p1 = new particle{};
    particle* p2 = new particle{};
    
    p1->x = -10.0;
    p1->Vy = -2.5;

    p2->x = 10.0;
    p2->Vy = 2.5;

    particle *particles[2] = {p1, p2};

    if (method == 2) {
        for (size_t i = 0; i < particles.size(); i++) {
                for (size_t j = i + 1; j < particles.size(); j++) {
                    euler_method(particles[i], particles[j]);
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