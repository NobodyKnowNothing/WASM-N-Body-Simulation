#include <iostream>
#include <cmath>

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

double inline particle_distance(particle *obj, particle *obj1) {
    return std::sqrt(std::pow(obj->x-obj1->x, 2) + std::pow(obj->y-obj1->y, 2));
}

double inline gravity_force(double mass1, double mass2, double pos1, double pos2, double distance) {
    if (distance == 0.0) return 0.0;
    return -1*G*mass1*mass2*(pos1-pos2)/std::pow(distance, 3); // Teach yourself how this was derived
}

void inline update_particle(particle *obj) {
    obj->Vx += (obj->Fx/obj->mass);
    obj->x += obj->Vx;

    obj->Vy += (obj->Fy/obj->mass);
    obj->y += obj->Vy;
}

int main() {
    std::cout << "START" << std::endl;

    particle* p1 = new particle{};
    particle* p2 = new particle{};
    
    p1->x = -10.0;
    p1->Vy = -1.0;

    p2->x = 10.0;
    p2->Vy = 1.0;

    particle *particles[2] = {p1, p2};

    for (int i = 0; i < 15; i++) {
        for (particle *j: particles) {
            for (particle *k: particles) {
                double r = particle_distance(j, k);

                j->Fx += gravity_force(j->mass, k->mass, j->x, k->x, r);
                j->Fy += gravity_force(j->mass, k->mass, j->y, k->y, r);
            }
            update_particle(j);
            std::cout << "X: " << j->x << ", Y: "<< j->y << std::endl;
        }
    }
    std::cout << "DONE" << std::endl;
}