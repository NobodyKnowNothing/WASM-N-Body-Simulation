#include <iostream>
#include <cmath>

struct particle {
    double x;
    double y;
    double Vx;
    double Vy;
    double Ax;
    double Ay;
    double Fx;
    double Fy;
    double mass;
    double radius;
};

const double G = 0.000000000066743;

double inline particle_distance(particle obj, particle obj1) {
    return std::sqrt(std::pow(obj.x-obj1.x, 2) + std::pow(obj.y-obj1.y, 2));
}

double inline gravity_force(double mass1, double mass2, double pos1, double pos2, double distance) {
    return -1*G*mass1*mass2*(pos1-pos2)/std::pow(distance, 3);
}

void update_particle(particle* obj) {
    *obj.x += *obj.Vx;
    *obj.Vx += *obj.Ax;


}

void main() {

}