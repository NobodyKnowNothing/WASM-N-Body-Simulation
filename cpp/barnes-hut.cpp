#include <cmath>
#include <span>
#include <algorithm>
#include "newtonian-dynamics.hpp"

struct qtnode {
    double CoMx = 0;
    double CoMy = 0;
    double cx = 0;
    double cy = 0;
    double totalMass = 0;
    double length = 0;
    qtnode* NW = nullptr;
    qtnode* NE = nullptr;
    qtnode* SW = nullptr;
    qtnode* SE = nullptr;
};

inline void init_qtroot(std::span<particle*> particles) {
    auto [min_x, max_x] = std::minmax_element(particles.begin(), particles.end(), [](const particle* a, const particle* b) {return a->x < b->x;});
    
    auto [min_y, max_y] = std::minmax_element(particles.begin(), particles.end(), [](const particle* a, const particle* b) {return a->y <    b->y;});

}