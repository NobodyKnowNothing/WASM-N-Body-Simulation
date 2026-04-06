#ifndef BARNES_HUT_H
#define BARNES_HUT_H

#include <cmath>
#include "newtonian-dynamics.hpp"

struct qtnode {
    double CoMx = 0;
    double CoMy = 0;
    double totalMass = 0;
    double length = 0;
    qtnode* NW = nullptr;
    qtnode* NE = nullptr;
    qtnode* SW = nullptr;
    qtnode* SE = nullptr;
};

#endif