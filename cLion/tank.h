//
// Created by vladi on 9/03/2026.
//

#ifndef CLION_TANK_H
#define CLION_TANK_H
#include "gas.h"
#include "liquid.h"
#include "solid.h"


class tank {
    gas gases[0];
    liquid liquids[0];
    solid solids[0];
    double dimensions[3]; // x,y,z of a box, assume box
};


#endif //CLION_TANK_H
