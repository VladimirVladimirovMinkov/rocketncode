//
// Created by Codex on 12/03/2026.
//

#ifndef CLION_VALVE_H
#define CLION_VALVE_H

#include "gas.h"

struct valve_flow_t {
    double mol_per_s;
    double kg_per_s;
    double area_m2;
    bool choked;
};

class valve {
public:
    valve();
    valve(const char *name, double diameter_m, double discharge_coeff = 0.9);

    const char *name() const;

    void set_opening(double fraction);
    double opening() const;

    double area() const;

    valve_flow_t molar_flow(const gas &a, const gas &b) const;

private:
    char name_[64];
    double diameter_m_;
    double discharge_coeff_;
    double opening_fraction_;
};

#endif //CLION_VALVE_H
