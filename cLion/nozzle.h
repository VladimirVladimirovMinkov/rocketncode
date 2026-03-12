//
// Created by Codex on 12/03/2026.
//

#ifndef CLION_NOZZLE_H
#define CLION_NOZZLE_H

#include "gas.h"

struct nozzle_flow_t {
    double mol_per_s;
    double kg_per_s;
    double area_m2;
    bool choked;
};

class nozzle {
public:
    nozzle();
    nozzle(const char *name, double throat_diameter_m, double discharge_coeff = 0.95);

    const char *name() const;
    double area() const;

    nozzle_flow_t molar_flow(const gas &chamber, const gas &ambient) const;

private:
    char name_[64];
    double throat_diameter_m_;
    double discharge_coeff_;
};

#endif //CLION_NOZZLE_H
