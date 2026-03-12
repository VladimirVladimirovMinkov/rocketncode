//
// Created by vladi on 3/03/2026.
//

#ifndef CLION_PIPE_H
#define CLION_PIPE_H

#include <array>

#include "tank.h"

class pipe {
public:
    static constexpr int MAX_SEGMENTS = 16;

    pipe();
    pipe(const char *name, int segments, double length_m, double diameter_m, const gas &seed_gas);

    const char *name() const;
    int segment_count() const;
    double length() const;
    double diameter() const;
    double area() const;

    tank &segment(int index);
    const tank &segment(int index) const;

    void init_uniform(const gas &seed_gas);
    void set_total_moles(double moles, const char *process = "ithm");
    void apply_molar_flow(int segment_index, double mol_per_s, double dt, const char *process = "ithm");

private:
    char name_[64];
    int segments_;
    double length_m_;
    double diameter_m_;
    double area_m2_;
    std::array<tank, MAX_SEGMENTS> seg_;
};

#endif //CLION_PIPE_H
