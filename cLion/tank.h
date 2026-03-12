//
// Created by vladi on 9/03/2026.
//

#ifndef CLION_TANK_H
#define CLION_TANK_H

#include <array>

#include "gas.h"
#include "liquid.h"
#include "solid.h"

class tank {
public:
    static constexpr int MAX_GASES = 8;
    static constexpr int MAX_LIQUIDS = 4;
    static constexpr int MAX_SOLIDS = 4;

    tank();
    tank(const char *name, double volume_m3);

    const char *name() const;
    double volume() const;

    int add_gas(const gas &g);
    int add_liquid(const liquid &l);
    int add_solid(const solid &s);

    int gas_count() const;
    int liquid_count() const;
    int solid_count() const;

    gas &gas_at(int index);
    const gas &gas_at(int index) const;

    double pressure() const;

    void print() const;

    void set_volume(double volume_m3, const char *process = "ithm");
    void set_temperature(double temperature_k, const char *process = "ivol");

private:
    char name_[64];
    double volume_m3_;
    std::array<gas, MAX_GASES> gases_;
    int gas_count_;
    std::array<liquid, MAX_LIQUIDS> liquids_;
    int liquid_count_;
    std::array<solid, MAX_SOLIDS> solids_;
    int solid_count_;
};

#endif //CLION_TANK_H
