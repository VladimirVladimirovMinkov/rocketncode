//
// Created by Codex on 4/03/2026.
//

#include "ambient.h"

#include <cstdio>
#include <cstring>

#include "gas.h"

ambient::ambient(const char *name, double pressure_pa, double temperature_k)
    : pressure_pa_(pressure_pa), temperature_k_(temperature_k) {
    std::strncpy(name_, name, sizeof(name_) - 1);
    name_[sizeof(name_) - 1] = '\0';
}

const char *ambient::name() const {
    return name_;
}

double ambient::pressure() const {
    return pressure_pa_;
}

double ambient::temperature() const {
    return temperature_k_;
}

void ambient::set_temperature(double temperature_k) {
    temperature_k_ = temperature_k;
}

void ambient::print() const {
    std::printf("%s (ambient): P=%.8g Pa, T=%.8g K\n", name_, pressure_pa_, temperature_k_);
}

void ambient::equalize_pressure(gas &sample, const char *process) const {
    sample.dPg(pressure_pa_, process);
}
