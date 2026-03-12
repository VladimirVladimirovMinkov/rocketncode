//
// Created by vladi on 9/03/2026.
//

#include "tank.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

tank::tank()
    : volume_m3_(0.0),
      gas_count_(0),
      liquid_count_(0),
      solid_count_(0) {
    std::memset(name_, 0, sizeof(name_));
}

tank::tank(const char *name, double volume_m3)
    : volume_m3_(volume_m3),
      gas_count_(0),
      liquid_count_(0),
      solid_count_(0) {
    std::strncpy(name_, name, sizeof(name_) - 1);
    name_[sizeof(name_) - 1] = '\0';
}

const char *tank::name() const {
    return name_;
}

double tank::volume() const {
    return volume_m3_;
}

int tank::add_gas(const gas &g) {
    if (gas_count_ >= MAX_GASES) {
        return -1;
    }
    gases_[gas_count_] = g;
    if (volume_m3_ > 0.0) {
        gases_[gas_count_].dVg(volume_m3_, "ithm");
    }
    return gas_count_++;
}

int tank::add_liquid(const liquid &l) {
    if (liquid_count_ >= MAX_LIQUIDS) {
        return -1;
    }
    liquids_[liquid_count_] = l;
    return liquid_count_++;
}

int tank::add_solid(const solid &s) {
    if (solid_count_ >= MAX_SOLIDS) {
        return -1;
    }
    solids_[solid_count_] = s;
    return solid_count_++;
}

int tank::gas_count() const {
    return gas_count_;
}

int tank::liquid_count() const {
    return liquid_count_;
}

int tank::solid_count() const {
    return solid_count_;
}

gas &tank::gas_at(int index) {
    if (gas_count_ == 0) {
        return gases_[0];
    }
    index = std::max(0, std::min(index, gas_count_ - 1));
    return gases_[index];
}

const gas &tank::gas_at(int index) const {
    if (gas_count_ == 0) {
        return gases_[0];
    }
    index = std::max(0, std::min(index, gas_count_ - 1));
    return gases_[index];
}

double tank::pressure() const {
    if (gas_count_ == 0 || volume_m3_ <= 0.0) {
        return 0.0;
    }
    double total = 0.0;
    for (int i = 0; i < gas_count_; ++i) {
        total += gases_[i].pressure_from_state(volume_m3_, gases_[i].g.Tk, gases_[i].g.ng);
    }
    return total;
}

void tank::print() const {
    std::printf("%s (tank): V=%.6g m^3 | gases=%d liquids=%d solids=%d | P=%.6g Pa\n",
                name_, volume_m3_, gas_count_, liquid_count_, solid_count_, pressure());
    for (int i = 0; i < gas_count_; ++i) {
        gases_[i].print();
    }
}

void tank::set_volume(double volume_m3, const char *process) {
    if (volume_m3 <= 0.0) {
        return;
    }
    volume_m3_ = volume_m3;
    for (int i = 0; i < gas_count_; ++i) {
        gases_[i].dVg(volume_m3_, process);
    }
}

void tank::set_temperature(double temperature_k, const char *process) {
    if (temperature_k <= 0.0) {
        return;
    }
    for (int i = 0; i < gas_count_; ++i) {
        gases_[i].dTk(temperature_k, process);
    }
}
