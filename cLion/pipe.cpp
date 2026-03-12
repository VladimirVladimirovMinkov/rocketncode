//
// Created by vladi on 3/03/2026.
//

#include "pipe.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace {
double clampd(const double value, const double lo, const double hi) {
    return std::max(lo, std::min(value, hi));
}
}

pipe::pipe()
    : segments_(0),
      length_m_(0.0),
      diameter_m_(0.0),
      area_m2_(0.0) {
    std::memset(name_, 0, sizeof(name_));
}

pipe::pipe(const char *name, int segments, double length_m, double diameter_m, const gas &seed_gas)
    : segments_(static_cast<int>(clampd(static_cast<double>(segments), 1.0, MAX_SEGMENTS))),
      length_m_(length_m),
      diameter_m_(diameter_m),
      area_m2_(0.0) {
    std::strncpy(name_, name, sizeof(name_) - 1);
    name_[sizeof(name_) - 1] = '\0';
    const double radius = 0.5 * diameter_m_;
    area_m2_ = M_PI * radius * radius;
    init_uniform(seed_gas);
}

const char *pipe::name() const {
    return name_;
}

int pipe::segment_count() const {
    return segments_;
}

double pipe::length() const {
    return length_m_;
}

double pipe::diameter() const {
    return diameter_m_;
}

double pipe::area() const {
    return area_m2_;
}

tank &pipe::segment(int index) {
    if (segments_ == 0) {
        return seg_[0];
    }
    index = std::max(0, std::min(index, segments_ - 1));
    return seg_[index];
}

const tank &pipe::segment(int index) const {
    if (segments_ == 0) {
        return seg_[0];
    }
    index = std::max(0, std::min(index, segments_ - 1));
    return seg_[index];
}

void pipe::init_uniform(const gas &seed_gas) {
    if (segments_ <= 0) {
        return;
    }
    const double segment_volume = (segments_ > 0 && length_m_ > 0.0)
        ? (area_m2_ * length_m_ / static_cast<double>(segments_))
        : 0.0;
    for (int i = 0; i < segments_; ++i) {
        seg_[i] = tank(name_, segment_volume);
        seg_[i].add_gas(seed_gas);
    }
}

void pipe::set_total_moles(double moles, const char *process) {
    if (segments_ <= 0) {
        return;
    }
    const double per_segment = moles / static_cast<double>(segments_);
    for (int i = 0; i < segments_; ++i) {
        if (seg_[i].gas_count() > 0) {
            seg_[i].gas_at(0).dng(per_segment, process);
        }
    }
}

void pipe::apply_molar_flow(int segment_index, double mol_per_s, double dt, const char *process) {
    if (segments_ <= 0) {
        return;
    }
    segment_index = std::max(0, std::min(segment_index, segments_ - 1));
    if (seg_[segment_index].gas_count() == 0) {
        return;
    }
    gas &g = seg_[segment_index].gas_at(0);
    const double n_new = std::max(0.0, g.g.ng + mol_per_s * dt);
    g.dng(n_new, process);
}

