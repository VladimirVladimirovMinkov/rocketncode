//
// Created by Codex on 12/03/2026.
//

#include "valve.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace {
double clampd(const double value, const double lo, const double hi) {
    return std::max(lo, std::min(value, hi));
}

double overlap_area_circles(const double radius, const double separation) {
    if (radius <= 0.0) {
        return 0.0;
    }
    const double d = std::max(0.0, separation);
    if (d >= 2.0 * radius) {
        return 0.0;
    }
    if (d <= 0.0) {
        return M_PI * radius * radius;
    }
    const double r2 = radius * radius;
    const double alpha = 2.0 * std::acos(d / (2.0 * radius));
    const double area = 0.5 * r2 * (alpha - std::sin(alpha)) * 2.0;
    return area;
}
}

valve::valve()
    : diameter_m_(0.0),
      discharge_coeff_(0.9),
      opening_fraction_(0.0) {
    std::memset(name_, 0, sizeof(name_));
}

valve::valve(const char *name, double diameter_m, double discharge_coeff)
    : diameter_m_(diameter_m),
      discharge_coeff_(discharge_coeff),
      opening_fraction_(0.0) {
    std::strncpy(name_, name, sizeof(name_) - 1);
    name_[sizeof(name_) - 1] = '\0';
}

const char *valve::name() const {
    return name_;
}

void valve::set_opening(double fraction) {
    opening_fraction_ = clampd(fraction, 0.0, 1.0);
}

double valve::opening() const {
    return opening_fraction_;
}

double valve::area() const {
    const double radius = 0.5 * diameter_m_;
    const double separation = 2.0 * radius * (1.0 - opening_fraction_);
    return overlap_area_circles(radius, separation);
}

valve_flow_t valve::molar_flow(const gas &a, const gas &b) const {
    valve_flow_t out{};
    out.mol_per_s = 0.0;
    out.kg_per_s = 0.0;
    out.area_m2 = area();
    out.choked = false;

    if (out.area_m2 <= 0.0) {
        return out;
    }

    const gas *up = &a;
    const gas *down = &b;
    double Pu = a.g.Pg;
    double Pd = b.g.Pg;
    int sign = 1;

    if (Pd > Pu) {
        std::swap(Pu, Pd);
        std::swap(up, down);
        sign = -1;
    }

    if (Pu <= 0.0 || up->g.Tk <= 0.0) {
        return out;
    }

    const double gamma = up->gamma(up->g.Tk);
    const double R = RCONST;
    const double T = up->g.Tk;
    const double pr = Pd / Pu;
    const double critical = std::pow(2.0 / (gamma + 1.0), gamma / (gamma - 1.0));

    double mdot = 0.0;
    if (pr <= critical) {
        out.choked = true;
        const double factor = std::pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (2.0 * (gamma - 1.0)));
        mdot = discharge_coeff_ * out.area_m2 * Pu * std::sqrt(gamma / (R * T)) * factor;
    } else {
        const double term = (std::pow(pr, 2.0 / gamma) - std::pow(pr, (gamma + 1.0) / gamma));
        if (term > 0.0) {
            mdot = discharge_coeff_ * out.area_m2 * Pu *
                   std::sqrt((2.0 * gamma / (R * T * (gamma - 1.0))) * term);
        }
    }

    const double mm_kg_per_mol = std::max(1.0e-12, up->g.mm / 1000.0);
    out.kg_per_s = sign * mdot;
    out.mol_per_s = sign * (mdot / mm_kg_per_mol);
    return out;
}
