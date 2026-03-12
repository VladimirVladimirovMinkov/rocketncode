//
// Created by Codex on 12/03/2026.
//

#include "nozzle.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace {
double clampd(const double value, const double lo, const double hi) {
    return std::max(lo, std::min(value, hi));
}
}

nozzle::nozzle()
    : throat_diameter_m_(0.0),
      discharge_coeff_(0.95) {
    std::memset(name_, 0, sizeof(name_));
}

nozzle::nozzle(const char *name, double throat_diameter_m, double discharge_coeff)
    : throat_diameter_m_(throat_diameter_m),
      discharge_coeff_(discharge_coeff) {
    std::strncpy(name_, name, sizeof(name_) - 1);
    name_[sizeof(name_) - 1] = '\0';
}

const char *nozzle::name() const {
    return name_;
}

double nozzle::area() const {
    const double r = 0.5 * throat_diameter_m_;
    return M_PI * r * r;
}

nozzle_flow_t nozzle::molar_flow(const gas &chamber, const gas &ambient) const {
    nozzle_flow_t out{};
    out.mol_per_s = 0.0;
    out.kg_per_s = 0.0;
    out.area_m2 = area();
    out.choked = false;

    if (out.area_m2 <= 0.0) {
        return out;
    }

    const double Pu = chamber.g.Pg;
    const double Pd = ambient.g.Pg;
    if (Pu <= 0.0 || chamber.g.Tk <= 0.0) {
        return out;
    }

    const double gamma = chamber.gamma(chamber.g.Tk);
    const double R = RCONST;
    const double T = chamber.g.Tk;
    const double pr = clampd(Pd / Pu, 0.0, 1.0);
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

    const double mm_kg_per_mol = std::max(1.0e-12, chamber.g.mm / 1000.0);
    out.kg_per_s = mdot;
    out.mol_per_s = mdot / mm_kg_per_mol;
    return out;
}
