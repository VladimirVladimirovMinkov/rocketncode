//
// Created by vladi on 9/03/2026.
//

#include "liquid.h"
#include "ambient.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

namespace {
constexpr double RCONST_LIQUID = 8.31446261815324;
constexpr double KTC_LIQUID = 273.15;
}

double solve_liquid(double P, double V, double n, double T) {
    // Pass a non-positive value for the unknown variable.
    if (P <= 0.0) return (n * RCONST_LIQUID * T) / V;
    if (V <= 0.0) return (n * RCONST_LIQUID * T) / P;
    if (n <= 0.0) return (P * V) / (RCONST_LIQUID * T);
    if (T <= 0.0) return (P * V) / (n * RCONST_LIQUID);
    return 0.0;
}

liquid::liquid(const char *name,
               double mm, double nl, double Vl,
               double Tk, double Pl, const double cg[18]) {
    std::strncpy(l.name, name, sizeof(l.name) - 1);
    l.name[sizeof(l.name) - 1] = '\0';
    l.mm = mm;
    l.nl = nl;
    l.Vl = Vl;
    l.Tk = Tk;
    l.Pl = Pl;

    if (cg) {
        for (int i = 0; i < 18; ++i) {
            l.cg[i] = cg[i];
        }
    } else {
        for (double &i : l.cg) {
            i = 0.0;
        }
    }
}

int liquid::eq_process(const char *a, const char *b) const {
    while (*a && *b) {
        if (std::tolower(static_cast<unsigned char>(*a)) !=
            std::tolower(static_cast<unsigned char>(*b))) {
            return 0;
        }
        ++a;
        ++b;
    }
    return *a == *b;
}

double liquid::mass() const {
    return l.nl * l.mm;
}

void liquid::print() const {
    std::printf("%s:\n", l.name);
    std::printf("Pressure:    Pa:  %.8g Atmo: %.8g Bar: %.8g Tor: %.8g Psi: %.8g\n",
                l.Pl, l.Pl / 101325.0, l.Pl / 100000.0, l.Pl / 133.322, l.Pl / 6894.76);
    std::printf("Moles:       n:   %.8g Kmol: %.8g\n", l.nl, l.nl / 1000.0);
    std::printf("Mass:        g:   %.8g   Kg: %.8g\n", mass(), mass() / 1000.0);
    std::printf("Volume:      m^3: %.8g    L: %.8g\n", l.Vl, l.Vl * 1000.0);
    std::printf("Temperature: K:   %.8g   C: %.8g  F: %.8g\n\n",
                l.Tk, l.Tk - KTC_LIQUID, (l.Tk - KTC_LIQUID) * 1.8 + 32.0);
}

void liquid::dVl(double V, const char *process) {
    if (l.Vl == V) return;

    if (eq_process(process, "ithm")) {
        l.Pl = (l.nl * RCONST_LIQUID * l.Tk) / V;
    } else if (eq_process(process, "ibar")) {
        l.Tk = (l.Pl * V) / (l.nl * RCONST_LIQUID);
    } else if (eq_process(process, "abat")) {
        const double vr = l.Vl / V;
        l.Pl = l.Pl * std::pow(vr, gamma());
        l.Tk = l.Tk * std::pow(vr, gamma() - 1.0);
    } else if (eq_process(process, "cut")) {
        l.nl /= l.Vl / V;
    }

    l.Vl = V;
}

void liquid::dPl(double P, const char *process) {
    if (l.Pl == P) return;

    const double pr = P / l.Pl;
    if (eq_process(process, "ithm")) {
        l.Vl = (l.nl * RCONST_LIQUID * l.Tk) / P;
    } else if (eq_process(process, "abat")) {
        l.Vl = l.Vl * std::pow(1.0 / pr, 1.0 / gamma());
        l.Tk = l.Tk * std::pow(pr, (gamma() - 1.0) / gamma());
    }

    l.Pl = P;
}

void liquid::dnl(double n, const char *process) {
    if (l.nl == n) return;

    if (eq_process(process, "ithm")) {
        l.Pl = (n * RCONST_LIQUID * l.Tk) / l.Vl;
    } else if (eq_process(process, "ibar")) {
        l.Vl = (n * RCONST_LIQUID * l.Tk) / l.Pl;
    }

    l.nl = n;
}

void liquid::dTk(double T, const char *process) {
    if (l.Tk == T) return;

    const double tr = T / l.Tk;
    if (eq_process(process, "ibar")) {
        l.Vl = (l.nl * RCONST_LIQUID * T) / l.Pl;
    } else if (eq_process(process, "abat")) {
        l.Vl = l.Vl * std::pow(1.0 / tr, 1.0 / (gamma() - 1.0));
        l.Pl = l.Pl * std::pow(tr, gamma() / (gamma() - 1.0));
    } else if (eq_process(process, "ivol")) {
        l.Pl = l.Pl * tr;
    }

    l.Tk = T;
}

void liquid::test_stepwise(double target_V, const char *process, int steps) {
    if (steps <= 0) {
        std::printf("steps must be > 0\n");
        return;
    }

    const double total_delta_V = target_V - l.Vl;
    const double step_size = total_delta_V / static_cast<double>(steps);

    for (int i = 0; i < steps; ++i) {
        const double new_v = l.Vl + step_size;
        dVl(new_v, process);
    }
}

double liquid::gamma() {
    int offset = (l.Tk < 1000.0) ? 0 : 9;

    const double Tinv = 1.0 / l.Tk;
    const double Tinv2 = Tinv * Tinv;
    const double T2 = l.Tk * l.Tk;
    const double T3 = T2 * l.Tk;
    const double T4 = T3 * l.Tk;

    const double cp_over_R = (l.cg[offset + 0] * Tinv2) +
                             (l.cg[offset + 1] * Tinv) +
                              l.cg[offset + 2] +
                             (l.cg[offset + 3] * l.Tk) +
                             (l.cg[offset + 4] * T2) +
                             (l.cg[offset + 5] * T3) +
                             (l.cg[offset + 6] * T4);

    return cp_over_R / (cp_over_R - 1.0);
}

void liquid::equalize_to(const ambient &outside, const char *process) {
    dPl(outside.pressure(), process);
}
