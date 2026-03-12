//
// Created by vladi on 9/03/2026.
//

#include "solid.h"
#include "ambient.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

namespace {
constexpr double RCONST_SOLID = 8.31446261815324;
constexpr double KTC_SOLID = 273.15;
}

double solve_solid(double P, double V, double n, double T) {
    // Pass a non-positive value for the unknown variable.
    if (P <= 0.0) return (n * RCONST_SOLID * T) / V;
    if (V <= 0.0) return (n * RCONST_SOLID * T) / P;
    if (n <= 0.0) return (P * V) / (RCONST_SOLID * T);
    if (T <= 0.0) return (P * V) / (n * RCONST_SOLID);
    return 0.0;
}

solid::solid() {
    std::memset(s.name, 0, sizeof(s.name));
    s.Ps = 0.0;
    s.Vs = 0.0;
    s.ns = 0.0;
    s.Tk = 0.0;
    s.mm = 0.0;
    for (double &coef : s.cs) {
        coef = 0.0;
    }
}

solid::solid(const char *name,
             double mm, double ns, double Vs,
             double Tk, double Ps, const double cg[18]) {
    std::strncpy(s.name, name, sizeof(s.name) - 1);
    s.name[sizeof(s.name) - 1] = '\0';
    s.mm = mm;
    s.ns = ns;
    s.Vs = Vs;
    s.Tk = Tk;
    s.Ps = Ps;

    if (cg) {
        for (int i = 0; i < 18; ++i) {
            s.cs[i] = cg[i];
        }
    } else {
        for (double &i : s.cs) {
            i = 0.0;
        }
    }
}

int solid::eq_process(const char *a, const char *b) const {
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

double solid::mass() const {
    return s.ns * s.mm;
}

void solid::print() const {
    std::printf("%s:\n", s.name);
    std::printf("Pressure:    Pa:  %.8g Atmo: %.8g Bar: %.8g Tor: %.8g Psi: %.8g\n",
                s.Ps, s.Ps / 101325.0, s.Ps / 100000.0, s.Ps / 133.322, s.Ps / 6894.76);
    std::printf("Moles:       n:   %.8g Kmol: %.8g\n", s.ns, s.ns / 1000.0);
    std::printf("Mass:        g:   %.8g   Kg: %.8g\n", mass(), mass() / 1000.0);
    std::printf("Volume:      m^3: %.8g    L: %.8g\n", s.Vs, s.Vs * 1000.0);
    std::printf("Temperature: K:   %.8g   C: %.8g  F: %.8g\n\n",
                s.Tk, s.Tk - KTC_SOLID, (s.Tk - KTC_SOLID) * 1.8 + 32.0);
}

void solid::dVs(double V, const char *process) {
    if (s.Vs == V) return;

    if (eq_process(process, "ithm")) {
        s.Ps = (s.ns * RCONST_SOLID * s.Tk) / V;
    } else if (eq_process(process, "ibar")) {
        s.Tk = (s.Ps * V) / (s.ns * RCONST_SOLID);
    } else if (eq_process(process, "abat")) {
        const double vr = s.Vs / V;
        s.Ps = s.Ps * std::pow(vr, gamma());
        s.Tk = s.Tk * std::pow(vr, gamma() - 1.0);
    } else if (eq_process(process, "cut")) {
        s.ns /= s.Vs / V;
    }

    s.Vs = V;
}

void solid::dPs(double P, const char *process) {
    if (s.Ps == P) return;

    const double pr = P / s.Ps;
    if (eq_process(process, "ithm")) {
        s.Vs = (s.ns * RCONST_SOLID * s.Tk) / P;
    } else if (eq_process(process, "abat")) {
        s.Vs = s.Vs * std::pow(1.0 / pr, 1.0 / gamma());
        s.Tk = s.Tk * std::pow(pr, (gamma() - 1.0) / gamma());
    }

    s.Ps = P;
}

void solid::dns(double n, const char *process) {
    if (s.ns == n) return;

    if (eq_process(process, "ithm")) {
        s.Ps = (n * RCONST_SOLID * s.Tk) / s.Vs;
    } else if (eq_process(process, "ibar")) {
        s.Vs = (n * RCONST_SOLID * s.Tk) / s.Ps;
    }

    s.ns = n;
}

void solid::dTk(double T, const char *process) {
    if (s.Tk == T) return;

    const double tr = T / s.Tk;
    if (eq_process(process, "ibar")) {
        s.Vs = (s.ns * RCONST_SOLID * T) / s.Ps;
    } else if (eq_process(process, "abat")) {
        s.Vs = s.Vs * std::pow(1.0 / tr, 1.0 / (gamma() - 1.0));
        s.Ps = s.Ps * std::pow(tr, gamma() / (gamma() - 1.0));
    } else if (eq_process(process, "ivol")) {
        s.Ps = s.Ps * tr;
    }

    s.Tk = T;
}

void solid::test_stepwise(double target_V, const char *process, int steps) {
    if (steps <= 0) {
        std::printf("steps must be > 0\n");
        return;
    }

    const double total_delta_V = target_V - s.Vs;
    const double step_size = total_delta_V / static_cast<double>(steps);

    for (int i = 0; i < steps; ++i) {
        const double new_v = s.Vs + step_size;
        dVs(new_v, process);
    }
}

double solid::gamma() {
    int offset = (s.Tk < 1000.0) ? 0 : 9;

    const double Tinv = 1.0 / s.Tk;
    const double Tinv2 = Tinv * Tinv;
    const double T2 = s.Tk * s.Tk;
    const double T3 = T2 * s.Tk;
    const double T4 = T3 * s.Tk;

    const double cp_over_R = (s.cs[offset + 0] * Tinv2) +
                             (s.cs[offset + 1] * Tinv) +
                              s.cs[offset + 2] +
                             (s.cs[offset + 3] * s.Tk) +
                             (s.cs[offset + 4] * T2) +
                             (s.cs[offset + 5] * T3) +
                             (s.cs[offset + 6] * T4);

    return cp_over_R / (cp_over_R - 1.0);
}

void solid::equalize_to(const ambient &outside, const char *process) {
    dPs(outside.pressure(), process);
}
