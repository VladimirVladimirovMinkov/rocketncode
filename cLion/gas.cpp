//
// Created by vladi on 3/03/2026.
//

#include "gas.h"
#include "ambient.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

double solve(double P, double V, double n, double T) {
    // Pass a non-positive value for the unknown variable.
    if (P <= 0.0) return (n * RCONST * T) / V;
    if (V <= 0.0) return (n * RCONST * T) / P;
    if (n <= 0.0) return (P * V) / (RCONST * T);
    if (T <= 0.0) return (P * V) / (n * RCONST);
    return 0.0;
}

gas::gas(const char *name,
         double mm, double ng, double Vg,
         double Tk, double Pg, const double cg[18]) {
    std::strncpy(g.name, name, sizeof(g.name) - 1);
    g.name[sizeof(g.name) - 1] = '\0';
    g.mm = mm;
    g.ng = ng;
    g.Vg = Vg;
    g.Tk = Tk;
    g.Pg = Pg;

    if (cg) {
        for (int i = 0; i < 18; ++i) {
            g.cg[i] = cg[i];
        }
    } else {
        for (double &i : g.cg) {
            i = 0.0;
        }
    }
}

int gas::eq_process(const char *a, const char *b) const {
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

double gas::mass() const {
    return g.ng * g.mm;
}

void gas::print() const {
    std::printf("%s:\n", g.name);
    std::printf("Pressure:    Pa:  %.8g Atmo: %.8g Bar: %.8g Tor: %.8g Psi: %.8g\n",
                g.Pg, g.Pg / 101325.0, g.Pg / 100000.0, g.Pg / 133.322, g.Pg / 6894.76);
    std::printf("Moles:       n:   %.8g Kmol: %.8g\n", g.ng, g.ng / 1000.0);
    std::printf("Mass:        g:   %.8g   Kg: %.8g\n", mass(), mass() / 1000.0);
    std::printf("Volume:      m^3: %.8g    L: %.8g\n", g.Vg, g.Vg * 1000.0);
    std::printf("Temperature: K:   %.8g   C: %.8g  F: %.8g\n\n",
                g.Tk, g.Tk - KTC, (g.Tk - KTC) * 1.8 + 32.0);
}

void gas::dVg(double V, const char *process) {
    if (g.Vg == V) return;

    if (eq_process(process, "ithm")) {
        g.Pg = (g.ng * RCONST * g.Tk) / V;
    } else if (eq_process(process, "ibar")) {
        g.Tk = (g.Pg * V) / (g.ng * RCONST);
    } else if (eq_process(process, "abat")) {
        const double vr = g.Vg / V;
        g.Pg = g.Pg * std::pow(vr, gamma());
        g.Tk = g.Tk * std::pow(vr, gamma() - 1.0);
    } else if (eq_process(process, "cut")) {
        g.ng /= g.Vg / V;
    }

    g.Vg = V;
}

void gas::dPg(double P, const char *process) {
    if (g.Pg == P) return;

    const double pr = P / g.Pg;
    if (eq_process(process, "ithm")) {
        g.Vg = (g.ng * RCONST * g.Tk) / P;
    } else if (eq_process(process, "abat")) {
        g.Vg = g.Vg * std::pow(1.0 / pr, 1.0 / gamma());
        g.Tk = g.Tk * std::pow(pr, (gamma() - 1.0) / gamma());
    }

    g.Pg = P;
}

void gas::dng(double n, const char *process) {
    if (g.ng == n) return;

    if (eq_process(process, "ithm")) {
        g.Pg = (n * RCONST * g.Tk) / g.Vg;
    } else if (eq_process(process, "ibar")) {
        g.Vg = (n * RCONST * g.Tk) / g.Pg;
    }

    g.ng = n;
}

void gas::dTk(double T, const char *process) {
    if (g.Tk == T) return;

    const double tr = T / g.Tk;
    if (eq_process(process, "ibar")) {
        g.Vg = (g.ng * RCONST * T) / g.Pg;
    } else if (eq_process(process, "abat")) {
        g.Vg = g.Vg * std::pow(1.0 / tr, 1.0 / (gamma() - 1.0));
        g.Pg = g.Pg * std::pow(tr, gamma() / (gamma() - 1.0));
    } else if (eq_process(process,"ivol")) {
        g.Pg = g.Pg * tr;
    }

    g.Tk = T;
}

void gas::test_stepwise(double target_V, const char *process, int steps) {
    if (steps <= 0) {
        std::printf("steps must be > 0\n");
        return;
    }

    const double total_delta_V = target_V - g.Vg;
    const double step_size = total_delta_V / static_cast<double>(steps);

    for (int i = 0; i < steps; ++i) {
        const double new_v = g.Vg + step_size;
        dVg(new_v, process);
    }
}

double gas::gamma() {
    int offset = (g.Tk < 1000.0) ? 0 : 9;

    // 2. Pre-calculate powers of T
    const double Tinv = 1.0 / g.Tk;
    const double Tinv2 = Tinv * Tinv;
    const double T2 = g.Tk * g.Tk;
    const double T3 = T2 * g.Tk;
    const double T4 = T3 * g.Tk;

    // 3. Calculate dimensionless Heat Capacity (Cp/R)
    // Using g.cg[offset + 0] through g.cg[offset + 6]
    const double cp_over_R = (g.cg[offset + 0] * Tinv2) +
                             (g.cg[offset + 1] * Tinv) +
                              g.cg[offset + 2] +
                             (g.cg[offset + 3] * g.Tk) +
                             (g.cg[offset + 4] * T2) +
                             (g.cg[offset + 5] * T3) +
                             (g.cg[offset + 6] * T4);

    // Solve Gamma from the ratio
    return cp_over_R / (cp_over_R - 1.0);
}

void gas::equalize_to(const ambient &outside, const char *process) {
    dPg(outside.pressure(), process);
}
