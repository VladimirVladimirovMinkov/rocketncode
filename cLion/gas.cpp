//
// Created by vladi on 3/03/2026.
//

#include "gas.h"
#include "ambient.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

namespace {
constexpr double MIN_POSITIVE = 1.0e-12;
constexpr double MIN_TEMPERATURE_K = 1.0;
constexpr double NASA_T_MIN = 200.0;
constexpr double NASA_T_MAX = 6000.0;
constexpr int EOS_SOLVE_STEPS = 96;

struct pr_params_t {
    double a_alpha;
    double b;
    bool valid;
};

double clampd(const double value, const double lo, const double hi) {
    return std::max(lo, std::min(value, hi));
}

bool build_pr_params(const gas_t &g, const double temperature_k, pr_params_t &out) {
    out = {0.0, 0.0, false};

    if (g.eos != gas_eos_t::peng_robinson) {
        return false;
    }

    if (g.tc <= 0.0 || g.pc <= 0.0 || temperature_k <= 0.0) {
        return false;
    }

    const double a = 0.45724 * (RCONST * RCONST) * (g.tc * g.tc) / g.pc;
    const double b = 0.07780 * RCONST * g.tc / g.pc;
    const double kappa = 0.37464 + 1.54226 * g.omega - 0.26992 * g.omega * g.omega;
    const double tr_sqrt = std::sqrt(std::max(temperature_k / g.tc, MIN_POSITIVE));
    const double alpha = std::pow(1.0 + kappa * (1.0 - tr_sqrt), 2.0);

    out.a_alpha = a * alpha;
    out.b = b;
    out.valid = true;
    return true;
}

double pressure_pr_from_vm(const pr_params_t &params, const double temperature_k, const double vm_m3_per_mol) {
    const double denom1 = vm_m3_per_mol - params.b;
    const double denom2 = (vm_m3_per_mol * vm_m3_per_mol) +
                          (2.0 * params.b * vm_m3_per_mol) -
                          (params.b * params.b);

    if (denom1 <= MIN_POSITIVE || denom2 <= MIN_POSITIVE) {
        return 1.0e12;
    }

    return (RCONST * temperature_k) / denom1 - params.a_alpha / denom2;
}

double dpressure_pr_dvm(const pr_params_t &params, const double temperature_k, const double vm_m3_per_mol) {
    const double denom1 = vm_m3_per_mol - params.b;
    const double denom2 = (vm_m3_per_mol * vm_m3_per_mol) +
                          (2.0 * params.b * vm_m3_per_mol) -
                          (params.b * params.b);

    if (denom1 <= MIN_POSITIVE || denom2 <= MIN_POSITIVE) {
        return -1.0e12;
    }

    const double term1 = -(RCONST * temperature_k) / (denom1 * denom1);
    const double term2 = params.a_alpha * (2.0 * vm_m3_per_mol + 2.0 * params.b) / (denom2 * denom2);
    return term1 + term2;
}

double solve_vm_pr(const pr_params_t &params,
                   const double pressure_pa,
                   const double temperature_k) {
    const double p = std::max(pressure_pa, MIN_POSITIVE);
    const double b_floor = params.b * 1.000001;
    double vm = std::max((RCONST * temperature_k) / p, b_floor * 1.01);

    for (int i = 0; i < EOS_SOLVE_STEPS; ++i) {
        const double f = pressure_pr_from_vm(params, temperature_k, vm) - p;
        const double dfdv = dpressure_pr_dvm(params, temperature_k, vm);

        if (std::abs(f) < std::max(1.0, p) * 1.0e-10) {
            return vm;
        }

        if (std::abs(dfdv) < 1.0e-18) {
            break;
        }

        double candidate = vm - f / dfdv;
        if (candidate <= b_floor || !std::isfinite(candidate)) {
            candidate = 0.5 * (vm + b_floor);
        }
        vm = candidate;
    }

    // Robust fallback bracket + bisection.
    double lo = b_floor;
    double hi = std::max(vm, (RCONST * temperature_k) / p);

    if (pressure_pr_from_vm(params, temperature_k, hi) > p) {
        for (int i = 0; i < 48; ++i) {
            hi *= 2.0;
            if (pressure_pr_from_vm(params, temperature_k, hi) <= p) {
                break;
            }
        }
    }

    for (int i = 0; i < EOS_SOLVE_STEPS; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double pmid = pressure_pr_from_vm(params, temperature_k, mid);
        if (pmid > p) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return 0.5 * (lo + hi);
}
}

double solve(double P, double V, double n, double T) {
    // Pass a non-positive value for the unknown variable.
    if (P <= 0.0) return (n * RCONST * T) / std::max(V, MIN_POSITIVE);
    if (V <= 0.0) return (n * RCONST * T) / std::max(P, MIN_POSITIVE);
    if (n <= 0.0) return (P * V) / (RCONST * std::max(T, MIN_POSITIVE));
    if (T <= 0.0) return (P * V) / (n * RCONST);
    return 0.0;
}

gas::gas() {
    std::memset(g.name, 0, sizeof(g.name));
    g.Pg = 0.0;
    g.Vg = 0.0;
    g.ng = 0.0;
    g.Tk = 0.0;
    g.mm = 0.0;
    g.eos = gas_eos_t::ideal;
    g.tc = 0.0;
    g.pc = 0.0;
    g.omega = 0.0;
    for (double &coef : g.cg) {
        coef = 0.0;
    }
}

gas::gas(const char *name,
         double mm, double ng, double Vg,
         double Tk, double Pg, const double cg[18])
    : gas(name, mm, ng, Vg, Tk, Pg, cg, gas_eos_t::ideal, 0.0, 0.0, 0.0) {
}

gas::gas(const char *name,
         double mm, double ng, double Vg,
         double Tk, double Pg, const double cg[18],
         gas_eos_t eos, double tc, double pc, double omega) {
    std::strncpy(g.name, name, sizeof(g.name) - 1);
    g.name[sizeof(g.name) - 1] = '\0';
    g.mm = mm;
    g.ng = ng;
    g.Vg = Vg;
    g.Tk = Tk;
    g.Pg = Pg;
    g.eos = eos;
    g.tc = tc;
    g.pc = pc;
    g.omega = omega;

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

double gas::cp_over_R() const {
    return cp_over_R(g.Tk);
}

double gas::cp_over_R(double temperature_k) const {
    const double t = clampd(temperature_k, NASA_T_MIN, NASA_T_MAX);
    const int offset = (t < 1000.0) ? 0 : 9;

    const double tinv = 1.0 / t;
    const double tinv2 = tinv * tinv;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double t4 = t3 * t;

    const double cp = (g.cg[offset + 0] * tinv2) +
                      (g.cg[offset + 1] * tinv) +
                      g.cg[offset + 2] +
                      (g.cg[offset + 3] * t) +
                      (g.cg[offset + 4] * t2) +
                      (g.cg[offset + 5] * t3) +
                      (g.cg[offset + 6] * t4);

    return std::max(cp, 1.001);
}

double gas::h_over_RT(double temperature_k) const {
    const double t = clampd(temperature_k, NASA_T_MIN, NASA_T_MAX);
    const int offset = (t < 1000.0) ? 0 : 9;

    const double tinv = 1.0 / t;
    const double tinv2 = tinv * tinv;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double t4 = t3 * t;
    const double lnT = std::log(t);

    return (-g.cg[offset + 0] * tinv2) +
           (g.cg[offset + 1] * lnT * tinv) +
           g.cg[offset + 2] +
           (g.cg[offset + 3] * t * 0.5) +
           (g.cg[offset + 4] * t2 / 3.0) +
           (g.cg[offset + 5] * t3 * 0.25) +
           (g.cg[offset + 6] * t4 * 0.2) +
           (g.cg[offset + 7] * tinv);
}

double gas::s_over_R(double temperature_k) const {
    const double t = clampd(temperature_k, NASA_T_MIN, NASA_T_MAX);
    const int offset = (t < 1000.0) ? 0 : 9;

    const double tinv = 1.0 / t;
    const double tinv2 = tinv * tinv;
    const double t2 = t * t;
    const double t3 = t2 * t;
    const double t4 = t3 * t;
    const double lnT = std::log(t);

    return (-0.5 * g.cg[offset + 0] * tinv2) +
           (-g.cg[offset + 1] * tinv) +
           (g.cg[offset + 2] * lnT) +
           (g.cg[offset + 3] * t) +
           (g.cg[offset + 4] * t2 * 0.5) +
           (g.cg[offset + 5] * t3 / 3.0) +
           (g.cg[offset + 6] * t4 * 0.25) +
           g.cg[offset + 8];
}

double gas::cp_molar(double temperature_k) const {
    return cp_over_R(temperature_k) * RCONST;
}

double gas::cv_molar(double temperature_k) const {
    return std::max(cp_molar(temperature_k) - RCONST, 1.0e-9);
}

double gas::enthalpy_molar(double temperature_k) const {
    return h_over_RT(temperature_k) * RCONST * temperature_k;
}

double gas::entropy_molar(double temperature_k) const {
    return s_over_R(temperature_k) * RCONST;
}

double gas::pressure_from_state(double volume_m3, double temperature_k, double moles) const {
    if (volume_m3 <= 0.0 || temperature_k <= 0.0 || moles <= 0.0) {
        return 0.0;
    }

    if (g.eos == gas_eos_t::ideal) {
        return (moles * RCONST * temperature_k) / volume_m3;
    }

    pr_params_t params{};
    if (!build_pr_params(g, temperature_k, params)) {
        return (moles * RCONST * temperature_k) / volume_m3;
    }

    const double vm = volume_m3 / moles;
    return pressure_pr_from_vm(params, temperature_k, vm);
}

double gas::volume_from_state(double pressure_pa, double temperature_k, double moles) const {
    if (pressure_pa <= 0.0 || temperature_k <= 0.0 || moles <= 0.0) {
        return 0.0;
    }

    if (g.eos == gas_eos_t::ideal) {
        return (moles * RCONST * temperature_k) / pressure_pa;
    }

    pr_params_t params{};
    if (!build_pr_params(g, temperature_k, params)) {
        return (moles * RCONST * temperature_k) / pressure_pa;
    }

    const double vm = solve_vm_pr(params, pressure_pa, temperature_k);
    return vm * moles;
}

double gas::temperature_from_state(double pressure_pa, double volume_m3, double moles) const {
    if (pressure_pa <= 0.0 || volume_m3 <= 0.0 || moles <= 0.0) {
        return g.Tk;
    }

    if (g.eos == gas_eos_t::ideal) {
        return (pressure_pa * volume_m3) / (moles * RCONST);
    }

    double t = std::max(g.Tk, MIN_TEMPERATURE_K);
    for (int i = 0; i < EOS_SOLVE_STEPS; ++i) {
        const double f = pressure_from_state(volume_m3, t, moles) - pressure_pa;
        if (std::abs(f) < std::max(1.0, pressure_pa) * 1.0e-9) {
            return t;
        }

        const double dt = std::max(1.0e-3, 1.0e-6 * t);
        const double f1 = pressure_from_state(volume_m3, t + dt, moles) - pressure_pa;
        const double dfdT = (f1 - f) / dt;
        if (std::abs(dfdT) < 1.0e-12) {
            break;
        }

        t -= f / dfdT;
        if (!std::isfinite(t) || t < MIN_TEMPERATURE_K) {
            t = MIN_TEMPERATURE_K;
        }
    }

    return t;
}

void gas::print() const {
    std::printf("%s:\n", g.name);
    std::printf("EOS:         %s\n", (g.eos == gas_eos_t::peng_robinson) ? "Peng-Robinson" : "Ideal");
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
    if (V <= 0.0) return;

    if (eq_process(process, "ithm")) {
        g.Pg = pressure_from_state(V, g.Tk, g.ng);
    } else if (eq_process(process, "ibar")) {
        g.Tk = temperature_from_state(g.Pg, V, g.ng);
    } else if (eq_process(process, "abat")) {
        const int steps = 64;
        const double v_start = g.Vg;
        for (int i = 0; i < steps; ++i) {
            const double f0 = static_cast<double>(i) / static_cast<double>(steps);
            const double f1 = static_cast<double>(i + 1) / static_cast<double>(steps);
            const double v0 = v_start + (V - v_start) * f0;
            const double v1 = v_start + (V - v_start) * f1;
            const double vr = std::max(v0 / std::max(v1, MIN_POSITIVE), MIN_POSITIVE);
            const double gm = gamma();
            g.Tk = g.Tk * std::pow(vr, gm - 1.0);
        }
        g.Pg = pressure_from_state(V, g.Tk, g.ng);
    } else if (eq_process(process, "cut")) {
        g.ng /= g.Vg / V;
    }

    g.Vg = V;
}

void gas::dPg(double P, const char *process) {
    if (g.Pg == P) return;
    if (P <= 0.0) return;

    const double p_start = g.Pg;
    if (eq_process(process, "ithm")) {
        g.Vg = volume_from_state(P, g.Tk, g.ng);
    } else if (eq_process(process, "abat")) {
        const int steps = 64;
        for (int i = 0; i < steps; ++i) {
            const double p0 = p_start + (P - p_start) * (static_cast<double>(i) / steps);
            const double p1 = p_start + (P - p_start) * (static_cast<double>(i + 1) / steps);
            const double pr = std::max(p1 / std::max(p0, MIN_POSITIVE), MIN_POSITIVE);
            const double gm = gamma();
            g.Tk = g.Tk * std::pow(pr, (gm - 1.0) / gm);
        }
        g.Vg = volume_from_state(P, g.Tk, g.ng);
    }

    g.Pg = P;
}

void gas::dng(double n, const char *process) {
    if (g.ng == n) return;
    if (n <= 0.0) return;

    if (eq_process(process, "ithm")) {
        g.Pg = pressure_from_state(g.Vg, g.Tk, n);
    } else if (eq_process(process, "ibar")) {
        g.Vg = volume_from_state(g.Pg, g.Tk, n);
    }

    g.ng = n;
}

void gas::dTk(double T, const char *process) {
    if (g.Tk == T) return;
    if (T <= 0.0) return;

    if (eq_process(process, "ibar")) {
        g.Vg = volume_from_state(g.Pg, T, g.ng);
    } else if (eq_process(process, "abat")) {
        const int steps = 64;
        const double t_start = g.Tk;
        for (int i = 0; i < steps; ++i) {
            const double t0 = t_start + (T - t_start) * (static_cast<double>(i) / steps);
            const double t1 = t_start + (T - t_start) * (static_cast<double>(i + 1) / steps);
            const double tr = std::max(t1 / std::max(t0, MIN_POSITIVE), MIN_POSITIVE);
            const double gm = gamma();
            g.Vg = g.Vg * std::pow(1.0 / tr, 1.0 / std::max(gm - 1.0, 1.0e-6));
            g.Tk = t1;
        }
        g.Pg = pressure_from_state(g.Vg, g.Tk, g.ng);
        return;
    } else if (eq_process(process, "ivol")) {
        g.Pg = pressure_from_state(g.Vg, T, g.ng);
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
    return gamma(g.Tk);
}

double gas::gamma(double temperature_k) const {
    const double cp = cp_over_R(temperature_k);
    const double denom = std::max(cp - 1.0, 1.0e-6);
    return clampd(cp / denom, 1.01, 3.0);
}

void gas::equalize_to(const ambient &outside, const char *process) {
    dPg(outside.pressure(), process);
}
