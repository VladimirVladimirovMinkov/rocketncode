//
// Created by vladi on 9/03/2026.
//

#include "liquid.h"
#include "ambient.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

namespace {
constexpr double KTC_LIQUID = 273.15;
constexpr double PA_PER_MPA = 1.0e6;
constexpr double MPA_PER_PA = 1.0e-6;
constexpr double MIN_THERMAL_SCALE = 0.05;
constexpr double MIN_COMPRESS_SCALE = 0.2;
constexpr double MAX_SOLVE_PRESSURE_PA = 1.0e9;
constexpr int PRESSURE_SOLVE_STEPS = 96;

double clampd(const double value, const double lo, const double hi) {
    return std::max(lo, std::min(value, hi));
}

double eval_cubic_cl(const double cl[6], const double p_mpa) {
    const double x = p_mpa - cl[0];
    return ((cl[5] * x + cl[4]) * x + cl[3]) * x + cl[2];
}

double eval_cubic_cl_derivative(const double cl[6], const double p_mpa) {
    const double x = p_mpa - cl[0];
    return (3.0 * cl[5] * x + 2.0 * cl[4]) * x + cl[3];
}
}

double solve_liquid(double P, double V, double n, double T) {
    // Retained for compatibility with gas-like interfaces.
    if (P <= 0.0 && V > 0.0) return 101325.0;
    if (V <= 0.0 && n > 0.0) return 3.0e-5 * n;
    if (n <= 0.0 && V > 0.0) return V / 3.0e-5;
    if (T <= 0.0) return KTC_LIQUID;
    return 0.0;
}

liquid::liquid() {
    std::memset(l.name, 0, sizeof(l.name));
    l.Pl = 101325.0;
    l.Vl = 0.0;
    l.nl = 0.0;
    l.Tk = 273.15;
    l.mm = 0.0;
    l.model.tref_k = 273.15;
    l.model.thermal_alpha_per_k = 0.0;
    l.model.compress_temp_per_k = 0.0;
    l.model.dead_space_fraction_ref = 0.0;
    l.model.dead_space_mm_ref = 1.0;
    l.model.vm_ref_m3_per_mol = 1.0e-6;
    for (double &coef : l.model.cl) {
        coef = 0.0;
    }
}

liquid::liquid(const char *name,
               double mm, double nl, double Vl,
               double Tk, double Pl, const liquid_model_t &model) {
    std::strncpy(l.name, name, sizeof(l.name) - 1);
    l.name[sizeof(l.name) - 1] = '\0';
    l.mm = mm;
    l.nl = nl;
    l.Tk = (Tk > 0.0) ? Tk : model.tref_k;
    l.Pl = (Pl > 0.0) ? Pl : 101325.0;
    l.model = model;

    if (Vl > 0.0) {
        l.Vl = Vl;
    } else {
        l.Vl = volume_from_pressure(l.Pl, l.Tk, l.nl);
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

double liquid::compression_fraction(double pressure_pa, double temperature_k) const {
    const double p_mpa = std::max(0.0, pressure_pa * MPA_PER_PA);
    const double p1 = l.model.cl[1];

    double c_ref = 0.0;
    if (p_mpa <= p1) {
        c_ref = eval_cubic_cl(l.model.cl, p_mpa);
    } else {
        const double c_end = eval_cubic_cl(l.model.cl, p1);
        const double slope_end = std::max(0.0, eval_cubic_cl_derivative(l.model.cl, p1));
        c_ref = c_end + slope_end * (p_mpa - p1);
    }

    const double t_scale = clampd(
        1.0 + l.model.compress_temp_per_k * (temperature_k - l.model.tref_k),
        MIN_COMPRESS_SCALE, 3.0);

    const double mm_ref = std::max(1.0e-9, l.model.dead_space_mm_ref);
    const double mm_ratio = std::max(1.0e-9, l.mm) / mm_ref;
    const double dead_space = clampd(
        l.model.dead_space_fraction_ref * std::cbrt(mm_ratio),
        0.02, 0.90);
    const double max_compression = 1.0 - dead_space;

    const double c_temp = c_ref * t_scale;
    return clampd(c_temp, 0.0, max_compression);
}

double liquid::volume_from_pressure(double pressure_pa, double temperature_k, double moles) const {
    if (moles <= 0.0) {
        return 0.0;
    }

    const double thermal_scale = clampd(
        1.0 + l.model.thermal_alpha_per_k * (temperature_k - l.model.tref_k),
        MIN_THERMAL_SCALE, 4.0);
    const double c = compression_fraction(pressure_pa, temperature_k);
    return moles * l.model.vm_ref_m3_per_mol * thermal_scale * (1.0 - c);
}

double liquid::pressure_from_volume(double volume_m3, double temperature_k, double moles) const {
    if (volume_m3 <= 0.0) {
        return MAX_SOLVE_PRESSURE_PA;
    }
    if (moles <= 0.0) {
        return 0.0;
    }

    const auto volume_at = [&](const double pressure_pa) {
        return volume_from_pressure(pressure_pa, temperature_k, moles);
    };

    const double v_at_zero = volume_at(0.0);
    if (volume_m3 >= v_at_zero) {
        return 0.0;
    }

    double lo = 0.0;
    double hi = 100.0 * PA_PER_MPA;
    double v_hi = volume_at(hi);
    while (volume_m3 < v_hi && hi < MAX_SOLVE_PRESSURE_PA) {
        lo = hi;
        hi *= 2.0;
        v_hi = volume_at(hi);
    }

    if (volume_m3 < v_hi) {
        return hi;
    }

    for (int i = 0; i < PRESSURE_SOLVE_STEPS; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double v_mid = volume_at(mid);
        if (v_mid > volume_m3) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return 0.5 * (lo + hi);
}

void liquid::print() const {
    const double c = compression_fraction(l.Pl, l.Tk);
    std::printf("%s:\n", l.name);
    std::printf("Pressure:    Pa:  %.8g Atmo: %.8g Bar: %.8g Tor: %.8g Psi: %.8g\n",
                l.Pl, l.Pl / 101325.0, l.Pl / 100000.0, l.Pl / 133.322, l.Pl / 6894.76);
    std::printf("Moles:       n:   %.8g Kmol: %.8g\n", l.nl, l.nl / 1000.0);
    std::printf("Mass:        g:   %.8g   Kg: %.8g\n", mass(), mass() / 1000.0);
    std::printf("Volume:      m^3: %.8g    L: %.8g\n", l.Vl, l.Vl * 1000.0);
    std::printf("Temperature: K:   %.8g   C: %.8g  F: %.8g\n",
                l.Tk, l.Tk - KTC_LIQUID, (l.Tk - KTC_LIQUID) * 1.8 + 32.0);
    std::printf("Compression: frac %.8g pct %.8g\n", c, c * 100.0);
    std::printf("Model vars:  Tref=%.4gK alpha=%.6g betaT=%.6g deadSpaceRef=%.6g mmRef=%.6g VmRef=%.6g\n\n",
                l.model.tref_k,
                l.model.thermal_alpha_per_k,
                l.model.compress_temp_per_k,
                l.model.dead_space_fraction_ref,
                l.model.dead_space_mm_ref,
                l.model.vm_ref_m3_per_mol);
}

void liquid::print_spline() const {
    std::printf("Spline equations for %s (fractional compression c):\n", l.name);
    std::printf("c(P,T) = clamp( c_ref(P) * (1 + betaT*(T-Tref)), 0, 1-dead_space(mm) )\n");
    std::printf("dead_space(mm) = clamp(deadSpaceRef * cbrt(mm/mmRef), 0.02, 0.90)\n");
    std::printf("  cl = {%.12g, %.12g, %.12g, %.12g, %.12g, %.12g}\n",
                l.model.cl[0], l.model.cl[1], l.model.cl[2],
                l.model.cl[3], l.model.cl[4], l.model.cl[5]);
    std::printf("  %.0f-%.0f MPa: x=P_MPa-%.0f | c_ref=%.12g + %.12g*x + %.12g*x^2 + %.12g*x^3\n",
                l.model.cl[0], l.model.cl[1], l.model.cl[0],
                l.model.cl[2], l.model.cl[3], l.model.cl[4], l.model.cl[5]);
    std::printf("\n");
}

void liquid::dVl(double V, const char *process) {
    if (l.Vl == V) return;
    if (V <= 0.0) return;

    const double prev_p = l.Pl;
    if (eq_process(process, "ithm")) {
        l.Pl = pressure_from_volume(V, l.Tk, l.nl);
    } else if (eq_process(process, "ibar")) {
        const double current_comp = compression_fraction(l.Pl, l.Tk);
        const double denom = l.nl * l.model.vm_ref_m3_per_mol * (1.0 - current_comp);
        if (denom > 0.0 && std::abs(l.model.thermal_alpha_per_k) > 1.0e-12) {
            l.Tk = l.model.tref_k + ((V / denom) - 1.0) / l.model.thermal_alpha_per_k;
        }
    } else if (eq_process(process, "abat")) {
        l.Pl = pressure_from_volume(V, l.Tk, l.nl);
        const double pr = (prev_p > 0.0) ? l.Pl / prev_p : 1.0;
        l.Tk = l.Tk * std::pow(std::max(pr, 1.0e-9), gamma() - 1.0);
    } else if (eq_process(process, "cut")) {
        l.nl /= l.Vl / V;
    }

    l.Vl = V;
}

void liquid::dPl(double P, const char *process) {
    if (l.Pl == P) return;
    if (P < 0.0) return;

    const double prev_p = l.Pl;
    const double pr = (l.Pl > 0.0) ? P / l.Pl : 1.0;
    if (eq_process(process, "ithm")) {
        l.Vl = volume_from_pressure(P, l.Tk, l.nl);
    } else if (eq_process(process, "abat")) {
        l.Vl = volume_from_pressure(P, l.Tk, l.nl);
        const double pr_local = (prev_p > 0.0) ? P / prev_p : pr;
        l.Tk = l.Tk * std::pow(std::max(pr_local, 1.0e-9), gamma() - 1.0);
    } else if (eq_process(process, "ibar")) {
        l.Vl = volume_from_pressure(P, l.Tk, l.nl);
    }

    l.Pl = P;
}

void liquid::dnl(double n, const char *process) {
    if (l.nl == n) return;
    if (n <= 0.0) return;

    if (eq_process(process, "ithm")) {
        l.Pl = pressure_from_volume(l.Vl, l.Tk, n);
    } else if (eq_process(process, "ibar")) {
        l.Vl = volume_from_pressure(l.Pl, l.Tk, n);
    }

    l.nl = n;
}

void liquid::dTk(double T, const char *process) {
    if (l.Tk == T) return;
    if (T <= 0.0) return;

    if (eq_process(process, "ibar")) {
        l.Vl = volume_from_pressure(l.Pl, T, l.nl);
    } else if (eq_process(process, "abat")) {
        l.Vl = volume_from_pressure(l.Pl, T, l.nl);
    } else if (eq_process(process, "ivol")) {
        l.Pl = pressure_from_volume(l.Vl, T, l.nl);
    } else if (eq_process(process, "ithm")) {
        l.Vl = volume_from_pressure(l.Pl, T, l.nl);
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
    const double c = compression_fraction(l.Pl, l.Tk);
    return 1.05 + 0.55 * c;
}

void liquid::equalize_to(const ambient &outside, const char *process) {
    dPl(outside.pressure(), process);
}
