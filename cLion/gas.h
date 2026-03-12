//
// Created by vladi on 3/03/2026.
//

#ifndef CLION_GAS_H
#define CLION_GAS_H

inline constexpr double RCONST = 8.31446261815324;
inline constexpr double KTC = 273.15;

class ambient;

enum class gas_eos_t {
    ideal = 0,
    peng_robinson = 1
};

struct gas_t {
    char name[64];
    double Pg;
    double Vg;
    double ng;
    double Tk;
    double mm;
    gas_eos_t eos;
    double tc;
    double pc;
    double omega;
    double cg[18]; // NASA9 coefficinets for both above and below 1000K for pressure and chemical response
};

double solve(double P, double V, double n, double T);

class gas {
public:
    gas_t g{};

    gas();

    gas(const char *name,
        double mm, double ng, double Vg,
        double Tk, double Pg, const double cg[18]);

    gas(const char *name,
        double mm, double ng, double Vg,
        double Tk, double Pg, const double cg[18],
        gas_eos_t eos, double tc, double pc, double omega);

    int eq_process(const char *a, const char *b) const;

    double mass() const;

    double cp_over_R() const;
    double cp_over_R(double temperature_k) const;
    double h_over_RT(double temperature_k) const;
    double s_over_R(double temperature_k) const;
    double cp_molar(double temperature_k) const;
    double cv_molar(double temperature_k) const;
    double enthalpy_molar(double temperature_k) const;
    double entropy_molar(double temperature_k) const;

    double pressure_from_state(double volume_m3, double temperature_k, double moles) const;

    double volume_from_state(double pressure_pa, double temperature_k, double moles) const;

    double temperature_from_state(double pressure_pa, double volume_m3, double moles) const;

    void print() const;

    void dVg(double V, const char *process);

    void dPg(double P, const char *process);

    void dTk(double T, const char *process);

    void dng(double n, const char *process);

    void test_stepwise(double target_V, const char *process, int steps);

    double gamma();
    double gamma(double temperature_k) const;

    void equalize_to(const ambient &outside, const char *process = "ithm");
};

#endif //CLION_GAS_H
