//
// Created by vladi on 3/03/2026.
//

#ifndef CLION_GAS_H
#define CLION_GAS_H

inline constexpr double RCONST = 8.31446261815324;
inline constexpr double KTC = 273.15;

class ambient;

struct gas_t {
    char name[64];
    double Pg;
    double Vg;
    double ng;
    double Tk;
    double mm;
    double cg[18];
};

double solve(double P, double V, double n, double T);

class gas {
public:
    gas_t g{};

    gas(const char *name,
        double mm, double ng, double Vg,
        double Tk, double Pg, const double cg[18]);

    int eq_process(const char *a, const char *b) const;

    double mass() const;

    void print() const;

    void dVg(double V, const char *process);

    void dPg(double P, const char *process);

    void dTk(double T, const char *process);

    void dng(double n, const char *process);

    void test_stepwise(double target_V, const char *process, int steps);

    double gamma();

    void equalize_to(const ambient &outside, const char *process = "ithm");
};

#endif //CLION_GAS_H
