//
// Created by vladi on 9/03/2026.
//

#ifndef CLION_LIQUID_H
#define CLION_LIQUID_H

class ambient;

struct liquid_t {
    char name[64];
    double Pl;
    double Vl;
    double nl;
    double Tk;
    double mm;
    double cg[18];
};

double solve_liquid(double P, double V, double n, double T);

class liquid {
public:
    liquid_t l{};

    liquid(const char *name,
           double mm, double nl, double Vl,
           double Tk, double Pl, const double cg[18]);

    int eq_process(const char *a, const char *b) const;

    double mass() const;

    void print() const;

    void dVl(double V, const char *process);

    void dPl(double P, const char *process);

    void dTk(double T, const char *process);

    void dnl(double n, const char *process);

    void test_stepwise(double target_V, const char *process, int steps);

    double gamma();

    void equalize_to(const ambient &outside, const char *process = "ithm");
};

#endif //CLION_LIQUID_H
