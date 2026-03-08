//
// Created by vladi on 9/03/2026.
//

#ifndef CLION_SOLID_H
#define CLION_SOLID_H

class ambient;

struct solid_t {
    char name[64];
    double Ps;
    double Vs;
    double ns;
    double Tk;
    double mm;
    double cg[18];
};

double solve_solid(double P, double V, double n, double T);

class solid {
public:
    solid_t s{};

    solid(const char *name,
          double mm, double ns, double Vs,
          double Tk, double Ps, const double cg[18]);

    int eq_process(const char *a, const char *b) const;

    double mass() const;

    void print() const;

    void dVs(double V, const char *process);

    void dPs(double P, const char *process);

    void dTk(double T, const char *process);

    void dns(double n, const char *process);

    void test_stepwise(double target_V, const char *process, int steps);

    double gamma();

    void equalize_to(const ambient &outside, const char *process = "ithm");
};

#endif //CLION_SOLID_H
