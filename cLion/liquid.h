//
// Created by vladi on 9/03/2026.
//

#ifndef CLION_LIQUID_H
#define CLION_LIQUID_H

class ambient;

struct liquid_model_t {
    double tref_k;
    double thermal_alpha_per_k;
    double compress_temp_per_k;
    double dead_space_fraction_ref;
    double dead_space_mm_ref;
    double vm_ref_m3_per_mol;
    // Single cubic spline over [p0,p1] MPa:
    // c_ref(P) = a + b*x + c*x^2 + d*x^3, x = P_MPa - p0
    // cl[0]=p0, cl[1]=p1, cl[2]=a, cl[3]=b, cl[4]=c, cl[5]=d
    double cl[6];
};

struct liquid_t {
    char name[64];
    double Pl;
    double Vl;
    double nl;
    double Tk;
    double mm;
    liquid_model_t model;
};

double solve_liquid(double P, double V, double n, double T);

class liquid {
public:
    liquid_t l{};

    liquid();

    liquid(const char *name,
           double mm, double nl, double Vl,
           double Tk, double Pl, const liquid_model_t &model);

    int eq_process(const char *a, const char *b) const;

    double mass() const;

    double compression_fraction(double pressure_pa, double temperature_k) const;

    double volume_from_pressure(double pressure_pa, double temperature_k, double moles) const;

    double pressure_from_volume(double volume_m3, double temperature_k, double moles) const;

    void print() const;

    void print_spline() const;

    void dVl(double V, const char *process);

    void dPl(double P, const char *process);

    void dTk(double T, const char *process);

    void dnl(double n, const char *process);

    void test_stepwise(double target_V, const char *process, int steps);

    double gamma();

    void equalize_to(const ambient &outside, const char *process = "ithm");
};

#endif //CLION_LIQUID_H
