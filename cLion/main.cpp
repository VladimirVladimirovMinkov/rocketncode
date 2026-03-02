/* pressureDupe.c
 *
 * direct translation of the python script pressure.py into C.
 * struct gas_t holds the state; a collection of helper functions
 * mimic the methods in the original class.  A small test
 * harness at the bottom exercises the same sequences of changes
 * that were in the Python version.
 */

#include <cstdio>
#include <tgmath.h>
#include <cstring>

const double RCONST = 8.31446261815324; /* exact value used in python */

typedef struct {
    char name[64];
    double Pg;
    double Vg;
    double ng;
    double Tk;
    double mm;
    double cg[4];
} gas_t;

double solve(double P, double V, double n, double T, double R) {
    /* return the missing variable in PV = nRT;
       pass non‑positive value for the unknown quantity */
    if (P <= 0.0) return (n * R * T) / V;
    if (V <= 0.0) return (n * R * T) / P;
    if (n <= 0.0) return (P * V) / (R * T);
    if (T <= 0.0) return (P * V) / (n * R);
    return 0.0;
}

static int eq_process(const char *a, const char *b) {
    /* case‑insensitive comparison of short process names */
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b))
            return 0;
        a++;
        b++;
    }
    return *a == *b;
}

void gas_init(gas_t *g, const char *name,
              double mm, double ng, double Vg,
              double Tk, double Pg, const double cg[4]) {
    strncpy(g->name, name, sizeof(g->name) - 1);
    g->name[sizeof(g->name) - 1] = '\0';
    g->mm = mm;
    g->ng = ng;
    g->Vg = Vg;
    g->Tk = Tk;
    g->Pg = Pg;
    if (cg) {
        for (int i = 0; i < 4; i++)
            g->cg[i] = cg[i];
    } else {
        /* default air coefficients */
        g->cg[0] = 28.11;
        g->cg[1] = 0.1967;
        g->cg[2] = 4.802;
        g->cg[3] = -1.966;
    }
}

double gas_mass(const gas_t *g) {
    return g->ng * g->mm;
}

void gas_print(const gas_t *g) {
    printf("%s:\n", g->name);
    printf("Pressure:    Pa:  %.8g Atmo: %.8g Bar: %.8g Tor: %.8g Psi: %.8g\n",
           g->Pg, g->Pg / 101325.0, g->Pg / 100000.0,
           g->Pg / 133.322, g->Pg / 6894.76);
    printf("Moles:       n:   %.8g Kmol: %.8g\n",
           g->ng, g->ng / 1000.0);
    printf("Mass:        g:   %.8g   Kg: %.8g\n",
           gas_mass(g), gas_mass(g) / 1000.0);
    printf("Volume:      m^3: %.8g    L: %.8g\n",
           g->Vg, g->Vg * 1000.0);
    printf("Temperature: K:   %.8g   °C: %.8g  °F: %.8g\n\n",
           g->Tk, g->Tk - 273.15, (g->Tk - 273.15) * 1.8 + 32);
}

void gas_dVg(gas_t *g, double Vg, const char *process) {
    if (g->Vg == Vg) return;
    if (eq_process(process, "ithm")) {
        g->Pg = (g->ng * RCONST * g->Tk) / Vg;
    } else if (eq_process(process, "ibar")) {
        g->Tk = (g->Pg * Vg) / (g->ng * RCONST);
    } else if (eq_process(process, "abat")) {
        double t = g->Tk / 1000.0;
        double c = g->cg[0] + g->cg[1] * t + g->cg[2] * t * t
                   + g->cg[3] * t * t * t;
        double gamma = c / (c - RCONST);
        double vr = g->Vg / Vg;
        g->Pg = g->Pg * pow(vr, gamma);
        g->Tk = g->Tk * pow(vr, gamma - 1.0);
    } else if (eq_process(process, "cut")) {
        g->ng /= g->Vg / Vg;
    }
    g->Vg = Vg;
}

void gas_dPg(gas_t *g, double Pg, const char *process) {
    if (g->Pg == Pg) return;
    double pr = Pg / g->Pg;
    if (eq_process(process, "ithm")) {
        g->Vg = (g->ng * RCONST * g->Tk) / Pg;
    } else if (eq_process(process, "abat")) {
        double t = g->Tk / 1000.0;
        double cp = g->cg[0] + g->cg[1] * t + g->cg[2] * t * t
                    + g->cg[3] * t * t * t;
        double gamma = cp / (cp - RCONST);
        g->Vg = g->Vg * pow(1.0 / pr, 1.0 / gamma);
        g->Tk = g->Tk * pow(pr, (gamma - 1.0) / gamma);
    }
    g->Pg = Pg;
}

void gas_dng(gas_t *g, double ng, const char *process) {
    if (g->ng == ng) return;
    if (eq_process(process, "ithm")) {
        g->Pg = (ng * RCONST * g->Tk) / g->Vg;
    } else if (eq_process(process, "ibar")) {
        g->Vg = (ng * RCONST * g->Tk) / g->Pg;
    }
    g->ng = ng;
}

void gas_dTk(gas_t *g, double Tk, const char *process) {
    if (g->Tk == Tk) return;
    double tr = Tk / g->Tk;
    if (eq_process(process, "ibar")) {
        g->Vg = (g->ng * RCONST * Tk) / g->Pg;
    } else if (eq_process(process, "abat")) {
        double t = g->Tk / 1000.0;
        double cp = g->cg[0] + g->cg[1] * t + g->cg[2] * t * t
                    + g->cg[3] * t * t * t;
        double gamma = cp / (cp - RCONST);
        g->Vg = g->Vg * pow(1.0 / tr, 1.0 / (gamma - 1.0));
        g->Pg = g->Pg * pow(tr, gamma / (gamma - 1.0));
    }
    g->Tk = Tk;
}

void test_stepwise(gas_t *g, double target_V, const char *process, int steps) {
    double total_delta_V = target_V - g->Vg;
    double step_size = total_delta_V / steps;

    printf("--- Step-wise Test: %s | %s | %d steps ---\n",
           g->name, process, steps);
    gas_print(g);

    for (int i = 0; i < steps; i++) {
        double new_v = g->Vg + step_size;
        gas_dVg(g, new_v, process);
    }

    printf("Final state after %d steps:\n", steps);
    gas_print(g);

    double expected_P = solve(0, g->Vg, g->ng, g->Tk, RCONST);
    double error_P = g->Pg - expected_P;

    printf("Ideal Gas Law Check (solve):\n");
    printf("  Calculated P: %10.2f Pa\n", g->Pg);
    printf("  Expected P:   %10.2f Pa\n", expected_P);
    printf("  Accuracy Error: %10.6f Pa\n", error_P);
    printf("--------------------------------------------------\n");
}

int main(void) {
    gas_t air, co2, n2o;
    double air_cg[4] = {28.11, 0.1967, 4.802, -1.966};
    double co2_cg[4] = {24.99, 55.18, -33.69, 7.948};
    double n2o_cg[4] = {26.85, 50.81, -31.42, 7.76};

    gas_init(&air, "air, dry", 28.967, 1.0, 0.022414, 273.15,
             101325.0, air_cg);
    gas_init(&co2, "CO2", 44.009, 1.0, 0.022414, 273.15,
             101325.0, co2_cg);
    gas_init(&n2o, "N2O", 44.013, 1.0, 0.022414, 273.15,
             101325.0, n2o_cg);

    gas_t *gases[3] = {&air, &co2, &n2o};

    printf("--- TESTING dVg (Volume) | Process: abat ---\n");
    for (int i = 0; i < 3; i++) {
        gas_t *g = gases[i];
        gas_print(g);
        gas_dVg(g, g->Vg / 2.0, "abat");
        gas_print(g);
        gas_dVg(g, g->Vg * 2.0, "abat");
        gas_print(g);
        printf("------------------------------\n");
    }

    printf("\n--- TESTING dPg (Pressure) | Process: abat ---\n");
    for (int i = 0; i < 3; i++) {
        gas_t *g = gases[i];
        gas_print(g);
        gas_dPg(g, g->Pg / 2.0, "abat");
        gas_print(g);
        gas_dPg(g, g->Pg * 2.0, "abat");
        gas_print(g);
        printf("------------------------------\n");
    }

    printf("\n--- TESTING dng (Moles) | Process: ithm ---\n");
    for (int i = 0; i < 3; i++) {
        gas_t *g = gases[i];
        gas_print(g);
        gas_dng(g, g->ng / 2.0, "ithm");
        gas_print(g);
        gas_dng(g, g->ng * 2.0, "ithm");
        gas_print(g);
        printf("------------------------------\n");
    }

    printf("\n--- TESTING dTk (Temperature) | Process: ibar ---\n");
    for (int i = 0; i < 3; i++) {
        gas_t *g = gases[i];
        gas_print(g);
        gas_dTk(g, g->Tk / 2.0, "ibar");
        gas_print(g);
        gas_dTk(g, g->Tk * 2.0, "ibar");
        gas_print(g);
        printf("------------------------------\n");
    }

    printf("\n--- TESTING THE 'CUT' (V/2, n/2) ---\n");
    for (int i = 0; i < 3; i++) {
        gas_t *g = gases[i];
        gas_print(g);
        gas_dVg(g, g->Vg / 2.0, "cut");
        gas_print(g);
        gas_dVg(g, g->Vg * 2.0, "cut");
        gas_print(g);
        printf("------------------------------\n");
    }

    gas_t air_step;
    gas_init(&air_step, "Air Step Test", 28.967, 1.0, 0.022414,
             273.15, 101325.0, air_cg);
    test_stepwise(&air_step, air_step.Vg / 2.0, "abat", 1);

    gas_t air_precise;
    gas_init(&air_precise, "Air High Precision", 28.967, 1.0, 0.022414,
             273.15, 101325.0, air_cg);
    test_stepwise(&air_precise, air_precise.Vg / 2.0, "abat", 100);

    return 0;
}
