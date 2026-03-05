#include <cstdio>
#include <cstring>
#include <ctgmath>

#include "gas.h"
#include "ambient.h"

int main() {
    // Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4
    // Indices 0-6: 200K-1000K | Indices 9-15: 1000K-6000K
    // 7,8,16,17 are used for reactions
    double o2_cg[18] = {
        -3.4255e4, 4.847e2, 1.119e0, 4.291e-3, -6.836e-7, -2.028e-9, 1.031e-12, -9.281e3, 1.113e1, // Low
        -1.0339e6, 2.296e3, -3.539e-1, 1.745e-3, -3.268e-7, 3.923e-11, -1.967e-15, -4.453e4, 4.490e0 // High
    };

    double n2_cg[18] = {
        2.2103e4, -3.818e2, 6.082e0, -8.530e-3, 1.384e-5, -9.625e-9, 2.519e-12, 7.108e2, -1.076e1, // Low
        5.8771e5, -2.239e3, 6.066e0, -6.139e-4, 1.491e-7, -1.923e-11, 1.061e-15, 1.283e4, -1.586e1 // High
    };

    double co2_cg[18] = {
        4.9436e4, -6.264e2, 5.301e0, 2.503e-3, -2.127e-7, -7.689e-10, 2.849e-13, -4.830e4, -5.700e-1, // Low
        1.1769e5, -1.788e3, 8.291e0, -9.223e-5, 4.863e-9, -1.891e-12, 6.330e-16, -3.908e4, -2.652e1 // High
    };

    double n2o_cg[18] = {
        3.5212e4, -4.621e2, 4.801e0, 6.102e-3, -4.102e-7, -3.125e-10, 5.121e-14, 8.210e3, 5.210e0, // Low
        1.2141e5, -1.892e3, 9.102e0, -1.021e-4, 5.210e-9, -2.012e-12, 7.112e-16, 9.812e3, -2.121e1 // High
    };

    double air_cg[18] = {
        -3.515e3, 1.580e2, 3.121e0, 8.112e-4, 9.825e-8, -2.819e-10, 1.103e-13, -1.012e3, 1.102e1, // Low
        1.121e5, -1.213e3, 4.212e0, 2.112e-4, -1.212e-8, 3.112e-13, -1.112e-17, -3.212e3, 2.121e0 // High
    };


    gas air("air, dry", 28.967, 1.0, 0.022414, 273.15, 101325.0, air_cg);
    gas co2("CO2", 44.009, 1.0, 0.022414, 273.15,101325.0, co2_cg);
    gas n2o("N2O", 44.013, 1.0, 0.022414, 273.15,101325.0, n2o_cg);
    gas o2("O2", 31.9988, 1, 0.022414, 273.15, 101325.0, o2_cg);
    gas n2("N2", 28.0134, 1, 0.022414, 273.15, 101325.0, n2_cg);
    ambient outside("outside", 101325.0, 293.15);

    gas gases[3] = {air, co2, n2o};

    printf("\n--- TESTING initial values ---\n");
    for (auto g: gases) {
        g.print();
    }

    printf("\n--- TESTING dVg (Pressure) | Process: abat ---\n");
    for (auto g: gases) {
        g.dVg(g.g.Vg / 2.0, "abat");
        g.dVg(g.g.Vg * 2.0, "abat");
        g.print();
    }

    printf("\n--- TESTING dPg (Pressure) | Process: abat ---\n");
    for (auto g: gases) {
        g.dPg(g.g.Pg * 2.0, "abat");
        g.print();
    }

    printf("\n--- TESTING dng (Moles) | Process: ithm ---\n");
    for (auto g: gases) {
        g.dng(g.g.ng / 2.0, "ithm");
        g.print();
    }

    printf("\n--- TESTING dTk (Temperature) | Process: ibar ---\n");
    for (auto g: gases) {
        g.dTk(g.g.Tk * 2.0, "ibar");
        g.print();
    }

    printf("\n--- TESTING THE 'CUT' (V/2, n/2) ---\n");
    for (auto g: gases) {
        g.dVg(g.g.Vg / 2.0, "cut");
        g.print();
    }

    printf("\n--- TESTING 60 STEP ABAT ---\n");
    for (auto g: gases) {
        g.test_stepwise(g.g.Vg * 2, "abat", 60);
        g.test_stepwise(g.g.Vg / 2, "abat", 60);
        g.dTk(273.15, "ivol");
        g.print();
    }

    printf("\n--- TESTING AMBIENT EQUALIZATION ---\n");
    outside.print();
    co2.dPg(202650.0, "ithm");
    co2.print();
    outside.equalize_pressure(co2, "ithm");
    co2.print();


    int error;
    scanf("%d", &error);
    return error;
}
