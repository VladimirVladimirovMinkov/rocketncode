#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

#include "ambient.h"
#include "gas.h"
#include "liquid.h"
#include "nozzle.h"
#include "pipe.h"
#include "tank.h"
#include "valve.h"

namespace {
constexpr double NASA_T_MIN = 200.0;
constexpr double NASA_T_MAX = 6000.0;

double clampd(const double value, const double lo, const double hi) {
    return std::max(lo, std::min(value, hi));
}

void print_nasa9_point(const gas &g, const double T) {
    const double dt = std::max(1.0e-3 * T, 1.0e-2);
    double t0 = T - dt;
    double t1 = T + dt;
    if (T < 1000.0) {
        t0 = clampd(t0, NASA_T_MIN + 1.0e-6, 999.0);
        t1 = clampd(t1, NASA_T_MIN + 1.0e-6, 999.0);
    } else {
        t0 = clampd(t0, 1000.1, NASA_T_MAX - 1.0e-6);
        t1 = clampd(t1, 1000.1, NASA_T_MAX - 1.0e-6);
    }
    if (t1 <= t0) {
        t1 = t0 + std::max(dt, 1.0e-3);
    }

    const double cp = g.cp_molar(T);
    const double h0 = g.enthalpy_molar(t0);
    const double h1 = g.enthalpy_molar(t1);
    const double cp_fd = (h1 - h0) / std::max(1.0e-12, (t1 - t0));
    const double s0 = g.entropy_molar(t0);
    const double s1 = g.entropy_molar(t1);
    const double dsdt = (s1 - s0) / std::max(1.0e-12, (t1 - t0));
    const double cp_over_t = cp / std::max(T, 1.0e-9);

    const double rel_cp = (cp_fd - cp) / std::max(cp, 1.0e-9);
    const double rel_s = (dsdt - cp_over_t) / std::max(cp_over_t, 1.0e-9);

    std::printf("  T=%.2f K | Cp=%.6g J/mol/K | dH/dT=%.6g | rel=%.3g | dS/dT=%.6g | rel=%.3g | gamma=%.6g\n",
                T, cp, cp_fd, rel_cp, dsdt, rel_s, g.gamma(T));
}

void run_nasa9_self_checks(const char *label, const gas &g) {
    std::printf("%s NASA9 self-checks:\n", label);
    print_nasa9_point(g, 300.0);
    print_nasa9_point(g, 1000.0);
    print_nasa9_point(g, 3000.0);

    const double t_lo = 999.9;
    const double t_hi = 1000.1;
    const double cp_lo = g.cp_molar(t_lo);
    const double cp_hi = g.cp_molar(t_hi);
    const double h_lo = g.enthalpy_molar(t_lo);
    const double h_hi = g.enthalpy_molar(t_hi);
    const double s_lo = g.entropy_molar(t_lo);
    const double s_hi = g.entropy_molar(t_hi);

    std::printf("  Continuity @1000K: dCp=%.6g dH=%.6g dS=%.6g\n\n",
                cp_hi - cp_lo, h_hi - h_lo, s_hi - s_lo);
}

void print_pipe_segments(const pipe &p, const char *label, int count) {
    const int segs = p.segment_count();
    if (segs <= 0) {
        std::printf("%s segments: none\n", label);
        return;
    }
    const int to_print = std::max(1, std::min(count, segs));
    std::printf("%s segments (first %d):\n", label, to_print);
    for (int i = 0; i < to_print; ++i) {
        std::printf("  seg %d:\n", i);
        p.segment(i).print();
    }
    if (segs > to_print) {
        std::printf("  seg %d (last):\n", segs - 1);
        p.segment(segs - 1).print();
    }
}

void print_flow_bar(const char *label, const double value, const double max_value) {
    const int steps = 5;
    const double denom = (max_value > 0.0) ? max_value : 1.0;
    int filled = static_cast<int>(std::round((value / denom) * steps));
    filled = std::max(0, std::min(filled, steps));
    char bar[6];
    for (int i = 0; i < steps; ++i) {
        bar[i] = (i < filled) ? '#' : '.';
    }
    bar[steps] = '\0';
    std::printf("  %-16s [%s] %.6g mol/s\n", label, bar, value);
}
}

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

    // Liquid-equivalent spline model:
    // c(P,T) = clamp(c_ref(P) * (1 + betaT*(T - Tref)), 0, 1 - dead_space(mm))
    // V(P,T,n) = n*Vm_ref*(1 + alpha*(T - Tref))*(1 - c(P,T))
    const liquid_model_t air_liquid_model = {
        87.0, 3.2e-3, 2.1e-3, 0.18, 32.0, 3.10e-5,
        {0.0, 100.0, 0.0, 0.000919596406543313, 8.5884285331187e-06, -5.80249396621072e-08}
    };

    const liquid_model_t co2_liquid_model = {
        216.6, 4.5e-3, 3.2e-3, 0.22, 44.0, 4.10e-5,
        {0.0, 100.0, 0.0, 0.00106091909565355, 1.97690657425188e-05, -1.23166656352509e-07}
    };

    const liquid_model_t n2o_liquid_model = {
        184.7, 3.9e-3, 3.0e-3, 0.24, 44.0, 3.70e-5,
        {0.0, 100.0, 0.0, 0.00109999312388692, 2.28878127771937e-05, -1.3780590108024e-07}
    };

    const liquid_model_t o2_liquid_model = {
        90.2, 1.5e-3, 2.2e-3, 0.20, 32.0, 2.80e-5,
        {0.0, 100.0, 0.0, 0.000885982527796691, 7.36573529715107e-06, -5.17178249479122e-08}
    };

    const liquid_model_t n2_liquid_model = {
        77.4, 2.8e-3, 2.5e-3, 0.19, 28.0, 3.47e-5,
        {0.0, 100.0, 0.0, 0.000941933975562325, 1.00374146502456e-05, -6.79849894451597e-08}
    };

    // Pure-component PR constants: Tc [K], Pc [Pa], acentric factor omega.
    gas air("air, dry", 28.967, 1.0, 0.022414, 273.15, 101325.0, air_cg,
            gas_eos_t::peng_robinson, 132.45, 3.77e6, 0.035);
    gas co2("CO2", 44.009, 1.0, 0.022414, 273.15, 101325.0, co2_cg,
            gas_eos_t::peng_robinson, 304.1282, 7.3773e6, 0.22394);
    gas n2o("N2O", 44.013, 1.0, 0.022414, 273.15, 101325.0, n2o_cg,
            gas_eos_t::peng_robinson, 309.57, 7.245e6, 0.162);
    gas o2("O2", 31.9988, 1, 0.022414, 273.15, 101325.0, o2_cg,
           gas_eos_t::peng_robinson, 154.58, 5.043e6, 0.0222);
    gas n2("N2", 28.0134, 1, 0.022414, 273.15, 101325.0, n2_cg,
           gas_eos_t::peng_robinson, 126.2, 3.3958e6, 0.0372);
    ambient outside("outside", 101325.0, 293.15);

    printf("\n--- NASA9 SELF-CHECKS (Cp/H/S identities + 1000K continuity) ---\n");
    run_nasa9_self_checks("Air", air);
    run_nasa9_self_checks("CO2", co2);
    run_nasa9_self_checks("N2O", n2o);
    run_nasa9_self_checks("O2", o2);
    run_nasa9_self_checks("N2", n2);

    gas gases[3] = {air, co2, n2o};

    liquid l_air("air, liquid-equivalent", 28.967, 1.0, 0.0, 87.0, 101325.0, air_liquid_model);
    liquid l_co2("CO2, liquid-equivalent", 44.009, 1.0, 0.0, 216.6, 101325.0, co2_liquid_model);
    liquid l_n2o("N2O, liquid-equivalent", 44.013, 1.0, 0.0, 184.7, 101325.0, n2o_liquid_model);
    liquid l_o2("O2, liquid-equivalent", 31.9988, 1.0, 0.0, 90.2, 101325.0, o2_liquid_model);
    liquid l_n2("N2, liquid-equivalent", 28.0134, 1.0, 0.0, 77.4, 101325.0, n2_liquid_model);

    liquid liquids[5] = {l_air, l_co2, l_n2o, l_o2, l_n2};

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

    printf("\n--- LIQUID-EQUIVALENT SPLINE COEFFICIENTS ---\n");
    for (auto &lf: liquids) {
        lf.print();
        lf.print_spline();
    }

    const double pressure_points_pa[3] = {10.0e6, 50.0e6, 100.0e6};
    printf("\n--- LIQUID-EQUIVALENT PRESSURE SWEEP (10/50/100 MPa) ---\n");
    for (auto &lf: liquids) {
        for (const double p: pressure_points_pa) {
            lf.dPl(p, "ithm");
            lf.print();
        }

        // Hold volume, raise temperature, recompute pressure for temperature sensitivity.
        lf.dTk(lf.l.Tk + 10.0, "ivol");
        lf.print();
    }

    printf("\n--- CHAIN SETUP (Tank -> Pipe -> Valve -> Pipe -> Chamber -> Nozzle) ---\n");
    const double supply_pressure_pa = 20.0e6;
    const double supply_temperature_k = 293.15;
    const double supply_moles = 1.0;
    const double pipe_length_m = 0.10;
    const double pipe_diameter_m = 0.01;
    const int pipe_segments = 5;
    const double chamber_volume_m3 = 0.002;
    const double nozzle_diameter_m = 0.008;

    gas n2o_supply("N2O supply", 44.013, supply_moles, 0.0, supply_temperature_k, supply_pressure_pa, n2o_cg);
    const double supply_volume_m3 = n2o_supply.volume_from_state(supply_pressure_pa, supply_temperature_k, supply_moles);
    tank supply_tank("N2O tank", supply_volume_m3);
    supply_tank.add_gas(n2o_supply);

    gas n2o_seed("N2O seed", 44.013, 0.0, 0.0, supply_temperature_k, 0.0, n2o_cg);
    pipe pipe_up("pipe_up", pipe_segments, pipe_length_m, pipe_diameter_m, n2o_seed);
    pipe pipe_down("pipe_down", pipe_segments, pipe_length_m, pipe_diameter_m, n2o_seed);

    valve main_valve("main_valve", pipe_diameter_m, 0.9);
    main_valve.set_opening(0.5);

    tank chamber("reaction_chamber", chamber_volume_m3);
    chamber.add_gas(n2o_seed);

    nozzle exit_nozzle("exit_nozzle", nozzle_diameter_m, 0.95);

    const double prime_moles = 0.05;
    pipe_up.set_total_moles(prime_moles, "ithm");
    supply_tank.gas_at(0).dng(supply_moles - prime_moles, "ithm");

    std::printf("tank V pipe segments V valve V pipe segments V tank V nozzle\n");
    supply_tank.print();
    std::printf("V\n");
    print_pipe_segments(pipe_up, "Upstream pipe", 2);
    std::printf("V\n");
    printf("Valve opening: %.1f%% area=%.6g m^2\n",
           main_valve.opening() * 100.0,
           main_valve.area());

    const double dt = 1.0e-3;
    const int steps = 100;
    const double omega = 2.0 * M_PI / static_cast<double>(steps - 1);
    double flow_sum = 0.0;
    double last_flow = 0.0;
    for (int i = 0; i < steps; ++i) {
        const double step = omega * static_cast<double>(i);
        main_valve.set_opening(std::sin(step));

        const valve_flow_t valve_flow = main_valve.molar_flow(
            pipe_up.segment(pipe_segments - 1).gas_at(0),
            pipe_down.segment(0).gas_at(0));
        double dn = valve_flow.mol_per_s * dt;
        const double n_up = pipe_up.segment(pipe_segments - 1).gas_at(0).g.ng;
        const double n_down = pipe_down.segment(0).gas_at(0).g.ng;
        if (dn > 0.0 && dn > n_up) {
            dn = n_up;
        } else if (dn < 0.0 && -dn > n_down) {
            dn = -n_down;
        }
        const double mol_per_s_limited = (dt > 0.0) ? (dn / dt) : 0.0;
        flow_sum += std::abs(mol_per_s_limited);
        last_flow = mol_per_s_limited;
        pipe_up.apply_molar_flow(pipe_segments - 1, -mol_per_s_limited, dt, "ithm");
        pipe_down.apply_molar_flow(0, mol_per_s_limited, dt, "ithm");
    }

    printf("After %.3g s x %d steps:\n", dt, steps);
    std::printf("Upstream pipe (last seg):\n");
    pipe_up.segment(pipe_segments - 1).print();
    std::printf("V\n");
    std::printf("Valve opening: %.1f%% area=%.6g m^2 last_flow=%.6g mol/s\n",
                main_valve.opening() * 100.0,
                main_valve.area(),
                last_flow);
    std::printf("V\n");
    print_pipe_segments(pipe_down, "Downstream pipe", 2);
    std::printf("V\n");
    std::printf("Chamber:\n");
    chamber.print();
    std::printf("V\n");
    std::printf("Nozzle area: %.6g m^2\n", exit_nozzle.area());
    const double avg_flow = flow_sum / static_cast<double>(steps);
    gas ambient_gas("ambient", 28.967, 1.0, 1.0, 293.15, 101325.0, air_cg);
    const nozzle_flow_t nozzle_flow = exit_nozzle.molar_flow(chamber.gas_at(0), ambient_gas);
    const double nozzle_mol = std::abs(nozzle_flow.mol_per_s);
    const double max_flow = std::max(avg_flow, nozzle_mol);

    std::printf("Flow rate graph (avg over %d steps):\n", steps);
    print_flow_bar("tank->pipe", avg_flow, max_flow);
    print_flow_bar("valve", avg_flow, max_flow);
    print_flow_bar("pipe->tank", avg_flow, max_flow);
    print_flow_bar("nozzle", nozzle_mol, max_flow);

    int error;
    scanf("%d", &error);
    return error;
}
