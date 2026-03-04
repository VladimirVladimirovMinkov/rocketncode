//
// Created by Codex on 4/03/2026.
//

#ifndef CLION_AMBIENT_H
#define CLION_AMBIENT_H

class gas;

class ambient {
public:
    ambient(const char *name, double pressure_pa, double temperature_k);

    const char *name() const;

    double pressure() const;

    double temperature() const;

    void set_temperature(double temperature_k);

    void print() const;

    void equalize_pressure(gas &sample, const char *process = "ithm") const;

private:
    char name_[64];
    double pressure_pa_;
    double temperature_k_;
};

#endif //CLION_AMBIENT_H
