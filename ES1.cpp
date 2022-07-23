//
// Created by legokol on 16.07.2022.
//

#include "ES1.hpp"

ES1::ES1(double timeStep, double step, double L, const std::vector<Particle> &particles) : timeStep_(timeStep),
                                                                                           step_(step), L_(L),
                                                                                           fourierTransform(step, L) {
    int n = static_cast<int>(L_ / step_ + 1);
    grid.resize(n);
    for (int i = 0; i < grid.size(); ++i) {
        grid[i] = Cell{i * step, 0, 0};
    }
    // Расчёт частот
    frequency = fourierTransform.frequency();
    frequency1.resize(frequency.size());
    frequency2.resize(frequency.size());
    for (int i = 0; i < frequency.size(); ++i) {
        if (frequency[i] == 0) {
            frequency1[i] = frequency[i];
            frequency2[i] = frequency[i];
        } else {
            frequency1[i] = frequency[i] * std::sin(frequency[i] * step_) / (frequency[i] * step_);
            frequency2[i] = frequency[i] * std::sin(frequency[i] * step_ / 2) / (frequency[i] * step_ / 2);
        }
    }

    // Расчёт плотности заряда и электрического поля
    weighting();
    // Пересчёт скорости из момента t = 0 в t = - dt / 2
    for (int i = 0; i < particles_.size(); ++i) {
        particles_[i].v -= particles_[i].qm * interpolateField(particles_[i]) * timeStep_ / 2;
    }
}

void ES1::weighting() {
    std::vector<double> rho(grid.size() - 1);
    // Расчёт плотности заряда
    grid[0].rho = 0;
    for (int i = 0; i < grid.size() - 1; ++i) {
        grid[i + 1].rho = 0;
        for (int j = 0; j < particles_.size() - 1; ++j) {
            if (particles_[j].x - grid[i].x < step_ && particles_[j].x > grid[i].x) {
                grid[i].rho += particles_[j].q * (grid[i + 1].x - particles_[i].x) / step_;
                grid[i + 1].rho += particles_[j].q * (particles_[i].x - grid[i].x) / step_;
            }
        }
        rho[i] = grid[i].rho;
    }
    grid.back().rho += grid[0].rho;
    grid[0].rho = grid.back().rho;
    rho[0] = grid[0].rho;
    // Расчёт потенциала и электрического поля
    std::vector<complexd> rhoImage = fourierTransform.transform(rho);
    std::vector<complexd> phiImage(rhoImage.size());
    std::vector<complexd> EImage(rhoImage.size());
    for (int i = 0; i < rhoImage.size(); ++i) {
        phiImage[i] = rhoImage[i] / (epsilon0 * frequency2[i] * frequency2[i]);
        EImage[i] = complexd(0, -frequency1[i]) * phiImage[i];
    }
    std::vector<complexd> phi = fourierTransform.inverse(phiImage);
    std::vector<complexd> E = fourierTransform.inverse(EImage);
    for (int i = 0; i < grid.size() - 1; ++i) {
        grid[i].phi = phi[i].real();
        grid[i].E = E[i].real();
    }
    grid.back().phi = grid[0].phi;
    grid.back().E = grid[0].E;
}

double ES1::interpolateField(const Particle &particle) {
    for (int i = 0; i < grid.size() - 1; ++i) {
        if (grid[i].x <= particle.x && grid[i + 1].x >= particle.x) {
            return (grid[i].E * (grid[i + 1].x - particle.x) + grid[i + 1].E * (particle.x - grid[i].x)) / step_;
        }
    }
}

void ES1::moveParticles() {
    for (int i = 0; i < particles_.size(); ++i) {
        // Обновление координаты
        particles_[i].x += particles_[i].v * timeStep_;
        // Обновление скорости
        particles_[i].v += particles_[i].qm * interpolateField(particles_[i]) * timeStep_;
        // Перенос частиц с учётом периодичности
        if (particles_[i].x > L_) {
            particles_[i].x -= L_;
        }
        if (particles_[i].x < 0) {
            particles_[i].x += L_;
        }
    }
}

void ES1::saveParticles(double time) const {
    std::ofstream particleWriter(std::to_string(time) + "/particles.txt");
//    particleWriter << time << '\n';
    for (int i = 0; i < particles_.size(); ++i) {
        particleWriter << particles_[i].m << ' ' << particles_[i].q << ' ' << particles_[i].x << ' ' << particles_[i].v
                       << '\n';
    }
}

void ES1::saveGrid(double time) const {
    std::ofstream gridWriter(std::to_string(time) + "/grid.txt");
//    gridWriter << time << '\n';
    for (int i = 0; i < grid.size(); ++i) {
        gridWriter << grid[i].x << ' ' << grid[i].rho << ' ' << grid[i].phi << ' ' << grid[i].E << '\n';
    }
    gridWriter.close();
}

void ES1::calc(double maxTime, int writeInterval) {
    double time = 0;
    std::filesystem::create_directory("0");
    saveParticles(time);
    saveGrid(time);
    int n = 0;
    time += timeStep_;
    n++;
    while (time < maxTime) {
        moveParticles();
        weighting();
        if (n == writeInterval) {
            saveParticles(time);
            saveGrid(time);
            n = 0;
        }
        time += timeStep_;
        n++;
    }
}
