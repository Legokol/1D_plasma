//
// Created by legokol on 16.07.2022.
//

#include "ES1.hpp"

ES1::ES1(double timeStep, double step, double L, const std::vector<Particle> &particles) : _timeStep(timeStep),
                                                                                           _step(step), _L(L),
                                                                                           fourierTransform(step, L) {
    int n = static_cast<int>(_L / _step + 1);
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
            frequency1[i] = frequency[i] * std::sin(frequency[i] * _step) / (frequency[i] * _step);
            frequency2[i] = frequency[i] * std::sin(frequency[i] * _step / 2) / (frequency[i] * _step / 2);
        }
    }

    // Расчёт плотности заряда и электрического поля
    weighting();
    // Пересчёт скорости из момента t = 0 в t = - dt / 2
    for (int i = 0; i < _particles.size(); ++i) {
        _particles[i].v -= _particles[i].qm * interpolateField(_particles[i]) * _timeStep / 2;
    }
}

void ES1::weighting() {
    std::vector<double> rho(grid.size() - 1);
    // Расчёт плотности заряда
    grid[0].rho = 0;
    for (int i = 0; i < grid.size() - 1; ++i) {
        grid[i + 1].rho = 0;
        for (int j = 0; j < _particles.size() - 1; ++j) {
            if (_particles[j].x - grid[i].x < _step && _particles[j].x > grid[i].x) {
                grid[i].rho += _particles[j].q * (grid[i + 1].x - _particles[i].x) / _step;
                grid[i + 1].rho += _particles[j].q * (_particles[i].x - grid[i].x) / _step;
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
    // TODO: Можно убрать вектор phiImage
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
            return (grid[i].E * (grid[i + 1].x - particle.x) + grid[i + 1].E * (particle.x - grid[i].x)) / _step;
        }
    }
}

void ES1::moveParticles() {
    for (int i = 0; i < _particles.size(); ++i) {
        // Обновление координаты
        _particles[i].x += _particles[i].v * _timeStep;
        // Обновление скорости
        _particles[i].v += _particles[i].qm * interpolateField(_particles[i]) * _timeStep;
        // Перенос частиц с учётом периодичности
        if (_particles[i].x > _L) {
            _particles[i].x -= _L;
        }
        if (_particles[i].x < 0) {
            _particles[i].x += _L;
        }
    }
}

void ES1::saveGrid(double time, std::ofstream &writer) const {
    writer << time << ',';
    for (int i = 0; i < grid.size(); ++i) {
        writer << grid[i].rho << ',' << grid[i].phi << ',' << grid[i].E << ',';
    }
    writer << std::endl;
}