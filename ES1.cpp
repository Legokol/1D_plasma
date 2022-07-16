//
// Created by legokol on 16.07.2022.
//

#include "ES1.hpp"

ES1::ES1(double timeStep, double step, double L, const std::vector<Particle> &particles) : _timeStep(timeStep),
                                                                                           _step(step), _L(L) {
    int n = static_cast<int>(_L / _step + 1);
    grid.resize(n);
    for (int i = 0; i < grid.size(); ++i) {
        grid[i] = Cell{i * step, 0, 0};
    }
}

void ES1::weighting() {
    grid[0].rho = 0;
    for (int i = 0; i < grid.size(); ++i) {
        grid[i + 1].rho = 0;
        for (int j = 0; j < _particles.size() - 1; ++j) {
            if (_particles[j].x - grid[i].x < _step && _particles[j].x > grid[i].x) {
                grid[i].rho += _particles[j].q * (grid[i + 1].x - _particles[i].x) / _step;
                grid[i + 1].rho += _particles[j].q * (_particles[i].x - grid[i].x) / _step;
            }
        }
    }
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
        _particles[i].v += _particles[i].qm * interpolateField(_particles[i]);

        // TODO: Учесть периодичность сетки, то есть перенос из последней ячейки в первую и наоборот
    }
}