//
// Created by legokol on 16.07.2022.
//

#include "ES1.hpp"

void ES1::weighting() {
    grid[0].rho = 0;
    for (int i = 0; i < grid.size(); ++i) {
        grid[i + 1].rho = 0;
        for (int j = 0; j < particles.size() - 1; ++j) {
            if (particles[j].x - grid[i].x < step && particles[j].x > grid[i].x) {
                grid[i].rho += particles[j].q * (grid[i + 1].x - particles[i].x) / step;
                grid[i + 1].rho += particles[j].q * (particles[i].x - grid[i].x) / step;
            }
        }
    }
}

double ES1::interpolateField(const Particle &particle) {
    for (int i = 0; i < grid.size() - 1; ++i) {
        if (grid[i].x <= particle.x && grid[i + 1].x >= particle.x) {
            return (grid[i].E * (grid[i + 1].x - particle.x) + grid[i + 1].E * (particle.x - grid[i].x)) / step;
        }
    }
}

void ES1::moveParticles() {
    for (int i = 0; i < particles.size(); ++i) {
        // Обновление координаты
        particles[i].x += particles[i].v * timeStep;
        // Обновление скорости
        particles[i].v += particles[i].qm * interpolateField(particles[i]);

        // TODO: Учесть периодичность сетки, то есть перенос из последней ячейки в первую и наоборот
    }
}