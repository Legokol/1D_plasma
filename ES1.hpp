//
// Created by legokol on 16.07.2022.
//

#ifndef INC_1D_PLASMA_ES1_HPP
#define INC_1D_PLASMA_ES1_HPP

#include <vector>
#include <cmath>

struct Particle {
    const double m; // Масса частицы
    const double q; // Заряд частицы
    const double qm = q / m; // Удельный заряд
    double x; // Координата частицы в момент времени t
    double v; // Скорость частицы в момент времени t - dt/2
};

class ES1 {
private:

    // Узел расчётной сетки
    struct Cell {
        double x; // Координата узла
        double rho; // Плотность заряда
        double E; // Проекция напряжённости электрического поля на ось x
    };

    double timeStep; // Шаг по времени
    double step; // Шаг сетки
    double L; // Суммарная длина сетки

    std::vector<Particle> particles; // Массив частиц
    std::vector<Cell> grid; // Расчётная сетка

    void weighting(); // Взвешивание

    double interpolateField(const Particle &particle); // Интерполяция электрического поля, действующего на частицу

    void moveParticles(); // Движение частиц
};

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


#endif //INC_1D_PLASMA_ES1_HPP
