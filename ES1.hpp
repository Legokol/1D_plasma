//
// Created by legokol on 16.07.2022.
//

#ifndef INC_1D_PLASMA_ES1_HPP
#define INC_1D_PLASMA_ES1_HPP

#include <vector>
#include <cmath>

#include "Constants.hpp"
#include "FourierTransform.hpp"

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

    const FourierTransform fourierTransform;

    std::vector<double> frequency; // Частоты преобразования Фурье
    std::vector<double> frequency1; // Частоты преобразования Фурье, домноженные на sin(k dx) / k dx
    std::vector<double> frequency2; // Частоты преобразования Фурье, домноженные на sin(k dx / 2) / (k dx / 2)

    const double _timeStep; // Шаг по времени
    const double _step; // Шаг сетки
    const double _L; // Суммарная длина сетки

    std::vector<Particle> _particles; // Массив частиц
    std::vector<Cell> grid; // Расчётная сетка

    void weighting(); // Взвешивание

    double interpolateField(const Particle &particle); // Интерполяция электрического поля, действующего на частицу

    void moveParticles(); // Движение частиц

public:
    ES1(double timeStep, double step, double L, const std::vector<Particle> &particles);
};

#endif //INC_1D_PLASMA_ES1_HPP
