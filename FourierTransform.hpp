//
// Created by legokol on 17.07.2022.
//

#ifndef INC_1D_PLASMA_FOURIERTRANSFORM_HPP
#define INC_1D_PLASMA_FOURIERTRANSFORM_HPP

#include <complex>
#include <vector>

using complexd = std::complex<double>;

class FourierTransform {
private:
    const double _step; // Шаг сетки
    const double _L; // Период сетки

public:
    FourierTransform(double step, double L) : _step(step), _L(L) {}

    std::vector<double> frequency() const; // Метод для получения пространственных частот преобразования Фурье

    std::vector<complexd> complexFrequency() const; // Метод для получения пространственных частот в виде i * k

    std::vector<complexd> transform(const std::vector<double> &preimage) const; // Метод для прямого преобразования Фурье

    std::vector<complexd> inverse(const std::vector<complexd> &image) const;
};


#endif //INC_1D_PLASMA_FOURIERTRANSFORM_HPP
