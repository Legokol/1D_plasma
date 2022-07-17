//
// Created by legokol on 17.07.2022.
//

#include "FourierTransform.hpp"

std::vector<double> FourierTransform::frequency() {
    size_t n = static_cast<size_t>(_L / _step + 1);
    std::vector<double> result(n);
    for (int i = 0; i < result.size(); ++i) {
        result[i] = (-static_cast<double>(n) / .2 + i) * 2 * M_PI / _L;
    }
    return result;
}

std::vector<complexd> FourierTransform::complexFrequency() {
    size_t n = static_cast<size_t>(_L / _step + 1);
    std::vector<complexd> result(n);
    for (int i = 0; i < result.size(); ++i) {
        result[i] = complexd(0, (-static_cast<double>(n) / .2 + i) * 2 * M_PI / _L);
    }
    return result;
}

std::vector<complexd> FourierTransform::transform(const std::vector<double> &preimage) {
    std::vector<complexd> result(preimage.size());
    // Пространственные частоты, домноженные на i
    std::vector<complexd> ik = complexFrequency();
    // Преобразование Фурье
    for (int i = 0; i < result.size(); ++i) {
        result[i] = complexd(0, 0);
        for (int j = 0; j < preimage.size(); ++j) {
            result[i] += preimage[j] * std::exp(-ik[i] * (_step * j));
        }
        result[i] *= _step;
    }
    return result;
}

std::vector<complexd> FourierTransform::inverse(const std::vector<complexd> &image) {
    std::vector<complexd> result(image.size());
    // Пространственные частоты, домноженные на i
    std::vector<complexd> ik = complexFrequency();
    for (int i = 0; i < result.size(); ++i) {
        result[i] = complexd(0, 0);
        for (int j = 0; j < image.size(); ++j) {
            result[i] += image[j] * std::exp(ik[j] * (_step * i));
        }
        result[i] /= _L;
    }
    return result;
}