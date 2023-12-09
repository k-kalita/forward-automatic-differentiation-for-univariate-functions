#include "../include/Jet.h"
#include <cstdio>
#include <stack>
#include <vector>
#include <stdexcept>

Jet::StackElem::StackElem(int f_der, int g_der, int order, double value) {
    this->f_der = f_der;
    this->g_der = g_der;
    this->order = order;
    this->value = value;
}

Jet::Jet(const Jet& other) {
    this->der_order = other.der_order;
    this->der = new double[other.der_order + 1];

    for (int i = 0; i <= der_order; i++)
        der[i] = other.der[i];
}

Jet::Jet(int der_order) {
    this->der_order = der_order;
    this->der = new double[der_order + 1];

    for (int i = 0; i <= der_order; i++)
        der[i] = 0;
}

Jet::Jet(int der_order, double value) {
    this->der_order = der_order;
    this->der = new double[der_order + 1];

    for (int i = 2; i <= der_order; i++)
        der[i] = 0;
    der[0] = value;
    der[1] = 1;
}

Jet Jet::operator/(const Jet& other) const {
    Jet result(der_order);
    Jet f = *this;
    Jet g = other;
    for (int i = 0; i <= der_order && f.der[0] != 0; i++) {
        result.der[i] = f.der[0] / g.der[0];
        f = f.diff() * g - f * g.diff();
        g = g * g;
    }
    return result;
}

Jet Jet::operator/(const double factor) const {
    Jet result(der_order);
    for (int i = 0; i <= der_order; i++)
        result.der[i] = der[i] / factor;
    return result;
}

Jet Jet::operator*(const Jet& other) const {
    Jet result(der_order);
    std::stack<StackElem> stack;
    StackElem temp = StackElem(0, 0, 0, der[0] * other.der[0]);
    stack.emplace(temp);

    while (!stack.empty()) {
        temp = stack.top();
        stack.pop();
        result.der[temp.order] += temp.value;
        if (temp.order < der_order - 1) {
            stack.emplace(
                    temp.f_der + 1,
                    temp.g_der,
                    temp.order + 1,
                    der[temp.f_der + 1] * other.der[temp.g_der]
            );
            stack.emplace(
                    temp.f_der,
                    temp.g_der + 1,
                    temp.order + 1,
                    der[temp.f_der] * other.der[temp.g_der + 1]
            );
        }
    }

    return result;
}

Jet Jet::operator*(const double factor) const {
    Jet result(der_order);
    for (int i = 0; i <= der_order; ++i)
        result.der[i] = der[i] * factor;
    return result;
}

Jet Jet::operator+(const Jet& other) const {
    Jet result(der_order);
    for (int i = 0; i <= der_order; i++)
        result.der[i] = der[i] + other.der[i];
    return result;
}

Jet Jet::operator+(const double other) const {
    Jet result(*this);
    result.der[0] += other;
    return result;
}

Jet Jet::operator-(const Jet& other) const {
    Jet result(der_order);
    for (int i = 0; i <= der_order; i++)
        result.der[i] = der[i] - other.der[i];
    return result;
}

Jet Jet::operator-(const double other) const {
    Jet result(*this);
    result.der[0] -= other;
    return result;
}

Jet Jet::operator-() const {
    Jet result(der_order);
    for (int i = 0; i <= der_order; i++)
        result.der[i] = -der[i];
    return result;
}

Jet& Jet::operator=(const Jet& other) {
    if (this == &other)
        return *this;
    for (int i = 0; i <= der_order; i++)
        der[i] = other.der[i];
    return *this;
}

[[nodiscard]] Jet Jet::pow(int power) const {
    Jet result(*this);
    if (power > 0)
        for (int i = 0; i < power - 1; i++)
            result = result * *this;
    else if (power < 0) {
        result = result.pow(-power);
        result = 1 / result;
    }
    else {
        result = Jet(der_order);
        result.der[0] = 1;
    }
    return result;
}

[[nodiscard]] Jet Jet::diff() const {
    Jet result(der_order);
    for (int i = 0; i < der_order; i++)
        result.der[i] = der[i + 1];
    result.der[der_order - 1] = 0;
    return result;
}

void Jet::print() {
    for (int i = 0; i <= der_order; i++)
        printf("%lf ", der[i]);
    printf("\n");
}

std::vector<double> Jet::getDerivatives() {
    std::vector<double> result;
    for (int i = 0; i <= der_order && der[i] != 0.0; i++)
        result.push_back(der[i]);
    return result;
}

double Jet::getDerivative(int order) {
    if (order < 0)
        throw std::invalid_argument("Derivative order cannot be negative.");
    if (order > der_order)
        throw std::invalid_argument("Derivative order is too high.");
    return der[order];
}

double Jet::val() {
    return der[0];
}

Jet operator/(double factor, const Jet& jet) {
    Jet result(jet.der_order);
    result.der[0] = factor;
    return result / jet;
}

Jet operator*(double factor, const Jet& jet) {
    Jet result(jet.der_order);
    for (int i = 0; i < jet.der_order; ++i)
        result.der[i] = jet.der[i] * factor;
    return result;
}

Jet operator+(const double other, const Jet& jet) {
    Jet result(jet);
    result.der[0] += other;
    return result;
}

Jet operator-(const double other, const Jet& jet) {
    Jet result(-jet);
    result.der[0] = other + result.der[0];
    return result;
}