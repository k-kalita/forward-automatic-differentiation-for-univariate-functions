#ifndef RECURSIVE_FAD_JET_H
#define RECURSIVE_FAD_JET_H

#include <vector>

/**
 * @brief Jet is a class that represents the value of a single-variable function and its derivatives
 * up to a certain order at a given point. It can be used to calculate the value of any rational
 * function and its derivatives up to a desired order at a given point using forward automatic
 * differentiation.
 *
 * @param der_order is the order of the highest derivative that will be calculated.
 * @param der is an array of doubles that stores under index i the value of the i-th derivative of
 * the function at the given point.
 */
class Jet {

private:

    /**
     * @brief StackElem is a struct used to temporarily store values of derivatives during Jet
     * multiplication.
     *
     * @param f_der is the order of the derivative of the first function in the product at the
     * current step of the multiplication.
     * @param g_der is the order of the derivative of the second function in the product at the
     * current step of the multiplication.
     * @param order is the order of the derivative of the product at the current step of the
     * multiplication.
     * @param value is the value of the derivative of the product at the current step of the
     * multiplication.
     */
    struct StackElem {
        int f_der;
        int g_der;
        int order;
        double value;

        StackElem(int f_der, int g_der, int order, double value);
    };

    int der_order;
    double* der;

public:

    // ------------------------------------- Constructors ------------------------------------- //

    /**
     * Copy constructor for Jet.
     *
     * @param other is the Jet object to be copied.
     */
    Jet(const Jet& other);

    /**
     * Constructs an empty Jet object
     *
     * @param der_order is the order of the highest derivative that will be calculated.
     */
    explicit Jet(int der_order);

    /**
     * Constructs a Jet object with a given value.
     *
     * @param der_order is the order of the highest derivative that will be calculated.
     * @param value is the value of the function at the given point.
     */
    explicit Jet(int der_order, double value);

    // --------------------------------------- Operators -------------------------------------- //

    /**
     * Calculates the Jet resulting from a division by another Jet.
     *
     * @param other is the Jet object to divide by.
     * @return the Jet resulting from the division.
     */
    Jet operator/(const Jet& other) const;

    /**
     * Calculates the Jet resulting from a division by a double.
     *
     * @param factor is the double to divide by.
     * @return the Jet resulting from the division.
     */
    Jet operator/(double factor) const;
    friend Jet operator/(double factor, const Jet& jet);

    /**
     * Calculates the Jet resulting from a multiplication by another Jet. Uses a stack to
     * iteratively calculate the derivatives of the product.
     *
     * @param other is the Jet object to multiply by.
     * @return the Jet resulting from the multiplication.
     */
    Jet operator*(const Jet& other) const;

    /**
     * Calculates a Jet resulting from multiplication of a Jet by a double.
     *
     * @param jet is the Jet object to multiply.
     * @param factor is the double to multiply by.
     * @return the Jet resulting from the multiplication.
     */
    Jet operator*(double factor) const;
    friend Jet operator*(double factor, const Jet& jet);

    /**
     * Calculates a Jet resulting from addition of another Jet.
     *
     * @param other is the Jet object to add.
     * @return the Jet resulting from the addition.
     */
    Jet operator+(const Jet& other) const;

    /**
     * Calculates a Jet resulting from addition of a double.
     *
     * @param other is the double to add.
     * @return the Jet resulting from the addition.
     */
    Jet operator+(double other) const;
    friend Jet operator+(double other, const Jet& jet);

    /**
     * Calculates a Jest resulting from subtraction of another Jet.
     *
     * @param other is the Jet object to subtract.
     * @return the Jet resulting from the subtraction.
     */
    Jet operator-(const Jet& other) const;

    /**
     * Calculates a Jest resulting from subtraction of a double.
     *
     * @param other is the double to subtract.
     * @return the Jet resulting from the subtraction.
     */
    Jet operator-(double other) const;
    friend Jet operator-(double other, const Jet& jet);

    /**
     * Calculates a Jet resulting from the negation of this Jet.
     *
     * @return the Jet resulting from the negation.
     */
    Jet operator-() const;

    /**
     * Assigns the value of another Jet to this Jet.
     *
     * @param other is the Jet object to assign.
     * @return
     */
    Jet& operator=(const Jet& other);

    // ------------------------------------- Other methods ------------------------------------ //

    /**
     * Calculates a Jet resulting from raising this Jet to a given power.
     *
     * @param power is the power to raise this Jet to.
     * @return the Jet resulting from the power.
     */
    [[nodiscard]] Jet pow(int power) const;

    /**
     * Differentiates the Jet object.
     *
     * @return the Jet of the function's first derivative.
     */
    [[nodiscard]] Jet diff() const;

    /**
     * Prints the values in the derivatives array of this Jet.
     */
    void print();

    /**
     * Returns the vector containing all non-zero derivatives of the function represented by this
     * Jet (up to the order of the highest derivative calculated).
     *
     * @return the vector containing all non-zero derivatives.
     */
    std::vector<double> getDerivatives();

    /**
     * Returns the value of the given derivative of the function represented by this Jet.
     *
     * @param order is the order of the derivative to be returned.
     * @return the value of the given derivative.
     * @throws std::invalid_argument if the given order is higher than the order of the highest
     * derivative calculated or if the given order is negative.
     */
    double getDerivative(int order);

    /**
     * Returns the value of the function represented by this Jet.
     *
     * @return the value of the function.
     */
    double val();
};

// Friend operators
Jet operator/(double factor, const Jet& jet);
Jet operator*(double factor, const Jet& jet);
Jet operator+(double other, const Jet& jet);
Jet operator-(double other, const Jet& jet);

#endif //RECURSIVE_FAD_JET_H