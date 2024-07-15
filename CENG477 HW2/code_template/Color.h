#ifndef __COLOR_H__
#define __COLOR_H__
#include <ostream>
#include <iostream>
class Color
{
public:
    double r, g, b;

    Color();
    Color(double r, double g, double b);
    Color(const Color &other);
    friend std::ostream &operator<<(std::ostream &os, const Color &c);
    Color operator+(const Color& other) const {
        return Color(r + other.r, g + other.g, b + other.b);
    }

    // Overload the - operator
    Color operator-(const Color& other) const {
        return Color(r - other.r, g - other.g, b - other.b);
    }
    Color& operator+=(const Color& other) {
        r += other.r;
        g += other.g;
        b += other.b;
        return *this;
    }
    Color& operator-=(const Color& other) {
        r -= other.r;
        g -= other.g;
        b -= other.b;
        return *this;
    }
    Color operator*(const Color& other) const {
        return Color(r * other.r, g * other.g, b * other.b);
    }
    Color operator*(const double& other) const {
        return Color(r * other, g * other, b * other);
    }
    Color operator/(const Color& other) const {
        return Color(r / other.r, g / other.g, b / other.b);
    }
};

#endif