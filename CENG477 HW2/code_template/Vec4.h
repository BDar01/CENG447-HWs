#ifndef __VEC4_H__
#define __VEC4_H__
#include "Vec3.h"
#define NO_COLOR -1
#include <ostream>
#include <iostream>
class Vec4
{
public:
    double x, y, z, t;
    int colorId;

    Vec4();
    Vec4(double x, double y, double z, double t);
    Vec4(double x, double y, double z, double t, int colorId);
    Vec4(Vec3& v, int cId);
    Vec4(const Vec4 &other);

    double getNthComponent(int n);

    friend std::ostream &operator<<(std::ostream &os, const Vec4 &v);

    Vec4 operator+(const Vec4& other) const {
        return Vec4(x + other.x, y + other.y, z + other.z, t + other.t);
    }

    // Overload the - operator
    Vec4 operator-(const Vec4& other) const {
        return Vec4(x - other.x, y - other.y, z - other.z, t - other.t);
    }
    Vec4& operator+=(const Vec4& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        t += other.t;
        return *this;
    }
    Vec4& operator-=(const Vec4& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        t -= other.t;
        return *this;
    }
    Vec4 operator*(const double& other) const {
        return Vec4(x * other, y * other, z * other, t * other);
    }
    Vec4 operator/(const double& other) const {
        return Vec4(x / other, y / other, z / other, t / other);
    }

};

#endif