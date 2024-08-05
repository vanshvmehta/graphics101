#ifndef RAY_H
#define RAY_H

#include "../helpers.h"

class Ray {
public:
    float3 o, d;
    Ray() : o(), d(float3(0.0f, 0.0f, 1.0f)) {}
    Ray(const float3& o, const float3& d) : o(o), d(d) {}
};

#endif //RAY_H
