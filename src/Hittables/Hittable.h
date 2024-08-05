#ifndef HITTABLE_H
#define HITTABLE_H

#include "../helpers.h"
#include "../Raytracer/Ray.h"
#include "HitInfo.h"

class Fluid;

using namespace linalg::aliases;

class Hittable {
public:
    int id = 0;
    bool shouldRasterize = true;
    bool hasGravity = false;
    float moveX = 0;
    float moveY = 0;
    virtual ~Hittable() = default;
    virtual bool hit(const Ray &ray, float rayTMin, float rayTMax, HitInfo &hit) const = 0;
    virtual void rasterize(const float4x4 &plm) = 0;
    virtual float3 getRandomPoint() {
        return float3(0.0f);
    }
    virtual void step(std::vector<Hittable*>& hittableObjects, std::vector<Fluid*>& fluids, int frameCount) {
        if (hasGravity == true) moveY = moveY + globalGravity.y;
    }
    virtual bool didCollide(float3 &center, float3 radius, float &moveX, float& moveY) = 0;
};

#endif //HITTABLE_H
