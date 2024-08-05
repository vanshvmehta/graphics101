#ifndef SPHERE_H
#define SPHERE_H

#include "Hittable.h"
#include "Fluid.h"

class Sphere : public Hittable {
    float3 center;
    float radiusX, radiusY, radiusZ; // Different radii for different axes
    Material* material;

public:
    ~Sphere() override = default;

    Sphere(const float3 &center, float radius, Material* mat)
        : center(center), radiusX(radius), radiusY(radius), radiusZ(radius), material(mat) {}

    bool hit(const Ray &ray, float rayTMin, float rayTMax, HitInfo &hit) const override {
        float3 oc = ray.o - center;
        float a = dot(ray.d / float3(radiusX, radiusY, radiusZ), ray.d / float3(radiusX, radiusY, radiusZ));
        float b = 2.0f * dot(oc / float3(radiusX, radiusY, radiusZ), ray.d / float3(radiusX, radiusY, radiusZ));
        float c = dot(oc / float3(radiusX, radiusY, radiusZ), oc / float3(radiusX, radiusY, radiusZ)) - 1;
        float discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            return false;
        } else {
            float sqrtd = sqrt(discriminant);
            float t = (-b - sqrtd) / (2.0f * a);
            if (t < rayTMin || t > rayTMax) {
                t = (-b + sqrtd) / (2.0f * a);
                if (t < rayTMin || t > rayTMax) return false;
            }

            hit.t = t;
            hit.P = ray.o + t * ray.d;
            hit.material = material;
            hit.T = float2(0.0f, 0.0f); // Placeholder for texture coordinates
            return true;
        }
    }

    void rasterize(const float4x4 &plm) override {
        float3 minPoint = center - float3(radiusX, radiusY, radiusZ);
        float3 maxPoint = center + float3(radiusX, radiusY, radiusZ);

        int minX = std::max(0, (int)(minPoint.x));
        int maxX = std::min(FrameBuffer.width - 1, (int)(maxPoint.x));
        int minY = std::max(0, (int)(minPoint.y));
        int maxY = std::min(FrameBuffer.height - 1, (int)(maxPoint.y));

        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                if (insideEllipse((float)x + 0.5f, (float)y + 0.5f)) {
                    float dx = (float)x + 0.5f - center.x;
                    float dy = (float)y + 0.5f - center.y;
                    float dz2 = radiusZ * radiusZ - dx * dx * (radiusZ * radiusZ / (radiusX * radiusX)) - dy * dy * (radiusZ * radiusZ / (radiusY * radiusY));
                    if (dz2 > 0) {
                        float dz = sqrt(dz2);
                        float depth = center.z - dz; // Calculate the closer z coordinate

                        if (depth < FrameBuffer.depth(x, y)) {
                            FrameBuffer.depth(x, y) = depth;

                            HitInfo hit;
                            hit.T = {0.0f, 0.0f}; // Texture coordinates are not used
                            hit.material = material;

                            float3 color = shade(hit, {0, 0, 0}, 0);
                            FrameBuffer.pixel(x, y) = color;
                        }
                    }
                }
            }
        }
    }

    float3 getRandomPoint() override {
        return float3(0.0f);
    }

    bool didCollide(float3& center, float3 radius, float& sphMoveX, float& sphMoveY) override {
        return false;
        // float3 diff = (this->center - center) / float3(radiusX, radiusY, radiusZ);
        // float distance = length(diff);
        // return distance <= 1.0f + radius / std::max({radiusX, radiusY, radiusZ});
    }

    void step(std::vector<Hittable*>& hittableObjects, std::vector<Fluid*>& fluids, int frameCount) override {
        Hittable::step(hittableObjects, fluids, frameCount);

        center.y += moveY;
        center.x += moveX;

        // Check collision with other hittable objects
        for (auto& object : hittableObjects) {
            if (object->id == id) continue;
            if (object->didCollide(center, float3{radiusX, radiusY, radiusZ}, moveX, moveY)) {
                // cout << frameCount << endl;
            }
        }

        // Check collision with fluids and adjust position
        for (auto& fluid : fluids) {
            if (fluid->checkCollision(center, float3{radiusX, radiusY, radiusZ})) {
                cout << "Fluid Collision" << endl;
                // Get fluid level at the sphere's x, z position
                float fluidLevel = fluid->getFluidLevel(center.x);
                cout << fluidLevel << endl;
                float sphereBase = center.y - radiusY;

                if (sphereBase < fluidLevel) {
                    // Adjust the sphere's position to float on the fluid
                    center.y = fluidLevel + radiusY;
                }
            }
        }
    }


    void squeeze(float amountX, float amountY, float amountZ) {
        cout << amountX << endl;
        radiusX = std::max(0.0f, radiusX - amountX);
        radiusY = std::max(0.0f, radiusY - amountY);
        radiusZ = std::max(0.0f, radiusZ - amountZ);

        center.y -= amountY;
        center.x -= amountX;
    }

private:
    bool insideEllipse(float px, float py) const {
        float dx = (px - center.x) / radiusX;
        float dy = (py - center.y) / radiusY;
        return dx * dx + dy * dy <= 1.0f; // This checks within the ellipse projection of the sphere on the xy-plane
    }
};

#endif //SPHERE_H
