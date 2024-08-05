#ifndef QUAD_H
#define QUAD_H

#include "Hittable.h"

class Quad: public Hittable {
    float3 corner, a, b, topRight, topLeft, N, w;
    float D;
    Material* material;

public:
    ~Quad() override = default;

    Quad(const float3 &corner, const float3 &a, const float3 &b, const float3 &topRight, const float3 &topLeft, Material* mat)
        : corner(corner), a(a), b(b), topRight(topRight), topLeft(topLeft), material(mat) {
        float3 crossAB = cross(a, b);
        N = normalize(crossAB);
        D = dot(corner, N);
        w = crossAB / length2(crossAB);
    }

    bool hit(const Ray &ray, float rayTMin, float rayTMax, HitInfo &hit) const override {
        float denom = dot(N, ray.d);

        // parallel case
        if (abs(denom) < 1e-8) return false;

        float t = (D - dot(N, ray.o)) / denom;
        if (t < rayTMin || t > rayTMax) return false;

        float3 P = ray.o + t * ray.d;
        float3 planarVec = P - corner;
        float alpha = dot(w, cross(a, planarVec));
        float beta = dot(w, cross(planarVec, b));

        if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1) return false;

        hit.t = t;
        hit.P = P;
        hit.T = float2(alpha, beta);
        hit.material = material;
        return true;
    }

    void rasterize(const float4x4 &plm) override {
        if (!shouldRasterize) return;
        // Define the four corners of the quad
        std::vector<float> quadXs = { corner.x, corner.x + a.x, topRight.x, topLeft.x };
        std::vector<float> quadYs = { corner.y, corner.y + a.y, topRight.y, topLeft.y };

        int minX = *std::min_element(quadXs.begin(), quadXs.end());
        int maxX = *std::max_element(quadXs.begin(), quadXs.end());
        int minY = *std::min_element(quadYs.begin(), quadYs.end());
        int maxY = *std::max_element(quadYs.begin(), quadYs.end());

        minX = std::max(minX, 0);
        maxX = std::min(maxX, FrameBuffer.width - 1);
        minY = std::max(minY, 0);
        maxY = std::min(maxY, FrameBuffer.height - 1);

        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                if (insideQuad(quadXs, quadYs, x + 0.5f, y + 0.5f)) {
                    float z = 0.0f; // Calculate appropriate depth value if needed
                    if (z < FrameBuffer.depth(x, y)) {
                        FrameBuffer.depth(x, y) = z;

                        HitInfo hit;
                        hit.T = float2(0, 0); // Calculate texture coordinates if needed
                        hit.material = material;

                        float3 color = shade(hit, {0, 0, 0}, 0);
                        FrameBuffer.pixel(x, y) = color;
                    }
                }
            }
        }
    }

    static bool insideQuad(const std::vector<float>& quadXs, const std::vector<float>& quadYs, float px, float py) {
        int windingNumber = 0;
        for (int i = 0; i < 4; ++i) {
            float currX = quadXs[i];
            float currY = quadYs[i];
            float nxtX = quadXs[(i + 1) % 4];
            float nxtY = quadYs[(i + 1) % 4];

            if (py > std::min(currY, nxtY)
                && py <= std::max(currY, nxtY)
                && px <= std::max(currX, nxtX)) {
                float xIntersect = (py - currY) * (nxtX - currX) / (nxtY - currY) + currX;

                if (currX == nxtX || px <= xIntersect) {
                    windingNumber++;
                }
                }
        }

        return windingNumber % 2 == 1;
    }

    bool didCollide(float3& center, float3 radius, float& sphMoveX, float& sphMoveY) override {
        float maxRadius = linalg::max(linalg::max(radius.x, radius.y), radius.z);
        // Compute the distance from the sphere center to the quad plane
        float distance = dot(N, center) - D;

        // Check if the sphere intersects the plane
        if (abs(distance) > maxRadius) return false;

        // Project the sphere center onto the plane
        float3 projectedCenter = center - distance * N;

        // Convert the quad vertices to local coordinates
        float3 localProj = projectedCenter - corner;
        float alpha = dot(localProj, a) / length2(a);
        float beta = dot(localProj, b) / length2(b);

        // Check if the projected point is inside the quad
        // sphMoveY = -1 * sphMoveY * 0.6;
        if (topLeft.y == topRight.y) {
            center.y = corner.y + radius.y;
            sphMoveY = 0;
        } else {
            center.x = corner.x + radius.x;
            sphMoveX = 0;
        }

        return alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1;
    }

    void step(std::vector<Hittable *>& hittableObjects, std::vector<Fluid*>& fluids, int frameCount) override {

    }

    float3 getRandomPoint() override {
        return float3(0.0f);
    }
};

#endif //QUAD_H
