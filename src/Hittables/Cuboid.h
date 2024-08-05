#ifndef CUBOID_H
#define CUBOID_H

#include "Hittable.h"

class Cuboid : public Hittable {
    float3 minCorner;
    float3 maxCorner;
    Material* material;

public:
    ~Cuboid() override = default;

    Cuboid(const float3 &minCorner, const float3 &maxCorner, Material* mat)
        : minCorner(minCorner), maxCorner(maxCorner), material(mat) {}

    bool hit(const Ray &ray, float rayTMin, float rayTMax, HitInfo &hit) const override {
        float tMin = (minCorner.x - ray.o.x) / ray.d.x;
        float tMax = (maxCorner.x - ray.o.x) / ray.d.x;
        if (tMin > tMax) std::swap(tMin, tMax);

        float tyMin = (minCorner.y - ray.o.y) / ray.d.y;
        float tyMax = (maxCorner.y - ray.o.y) / ray.d.y;
        if (tyMin > tyMax) std::swap(tyMin, tyMax);

        if ((tMin > tyMax) || (tyMin > tMax))
            return false;

        if (tyMin > tMin)
            tMin = tyMin;

        if (tyMax < tMax)
            tMax = tyMax;

        float tzMin = (minCorner.z - ray.o.z) / ray.d.z;
        float tzMax = (maxCorner.z - ray.o.z) / ray.d.z;
        if (tzMin > tzMax) std::swap(tzMin, tzMax);

        if ((tMin > tzMax) || (tzMin > tMax))
            return false;

        if (tzMin > tMin)
            tMin = tzMin;

        if (tzMax < tMax)
            tMax = tzMax;

        if ((tMin < rayTMax) && (tMax > rayTMin)) {
            hit.t = tMin;
            hit.P = ray.o + tMin * ray.d;
            hit.material = material;
            hit.T = float2(0.0f, 0.0f); // Placeholder for texture coordinates
            return true;
        }
        return false;
    }

    void rasterize(const float4x4 &plm) override {
        if (!shouldRasterize) return;

        float3 vertices[8] = {
            minCorner,
            {minCorner.x, minCorner.y, maxCorner.z},
            {minCorner.x, maxCorner.y, minCorner.z},
            {minCorner.x, maxCorner.y, maxCorner.z},
            {maxCorner.x, minCorner.y, minCorner.z},
            {maxCorner.x, minCorner.y, maxCorner.z},
            {maxCorner.x, maxCorner.y, minCorner.z},
            {maxCorner.x, maxCorner.y, maxCorner.z}
        };

        int faces[6][4] = {
            {0, 1, 5, 4}, // Front face
            {2, 3, 7, 6}, // Back face
            {0, 2, 6, 4}, // Left face
            {1, 3, 7, 5}, // Right face
            {0, 1, 3, 2}, // Bottom face
            {4, 5, 7, 6}  // Top face
        };

        for (const auto& face : faces) {
            std::vector<float> xCoords = {vertices[face[0]].x, vertices[face[1]].x, vertices[face[2]].x, vertices[face[3]].x};
            std::vector<float> yCoords = {vertices[face[0]].y, vertices[face[1]].y, vertices[face[2]].y, vertices[face[3]].y};
            std::vector<float> zCoords = {vertices[face[0]].z, vertices[face[1]].z, vertices[face[2]].z, vertices[face[3]].z};

            int minX = std::max(0, static_cast<int>(*std::min_element(xCoords.begin(), xCoords.end())));
            int maxX = std::min(FrameBuffer.width - 1, static_cast<int>(*std::max_element(xCoords.begin(), xCoords.end())));
            int minY = std::max(0, static_cast<int>(*std::min_element(yCoords.begin(), yCoords.end())));
            int maxY = std::min(FrameBuffer.height - 1, static_cast<int>(*std::max_element(yCoords.begin(), yCoords.end())));

            for (int y = minY; y <= maxY; ++y) {
                for (int x = minX; x <= maxX; ++x) {
                    if (insideQuad(xCoords, yCoords, static_cast<float>(x) + 0.5f, static_cast<float>(y) + 0.5f)) {
                        float depth = interpolateDepth(xCoords, yCoords, zCoords, static_cast<float>(x) + 0.5f, static_cast<float>(y) + 0.5f);
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
        float maxRadius = linalg::max(linalg::max(radius.x, radius.y), radius.z);
        float3 closestPoint = clamp(center, minCorner, maxCorner);
        if (length(closestPoint - center) <= maxRadius) {
            if (sphMoveX == 0) center.x += 4*moveX;
            sphMoveX = 2*moveX;
            return true;
        }

        return false;
    }

    void step(std::vector<Hittable*>& hittableObjects, std::vector<Fluid*>& fluids, int frameCount) override {
        minCorner.x += moveX;
        maxCorner.x += moveX;
    }

private:
    static bool insideQuad(const std::vector<float>& xCoords, const std::vector<float>& yCoords, float px, float py) {
        // Implement your point-in-polygon algorithm here (e.g., winding number algorithm)
        // This is a placeholder for the actual implementation
        int windingNumber = 0;
        for (int i = 0; i < 4; ++i) {
            float currX = xCoords[i];
            float currY = yCoords[i];
            float nxtX = xCoords[(i + 1) % 4];
            float nxtY = yCoords[(i + 1) % 4];

            if (py > std::min(currY, nxtY)
                && py <= std::max(currY, nxtY)
                && px <= std::max(currX, nxtX)) {
                float xIntersect = (py - currY) * (nxtX - currX) / (nxtY - currY) + currX;
                if (currX == nxtX || px <= xIntersect)
                    windingNumber += (currY <= py) ? 1 : -1;
            }
        }
        return windingNumber != 0;
    }

    static float interpolateDepth(const std::vector<float>& xCoords, const std::vector<float>& yCoords, const std::vector<float>& zCoords, float px, float py) {
        // Implement your depth interpolation here
        // This is a placeholder for the actual implementation
        float totalArea = 0.0f;
        float depth = 0.0f;
        for (int i = 0; i < 4; ++i) {
            float x1 = xCoords[i];
            float y1 = yCoords[i];
            float z1 = zCoords[i];
            float x2 = xCoords[(i + 1) % 4];
            float y2 = yCoords[(i + 1) % 4];
            float z2 = zCoords[(i + 1) % 4];
            float area = 0.5f * ((x1 - px) * (y2 - py) - (x2 - px) * (y1 - py));
            totalArea += area;
            depth += area * (z1 + z2) / 2;
        }
        return depth / totalArea;
    }
};

#endif // CUBOID_H
