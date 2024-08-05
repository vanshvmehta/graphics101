#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include "Hittable.h"

class Fluid {
public:
    Fluid(float startX, float startY, int width, int height, float density)
    : startX(startX), startY(startY), width(width), height(height), density(density),
      u(width * height, 0.0f), v(width * height, 0.0f),
      u_prev(width * height, 0.0f), v_prev(width * height, 0.0f),
      p(width * height, 0.0f), div(width * height, 0.0f),
      s(width * height, 0.0f) {
        initializeDensity();
        initializeVelocity();
    }

    void initializeDensity() {
        int centerX = static_cast<int>(startX);
        int centerY = static_cast<int>(startY);
        int boxWidth = 450;  // Define the size of the box
        int boxHeight = 450;

        int fluidWidth = 100;  // Width of initial fluid block
        int fluidHeight = 100; // Height of initial fluid block

        int fluidStartX = centerX - fluidWidth / 2;
        int fluidStartY = centerY - fluidHeight / 2;

        // Add some initial velocity perturbations
        for (int x = centerX - boxWidth / 2; x < centerX + boxWidth / 2; ++x) {
            for (int y = centerY - boxHeight / 2; y < centerY + boxHeight / 2; ++y) {
                int idx = index(x, y);
                s[idx] = 0;
                if (x >= fluidStartX && x < fluidStartX + fluidWidth &&
                    y >= fluidStartY && y < fluidStartY + fluidHeight) {
                    s[idx] = density;
                }
            }
        }
    }


    void initializeVelocity() {
        int centerX = static_cast<int>(startX);
        int centerY = static_cast<int>(startY);
        int boxWidth = 100;  // Define the size of the box
        int boxHeight = 100;

        // Add some initial velocity perturbations
        for (int x = centerX - boxWidth / 2; x < centerX + boxWidth / 2; ++x) {
            for (int y = centerY - boxHeight / 2; y < centerY + boxHeight / 2; ++y) {
                int idx = index(x, y);
                u[idx] = 0.1f * (PCG32::rand() / RAND_MAX - 0.5f);  // Random small initial velocity
                v[idx] = 0.1f * (PCG32::rand() / RAND_MAX - 0.5f);
            }
        }
    }

    void step(std::vector<Hittable*> &hittableObjects, float dt) {
        applyGravity(dt);
        addSource(dt);
        advect(dt);
        handleCollisions(hittableObjects);
        project(dt, hittableObjects);
        setBoundaryConditions(hittableObjects);
    }

    void handleCollisions(const std::vector<Hittable*>& hittableObjects) {
        for (int x = 1; x < width - 1; ++x) {
            for (int y = 1; y < height - 1; ++y) {
                int idx = index(x, y);
                float3 position = float3(x, y, 0);

                for (const auto& object : hittableObjects) {
                    if (object->didCollide(position, float3{1, 1, 0}, u[idx], v[idx])) {
                        // Reverse the velocity component upon collision
                        // Adjust the velocity vector based on the object's surface normal or behavior
                        u[idx] = -u[idx]; // Example: reverse x velocity
                        v[idx] = -v[idx]; // Example: reverse y velocity
                    }
                }
            }
        }
    }


    void applyGravity(float dt) {
        // Apply gravity to the vertical velocity component
        for (int i = 0; i < width * height; ++i) {
            v[i] += dt * globalGravity.y;  // gravity acceleration
        }
    }

    void addSource(float dt) {
        int centerX = static_cast<int>(startX);
        int centerY = static_cast<int>(startY);
        int boxWidth = 450;  // Define the size of the box
        int boxHeight = 450;

        for (int x = centerX - boxWidth / 2; x < centerX + boxWidth / 2; ++x) {
            for (int y = centerY - boxHeight / 2; y < centerY + boxHeight / 2; ++y) {
                int idx = index(x, y);
                u[idx] += dt * u_prev[idx];
                v[idx] += dt * v_prev[idx];
                s[idx] += dt * density;
            }
        }
    }

    void advect(float dt) {
        std::vector<float> u_new = u;
        std::vector<float> v_new = v;
        std::vector<float> s_new = s;

        for (int x = 1; x < width - 1; ++x) {
            for (int y = 1; y < height - 1; ++y) {
                int idx = index(x, y);
                float x0 = x - dt * u[idx];
                float y0 = y - dt * v[idx];

                if (x0 < 0.5f) x0 = 0.5f;
                if (x0 > width - 1.5f) x0 = width - 1.5f;
                int i0 = static_cast<int>(x0);
                int i1 = i0 + 1;

                if (y0 < 0.5f) y0 = 0.5f;
                if (y0 > height - 1.5f) y0 = height - 1.5f;
                int j0 = static_cast<int>(y0);
                int j1 = j0 + 1;

                float s1 = x0 - i0;
                float s0 = 1.0f - s1;
                float t1 = y0 - j0;
                float t0 = 1.0f - t1;

                s_new[idx] = s0 * (t0 * s[index(i0, j0)] + t1 * s[index(i0, j1)]) +
                             s1 * (t0 * s[index(i1, j0)] + t1 * s[index(i1, j1)]);

                u_new[idx] = s0 * (t0 * u[index(i0, j0)] + t1 * u[index(i0, j1)]) +
                             s1 * (t0 * u[index(i1, j0)] + t1 * u[index(i1, j1)]);

                v_new[idx] = s0 * (t0 * v[index(i0, j0)] + t1 * v[index(i0, j1)]) +
                             s1 * (t0 * v[index(i1, j0)] + t1 * v[index(i1, j1)]);
            }
        }

        u = u_new;
        v = v_new;
        s = s_new;
    }

    void project(float dt, const std::vector<Hittable*>& hittableObjects) {
        float h = 1.0f / std::max(width, height);
        for (int x = 1; x < width - 1; ++x) {
            for (int y = 1; y < height - 1; ++y) {
                int idx = index(x, y);
                div[idx] = -0.5f * h * (u[index(x + 1, y)] - u[index(x - 1, y)] +
                                        v[index(x, y + 1)] - v[index(x, y - 1)]);
                p[idx] = 0;
            }
        }
        setBoundaryConditions(hittableObjects);
        for (int k = 0; k < 20; ++k) {
            for (int x = 1; x < width - 1; ++x) {
                for (int y = 1; y < height - 1; ++y) {
                    int idx = index(x, y);
                    p[idx] = (div[idx] + p[index(x - 1, y)] + p[index(x + 1, y)] +
                              p[index(x, y - 1)] + p[index(x, y + 1)]) / 4;
                }
            }
            setBoundaryConditions(hittableObjects);
        }
        for (int x = 1; x < width - 1; ++x) {
            for (int y = 1; y < height - 1; ++y) {
                int idx = index(x, y);
                u[idx] -= 0.5f * (p[index(x + 1, y)] - p[index(x - 1, y)]) / h;
                v[idx] -= 0.5f * (p[index(x, y + 1)] - p[index(x, y - 1)]) / h;
            }
        }
        setBoundaryConditions(hittableObjects);
    }


    void setBoundaryConditions(const std::vector<Hittable*>& hittableObjects) {
        for (int x = 1; x < width - 1; ++x) {
            for (int y = 1; y < height - 1; ++y) {
                int idx = index(x, y);

                // Check for collisions with boundary and objects
                float3 position = float3(x, y, 0);

                for (const auto& object : hittableObjects) {
                    if (object->didCollide(position, float3{1, 1, 0}, u[idx], v[idx])) {
                        // Reverse or modify the velocity to simulate a bounce
                        u[idx] = -u[idx];
                        v[idx] = -v[idx];
                    }
                }
            }
        }
    }

    void rasterize() const {
        for (int x = 75; x < 675; ++x) {
            for (int y = 75; y < 400; ++y) {
                int idx = index(x, y);

                // Set to blue, only render where there's fluid
                if (idx < s.size() && s[idx] > 0.01f) {  // Threshold for showing fluid
                    // cout << "Rasterizing" << endl;
                    float densityValue = s[idx];
                    float normalizedDensity = linalg::clamp(densityValue / density, 0.0f, 1.0f);

                    // Set color based solely on density for consistent coloring
                    float3 color = float3(0, 0, normalizedDensity);

                    FrameBuffer.pixel(x, y) = color;
                }
            }
        }
    }

    const std::vector<float>& getDensityField() const {
        return s;
    }

    // Check if the sphere is colliding with the fluid
    bool checkCollision(const float3& sphereCenter, const float3& sphereRadius) const {
        float fluidLevel = getFluidLevel(sphereCenter.x);

        // Check if the sphere's center is below the fluid level
        if (sphereCenter.y - sphereRadius.y < fluidLevel) {
            return true; // Collision detected
        }

        return false; // No collision
    }


    // Get the fluid level at a specific x position along the x-axis
    float getFluidLevel(float x) const {
        int maxY = 0;
        for (int x = 1; x < width - 1; ++x) {
            for (int y = 1; y < height - 1; ++y) {
                int idx = index(x, y);

                if (idx < s.size() && s[idx] > 0.05f) {
                    maxY = std::max(y, maxY);
                }
            }
        }

        cout << maxY << endl;
        return (float)maxY;
    }



private:
    float startX, startY;
    int width, height;
    float density;
    std::vector<float> u, v; // velocity field
    std::vector<float> u_prev, v_prev;
    std::vector<float> p, div; // pressure and divergence
    std::vector<float> s; // smoke/density field

    int index(int x, int y) const {
        return x + y * width;
    }
};


#endif //FLUID_H
