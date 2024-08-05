#ifndef ANIMATIONS_H
#define ANIMATIONS_H

#include "Scene.h"
#include "../helpers.h"
#include "../Hittables/Cuboid.h"
#include "../Hittables/Fluid.h"
#include "../Hittables/Quad.h"
#include "../Hittables/Sphere.h"
#include "../Hittables/TriangleMesh.h"


class Animations {
private:
    Material* red = new Material();
    Material* white = new Material();
    Material* green = new Material();
    TriangleMesh mesh;
    PointLightSource light;
    std::vector<Sphere*> spheres;
    std::vector<Quad*> quads;
    std::vector<Cuboid*> cuboids;
    std::vector<Fluid*> fluids;

    int startSphereFallingFrame = 0;
    float fallingFactor = 0;

    int framesToStretch = 13;
    int framesToFall = 10;

    int startSqueezingFrame = startSphereFallingFrame + framesToFall;
    float squeezingFactor = 1.5;

    int startStretchingFrame = startSqueezingFrame + framesToStretch;
    float stretchingFactor = 1.5;

    int startNormalizingFrame = startStretchingFrame + framesToStretch * 2;

    int framesToMoveCuboid = 25;
    float cuboidMovingFactor = 10;
    int startCuboidMovingFrame = startNormalizingFrame + framesToStretch;
    int startCuboidRetractingFrame = startCuboidMovingFrame + framesToMoveCuboid;

    int startNormlizingSphereAgain = startCuboidRetractingFrame + framesToMoveCuboid;

    int framesToFlowFluid = 20;
    int startFluidFrame = startNormlizingSphereAgain + framesToMoveCuboid;
    int endFluidFlowFrame = startFluidFrame + framesToFlowFluid;

    int endFrame = endFluidFlowFrame + 100;
public:
    Animations() {
        red->Ka = float3(.65, .05, .05);
        white->Ka = float3(.73, .73, .73);
        green->Ka = float3(.12, .45, .15);
    };
    ~Animations() = default;

    void setupScene(int argc, const char* argv[]) {
        // Add Lighting
        light.position = float3(3.0f, 3.0f, 3.0f);
        light.wattage = float3(1000.0f, 1000.0f, 1000.0f);
        globalScene.addLight(&light);

        // Add the Cornellbox Object
        if (argc > 1) mesh.load(argv[1]);
        // globalScene.addObject(&mesh);

        // Create the main sphere
        auto sphere = new Sphere(float3(375, 600, 50), 75.0f, white);
        sphere->id = 1;
        sphere->hasGravity = true;
        // globalScene.addHittableObject(sphere);
        spheres.push_back(sphere);

        // Create the Floor Quad for collisions (not visible)
        float sideLength = 600;
        float3 bottomLeft = float3(75, 75, 0);
        float xOffset = 20;
        float yOffset = 20;
        auto floor = new Quad(bottomLeft,
                              float3(sideLength, 0, 0),
                              float3(xOffset, yOffset, sideLength),
                              bottomLeft + float3(sideLength - xOffset, yOffset, sideLength),
                              bottomLeft + float3(xOffset, yOffset, sideLength),
                              green);
        floor->shouldRasterize = false;
        floor->id = 2;
        globalScene.addHittableObject(floor);

        // Left wall
        auto leftWall = new Quad(bottomLeft,
                                 float3(0, sideLength, 0),
                                 float3(0, 0, sideLength),
                                 bottomLeft + float3(xOffset, sideLength, 0),
                                 bottomLeft + float3(xOffset, 0, sideLength),
                                 red);
        leftWall->shouldRasterize = false;
        globalScene.addHittableObject(leftWall);

        // Right wall
        float3 rightBottomLeft = bottomLeft + float3(sideLength, 0, 0);
        auto rightWall = new Quad(rightBottomLeft,
                                  float3(0, sideLength, 0),
                                  float3(0, 0, sideLength),
                                  rightBottomLeft + float3(-xOffset, sideLength, 0),
                                  rightBottomLeft + float3(-xOffset, 0, sideLength),
                                  red);
        rightWall->shouldRasterize = false;
        globalScene.addHittableObject(rightWall);

        // Create the cuboid
        auto cuboid = new Cuboid(float3(525, 75, 0), float3(600, 275, 0), green);
        cuboid->id = 3;
        cuboid->shouldRasterize = false;
        globalScene.addHittableObject(cuboid);
        cuboids.push_back(cuboid);

        auto fluid = new Fluid(375, 375, 300, 300, 1.0f);
        fluids.push_back(fluid);
        globalScene.addFluid(fluid);
    }

    void step(int frameCount) const {
        globalScene.step(globalFrameCount);

        if (frameCount >= 40) globalScene.addHittableObject(spheres[0]);

        // if (frameCount == startSphereFallingFrame) {
        //     spheres[0]->moveY = fallingFactor;
        // } else if (frameCount >= startSqueezingFrame && frameCount < startStretchingFrame) {
        //     spheres[0]->squeeze(0, squeezingFactor, 0);
        // } else if (frameCount >= startStretchingFrame && frameCount < startNormalizingFrame) {
        //     spheres[0]->squeeze(0, -stretchingFactor, 0);
        // } else if (frameCount >= startNormalizingFrame && frameCount < startCuboidMovingFrame) {
        //     spheres[0]->squeeze(0, stretchingFactor, 0);
        // } else if (frameCount >= startCuboidMovingFrame && frameCount < startCuboidRetractingFrame) {
        //     cuboids[0]->shouldRasterize = true;
        //     cuboids[0]->moveX = -cuboidMovingFactor;
        // } else if (frameCount >= startCuboidRetractingFrame && frameCount < startNormlizingSphereAgain) {
        //     spheres[0]->squeeze(squeezingFactor * (float)0.75, 0, 0);
        //     cuboids[0]->moveX = cuboidMovingFactor;
        // } else if (frameCount >= startNormlizingSphereAgain && frameCount < startFluidFrame) {
        //     cout << "Changing Squeeze" << endl;
        //     cuboids[0]->moveX = 0;
        //     cuboids[0]->shouldRasterize = false;
        //     spheres[0]->squeeze(squeezingFactor * (float)-0.75, 0, 0);
        // }
    }
};

#endif //ANIMATIONS_H
