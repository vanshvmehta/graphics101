#include "cs488.h"
#include "Hittables/Quad.h"
#include "Hittables/Sphere.h"
#include "Scenes/Animations.h"
CS488Window CS488;

 
// draw something in each frame
static void draw() {
    for (int j = 0; j < globalHeight; j++) {
        for (int i = 0; i < globalWidth; i++) {
            // FrameBuffer.pixel(i, j) = float3(PCG32::rand()); // noise
            FrameBuffer.pixel(i, j) = float3(0.5f * (cos((i + globalFrameCount) * 0.1f) + 1.0f)); // moving cosine
        }
    }
}

// setting up lighting
static PointLightSource light;
static void setupLightSource() {
    //light.position = float3(0.5f, 4.0f, 1.0f); // use this for sponza.obj
    light.position = float3(3.0f, 3.0f, 3.0f);
    light.wattage = float3(1000.0f, 1000.0f, 1000.0f);
    globalScene.addLight(&light);
}



// ======== you probably don't need to modify below in A1 to A3 ========
// loading .obj file from the command line arguments
static TriangleMesh mesh;
static void setupScene(int argc, const char* argv[]) {
    if (argc > 1) {
        bool objLoadSucceed = mesh.load(argv[1]);
        if (!objLoadSucceed) {
            printf("Invalid .obj file.\n");
            printf("Making a single triangle instead.\n");
            mesh.createSingleTriangle();
        }
    } else {
        printf("Specify .obj file in the command line arguments. Example: CS488.exe cornellbox.obj\n");
        printf("Making a single triangle instead.\n");
        mesh.createSingleTriangle();
    }
    globalScene.addObject(&mesh);

    float3 sphereCenter = float3(375, 375, 50); // Place the sphere 5 units away along the z-axis

    // Define the radius of the sphere
    float sphereRadius = 75.0f;

    auto red = new Material();
    red->Ka = float3(.65, .05, .05);

    auto white = new Material();
    white->Ka = float3(.73, .73, .73);

    auto green = new Material();
    green->Ka = float3(.12, .45, .15);

    // Create the sphere object
    auto sphere = new Sphere(sphereCenter, sphereRadius, white);
    sphere->moveY = -5;
    globalScene.addHittableObject(sphere);

    float sideLength = 500;

    // Create the walls of the cube
    float3 bottomLeft = float3(150, 75, 0);

    // Define offsets for trapezoid adjustments
    float xOffset = 1;
    float yOffset = 1;

    // Right wall
    auto floor = new Quad(bottomLeft,
                          float3(sideLength, 0, 0),
                          float3(xOffset, yOffset, sideLength),
                          bottomLeft + float3(sideLength - xOffset, yOffset, sideLength),
                          bottomLeft + float3(xOffset, yOffset, sideLength),
                          green);
    floor->shouldRasterize = false;
    globalScene.addHittableObject(floor);
}

static void A1(int argc, const char* argv[]) {
    setupScene(argc, argv);
    setupLightSource();
    globalRenderType = RENDER_RASTERIZE;
}


int main(int argc, const char* argv[]) {
    auto animation = new Animations();
    animation->setupScene(argc, argv);

    CS488.start(animation);
}
