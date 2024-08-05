cmake_minimum_required(VERSION 3.16)
project(CS488)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(external/glew)
add_subdirectory(external/glfw)

file(GLOB SRC_FILES
        src/*.h
        src/*.cpp)

add_executable(CS488
        ${SRC_FILES}
        src/helpers.h
        src/helpers.h
        src/Image.h
        src/Hittables/TriangleMesh.h
        src/RayTracer/Ray.h
        src/Hittables/HitInfo.h
        src/Scenes/Scene.h
        src/Hittables/Hittable.h
        src/Hittables/Sphere.h
        src/Hittables/Quad.h
        src/Scenes/Animations.h
        src/Hittables/Cuboid.h
        src/Hittables/Fluid.h)

target_include_directories(CS488 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/glfw/include)
target_include_directories(CS488 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/external/glew/include)
target_link_libraries(CS488 glfw libglew_static)
