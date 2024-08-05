#ifndef SCENE_H
#define SCENE_H

#include "../helpers.h"
#include "../Hittables/Fluid.h"
#include "../Hittables/Hittable.h"
#include "../Hittables/TriangleMesh.h"

class PointLightSource {
public:
	float3 position, wattage;
};

// scene definition
class Scene {
public:
	std::vector<Hittable*> hittableObjects;
	std::vector<TriangleMesh*> objects;
	std::vector<Fluid*> fluids;
	std::vector<PointLightSource*> pointLightSources;
	std::vector<BVH> bvhs;

	void addHittableObject(Hittable* obj) {
		hittableObjects.push_back(obj);
	}
	void addObject(TriangleMesh* pObj) {
		objects.push_back(pObj);
	}
	void addFluid(Fluid* obj) {
		fluids.push_back(obj);
	}
	void addLight(PointLightSource* pObj) {
		pointLightSources.push_back(pObj);
	}

	void step(int frameCount) {
		// const float4x4 pm = perspectiveMatrix(globalFOV, globalAspectRatio, globalDepthMin, globalDepthMax);
		// const float4x4 lm = lookatMatrix(globalEye, globalLookat, globalUp);
		// const float4x4 plm = mul(pm, lm);

		for (auto hittable: hittableObjects) {
			hittable->step(hittableObjects, fluids, frameCount);
		}

		for (auto fluid: fluids) {
			fluid->step(hittableObjects, 0.1f);
		}
	}

	void preCalc() {
		bvhs.resize(objects.size());
		for (int i = 0; i < objects.size(); i++) {
			objects[i]->preCalc();
			bvhs[i].build(objects[i]);
		}
	}

	// ray-scene intersection
	bool intersect(HitInfo& minHit, const Ray& ray, float tMin = 0.0f, float tMax = FLT_MAX) const {
		bool hit = false;
		HitInfo tempMinHit;
		minHit.t = FLT_MAX;

		for (int i = 0, i_n = (int)objects.size(); i < i_n; i++) {
			//if (objects[i]->bruteforceIntersect(tempMinHit, ray, tMin, tMax)) { // for debugging
			if (bvhs[i].intersect(tempMinHit, ray, tMin, tMax)) {
				if (tempMinHit.t < minHit.t) {
					hit = true;
					minHit = tempMinHit;
				}
			}
		}

		for (auto hittable: hittableObjects) {
			if (hittable->hit(ray, tMin, tMax, tempMinHit)) {
				if (tempMinHit.t < minHit.t) {
					hit = true;
					minHit = tempMinHit;
				}
			}
		}

		return hit;
	}

	// camera -> screen matrix (given to you for A1)
	float4x4 perspectiveMatrix(float fovy, float aspect, float zNear, float zFar) const {
		float4x4 m;
		const float f = 1.0f / (tan(fovy * DegToRad / 2.0f));
		m[0] = { f / aspect, 0.0f, 0.0f, 0.0f };
		m[1] = { 0.0f, f, 0.0f, 0.0f };
		m[2] = { 0.0f, 0.0f, (zFar + zNear) / (zNear - zFar), -1.0f };
		m[3] = { 0.0f, 0.0f, (2.0f * zFar * zNear) / (zNear - zFar), 0.0f };

		return m;
	}

	// model -> camera matrix (given to you for A1)
	static float4x4 lookatMatrix(const float3& _eye, const float3& _center, const float3& _up) {
		// transformation to the camera coordinate
		float4x4 m;
		const float3 f = normalize(_center - _eye);
		const float3 upp = normalize(_up);
		const float3 s = normalize(cross(f, upp));
		const float3 u = cross(s, f);

		m[0] = { s.x, s.y, s.z, 0.0f };
		m[1] = { u.x, u.y, u.z, 0.0f };
		m[2] = { -f.x, -f.y, -f.z, 0.0f };
		m[3] = { 0.0f, 0.0f, 0.0f, 1.0f };
		m = transpose(m);

		// translation according to the camera location
		const float4x4 t = float4x4{ {1.0f, 0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f, 0.0f}, { -_eye.x, -_eye.y, -_eye.z, 1.0f} };

		m = mul(m, t);
		return m;
	}

	// rasterizer
	void Rasterize() const {
		// ====== implement it in A1 ======
		// fill in plm by a proper matrix
		const float4x4 pm = perspectiveMatrix(globalFOV, globalAspectRatio, globalDepthMin, globalDepthMax);
		const float4x4 lm = lookatMatrix(globalEye, globalLookat, globalUp);
		const float4x4 plm = mul(pm, lm);

		FrameBuffer.clear();
		for (auto object : objects) {
			for (int k = 0, k_n = (int)object->triangles.size(); k < k_n; k++) {
				object->rasterizeTriangle(object->triangles[k], plm);
			}
		}

		for (auto hittable: hittableObjects) {
			hittable->rasterize(plm);
		}

		for (auto fluid: fluids) {
			fluid->rasterize();
		}
	}

	// eye ray generation (given to you for A2)
	static Ray eyeRay(int x, int y) {
		// compute the camera coordinate system
		const float3 wDir = normalize(float3(-globalViewDir));
		const float3 uDir = normalize(cross(globalUp, wDir));
		const float3 vDir = cross(wDir, uDir);

		// compute the pixel location in the world coordinate system using the camera coordinate system
		// trace a ray through the center of each pixel
		const float imPlaneUPos = (x + 0.5f) / float(globalWidth) - 0.5f;
		const float imPlaneVPos = (y + 0.5f) / float(globalHeight) - 0.5f;

		const float3 pixelPos = globalEye + float(globalAspectRatio * globalFilmSize * imPlaneUPos) * uDir + float(globalFilmSize * imPlaneVPos) * vDir - globalDistanceToFilm * wDir;

		return Ray(globalEye, normalize(pixelPos - globalEye));
	}

	// ray tracing (you probably don't need to change it in A2)
	void Raytrace() const {
		FrameBuffer.clear();

		// loop over all pixels in the image
		for (int j = 0; j < globalHeight; ++j) {
			for (int i = 0; i < globalWidth; ++i) {
				const Ray ray = eyeRay(i, j);
				HitInfo hitInfo;
				if (intersect(hitInfo, ray)) {
					FrameBuffer.pixel(i, j) = shade(hitInfo, -ray.d);
				} else {
					FrameBuffer.pixel(i, j) = float3(0.0f);
				}
			}

			// show intermediate process
			if (globalShowRaytraceProgress) {
				constexpr int scanlineNum = 64;
				if ((j % scanlineNum) == (scanlineNum - 1)) {
					glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, globalWidth, globalHeight, GL_RGB, GL_FLOAT, &FrameBuffer.pixels[0]);
					glRecti(1, 1, -1, -1);
					glfwSwapBuffers(globalGLFWindow);
					printf("Rendering Progress: %.3f%%\r", j / float(globalHeight - 1) * 100.0f);
					fflush(stdout);
				}
			}
		}
	}

};
static Scene globalScene;

#endif //SCENE_H
