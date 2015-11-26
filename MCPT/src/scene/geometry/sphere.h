#pragma once

#include <scene/geometry/geometry.h>

//A sphere is assumed to have a radius of 0.5 and a center of <0,0,0>.
//These attributes can be altered by applying a transformation matrix to the sphere.
class Sphere : public Geometry
{
public:
    Intersection GetIntersection(Ray r);
    virtual glm::vec2 GetUVCoordinates(const glm::vec3 &point);
    virtual glm::vec3 ComputeNormal(const glm::vec3 &P);
    virtual Intersection GetRandISX(float rand1, float rand2, const glm::vec3 &isx_normal);
    virtual float RayPDF(const Intersection &isx, const Ray &ray);
    void create();

    virtual void ComputeArea();
};
