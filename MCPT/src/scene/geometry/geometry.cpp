#include <scene/geometry/geometry.h>

float Geometry::RayPDF(const Intersection &isx, const Ray &ray)
{
    //TODO
    //The isx passed in was tested ONLY against us (no other scene objects), so we test if NULL
    //rather than if != this.
    if(isx.object_hit == NULL)
    {
        return 0;
    }
    //Add more here
    //find the intersection of the ray and myself(light)
    //intersecting again assures that i get the closer point
    //Intersection with_me = GetIntersection(ray);
    float r (glm::distance(ray.origin, isx.point)); //distance from ray to intersection
    float cos_theta(glm::dot(isx.normal, -ray.direction)); //angle between ray and surface normal
    float pdf_SA(r*r/cos_theta/area);

    return pdf_SA;
}

Intersection Geometry::GetRandISX(float rand1, float rand2, const glm::vec3 &isx_normal)
{
    //this is fudge just to get things to compile
    //override this function in individual geometry plz
    Intersection mike;
    return mike;
}

Ray Geometry::GetRandRay(float rand1, float rand2, const Intersection &isx)
{
    //ray goes away from the intersection
    glm::vec3 dir;

    dir.z = rand1;
    dir.x = glm::cos(TWO_PI*rand2)*glm::sqrt(1 - rand1*rand1);
    dir.y = glm::sin(TWO_PI*rand2)*glm::sqrt(1 - rand1*rand1);

    //transform ray into world space
    glm::vec3 dir_world = glm::vec3(transform.T()*glm::vec4(dir, 0));

    return Ray(isx.point, dir_world);
}
