#include "photonmap.h"




PhotonNode::PhotonNode()
{
   this->leftChild = NULL;
   this->rightChild = NULL;
//   photo_area = {0,0,0,0,0,0,0,0,0,0};
}

PhotonMap::PhotonMap():num_bounces(3),mersenne_generator(rand()), unif_distribution(0.0f,1.0f),max_photons(1000)
{
   scene = NULL;
   intersection_engine = NULL;
}

PhotonNode* PhotonMap::createDirectPhotonMap(PhotonNode* root, QList <Geometry*> &scene_geom, QList<Photon*> &photons )
{
    /*
     *  random Ray from light to somewhere on the scene
     * if it intersects between light and diffuse put in Direct Photon Map
     * Indirect Photon map ~ global illumination for all surfaces (no idea what that means)
     * Radiance Map - photons that represent a computed average radience of nearby photon hits in a neighbor hood
     * Caustic Map - Map of Photons that hit at least one specular like surface
     * These are K-d Tree
     * Can store photon multiple times
     * Photons with a surface are computed with BRDF/ BTDF/ Absorbed (decided with russian roulette)
     */
    for (Geometry* Lite: lights)

    {

     }
    return NULL;
}
//Ray PhotonMap::createRandRay(float u1, float u2, float l1, float l2)
//{

//}
// flip this around so geometry makes ray to light
Intersection PhotonMap::RandomSampleLight()
{
    int light_choice = 0;
    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);

    light_choice = rand()%scene->lights.count();
    Geometry* light_obj = scene->lights.at(light_choice);
    light_obj->ComputeNormal()
//    Intersection light_isx = light_obj->GetRandISX(rand1, rand2, light_obj->)

    //return isx;
}

PhotonNode PhotonMap::placePhotons(const Intersection &light_isx)
{

    Photon P;
    glm::vec3 wiW;
    int geom_choice = 0;
    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    geom_choice = rand()%scene->objects.count();
    Geometry* scene_obj = scene->objects.at(geom_choice);
    Intersection obj_isx = scene_obj->GetRandISX(rand1, rand2,light_isx.normal);

    wiW = obj_isx.point - light_isx.point;
    obj_isx.t = glm::length(wiW);
    wiW = glm::normalize(wiW);

    Ray  photon_ray = Ray(light_isx.point, wiW); // direction from light source to point on scene
    QList<Intersection>obstruction_test = intersection_engine->GetAllIntersections(photon_ray);

    for (Intersection scene_isx: obstruction_test)
    {
        if(scene_isx.object_hit == *scene_obj)
        {
            P = createPhoton(scene_isx.point,light_isx.object_hit->material->intensity/1000, photon_ray.direction);

        }
    }
}
Photon PhotonMap::createPhoton(const glm::vec3 point, const glm::vec3 &energy, float alpha, const glm::vec3 direction)
{
    Photon P;
    P.alpha = alpha;
    P.pos = point;
    P.energy = energy;
    P.dir = direction;
    return P;
}

