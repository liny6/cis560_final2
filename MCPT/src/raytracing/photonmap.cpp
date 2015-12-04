#include "photonmap.h"

glm::vec3 ComponentMult(const glm::vec3 &a, const glm::vec3 &b)
{
    return glm::vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}
float MIS(float f_PDF, float g_PDF)
 {
     return glm::pow((1.0f*f_PDF), 2.0f)/(glm::pow(1.0f*f_PDF, 2.0f) + glm::pow(1.0f*g_PDF, 2.0f));
 }

//~~~~~~~~Sorting  Functions ~~~~~~~~~~~//
bool zLessThan(const Photon &P1, const Photon &P2)
{
    return  P1.pos.z < P2.pos.z;

}

bool yLessThan (const Photon &P1, const Photon &P2)
{

  return  P1.pos.y < P2.pos.y;
}

bool xLessThan(const Photon &P1, const Photon &P2)
{
      return  P1.pos.x < P2.pos.x;

}

PhotonNode::PhotonNode()
{
   this->leftChild = NULL;
   this->rightChild = NULL;
//   photo_area = {0,0,0,0,0,0,0,0,0,0};
}

PhotonMap::PhotonMap():num_bounces(3),mersenne_generator(rand()), unif_distribution(0.0f,1.0f),max_photons(10000), num_photons(10000)
{
   scene = NULL;
   intersection_engine = NULL;
}
PhotonNode* PhotonMap::createIndirectPhotonMap(PhotonNode* root,QList<Photon*> &photon_list)
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
//~~~~~~Gathering Photons~~~~~~~//
PhotonNode* PhotonMap::gatherIndirectPhotons(QList<Photon> &indirect_photon_list, unsigned char XYZ) // return a node
{
    // find distance between nodes and have list sorted by distance from first item in list
    PhotonNode* node = new PhotonNode();
    for(int i = 0; i < PHOTONPERNODE; i++)
    {
        node->photo_area[i] =
    }

}

QList<Photon> PhotonMap::shootIndirectPhotons()
{
    QList<Photon> indirect_photons;
    Intersection light_isx;
    for (int i = 0; i < max_photons; i++)
    {
        light_isx = RandomSampleLight();
        placeIndirectPhoton(light_isx, indirect_photons);
    }
    return indirect_photons;

}
//~~~~~~Photon Level~~~~~~~~~~//
void PhotonMap::placeIndirectPhoton(const Intersection &light_isx, QList<Photon> &photon_list)
{
   // randomly choose geometry in the scene
    int geom_choice = 0;
    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    geom_choice = rand()%scene->objects.count();
    Geometry* scene_obj = scene->objects.at(geom_choice);
    Ray  photon_ray = scene_obj->GetRandRay(rand1,rand2, light_isx); // direction from light source to point on scene
   //~~~~~~~~~~~~~~~~~
   // using ray from random geometry isx and light isx find what the light actually hits
    Intersection obj_isx= intersection_engine->GetIntersection(photon_ray);
    Photon P;


   //DO RUSSIAN ROULETTE NEEDS CLEAN UP 2 ZA MAYZX
    float russian = unif_distribution(mersenne_generator);
    int depth = 0;
    light_isx.point = light_isx.point + 0.0001f; // avoid shadow acne
    float pdf_light_temp;
    float pdf_bxdf_temp;
    glm::vec3 alpha;
    glm::vec3 anew;
    alpha = glm::abs(glm::dot(light_isx.normal, photon_ray.direction) *L_energy)/(pdf_bxdf_temp*pdf_light_temp);
    float throughput;
    float throughput2;
   while ( depth < max_depth  )
   {
       // find max rgb value
       throughput = glm::max(glm::max(alpha.r, alpha.g), alpha.b);

       glm::vec3 L_energy = getPhotonEnergy(light_isx, obj_isx, photon_ray, wi_ret, pdf_bxdf_temp);
       pdf_light_temp = light_isx.object_hit->RayPDF(light_isx, photon_ray);

       P = createPhoton(obj_isx.point, L_temp, light_isx.object_hit->material->intensity, photon_ray.direction);
       photon_list.append(P);
       if(pdf_light_temp <=0.0f)
       {
           return;
       }
        anew = alpha*fr* glm::abs(glm::dot(wi_ret, obj.normal))/pdf_bxdf_temp;
        throughput2 = glm::max(glm::max(anew.r, anew.g), anew.b);
        float continueProb = glm::min(1.0f , throughput2/throughput);
        if( russian > continueProb)
        {
            return;
        }
        alpha = anew/ continueProb;

       //update sampler ray

       //find the next intersection and update accordingly

       wi_temp = -wo_temp;
        photon_ray = Ray(obj.point, wo_temp);
       obj_isx = intersection_engine->GetIntersection(photon_ray);

       depth++;
       russian = unif_distribution(mersenne_generator);
   }

}
Intersection PhotonMap::RandomSampleLight()
{
    int light_choice = 0;
    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    light_choice = rand()%scene->lights.count();
    Geometry* light_obj = scene->lights.at(light_choice);
    Intersection light_isx = light_obj->GetRandISX(rand1, rand2, glm::vec3(0));

    return light_isx;
}
glm::vec3 PhotonMap::getPhotonEnergy(const Intersection & light_isx, Intersection &obj, Ray photon_ray, glm::vec3 & wi_ret ,float & pdf_bxdf)
{
    float rand3,rand4, W,pdf_light_temp;
     rand3 = unif_distribution(mersenne_generator);
     rand4 = unif_distribution(mersenne_generator);
    glm::vec3 L_temp;
    pdf_light_temp = light_isx.object_hit->RayPDF(light_isx, photon_ray);
    glm::vec3  brdf_energy_temp =obj.object_hit->material->SampleAndEvaluateScatteredEnergy(obj ,photon_ray.direction,wi_ret,
                                                                                 pdf_bxdf,rand3,rand4);
    W = MIS(pdf_bxdf,pdf_light_temp);
    L_temp = brdf_energy_temp;
    L_temp = ComponentMult(ComponentMult(L_temp, obj.object_hit->material->base_color), obj.texture_color);
    L_temp = L_temp*W/pdf_bxdf*glm::abs(glm::dot(obj.normal, wi_ret));

   return L_temp;
}
Photon PhotonMap::createPhoton(const glm::vec3 point, const glm::vec3 &energy, float alpha, const glm::vec3 direction)
{
    Photon P;
    P.alpha = alpha/ num_photons;
    P.pos = point;
    P.energy = energy;
    P.dir = direction;
    return P;
}
bool PhotonMap::isCaustic(const Geometry* scene_obj)
{
    if (scene_obj->material->bxdfs.first()->name == "lambert1")
    {
            return false;
    }
            return true;
}

