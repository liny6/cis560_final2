#pragma once

#include <src/raytracing/intersection.h>
#include <src/la.h>
#include <scene/scene.h>
#include <src/raytracing/integrator.h>
#include <random>
#define  PHOTONAREA  100
#define PHOTONPERNODE 5
//Photon//
typedef struct __attribute__((__packed__))
{
    glm::vec3 pos;
    unsigned char XYZ; // KD splitting plane
    glm::vec3 dir;
    unsigned char photon_type; // {0 - direct, 1- indirect, 2 - whatever this thing is}
    glm::vec3 energy;
    float alpha; // Power of light source / number of emitted photons
} Photon;

//Photon node//
class PhotonNode
{
public:
    PhotonNode();
    PhotonNode* leftChild;
    PhotonNode* rightChild;
    Photon photo_area[PHOTONPERNODE];
//    QString name;

private:
};
//Creating Map//
class PhotonMap
{
public:
    // MAP LEVEL//
    PhotonMap();
    PhotonNode* createCausticPhotonMap();
    PhotonNode* createIndirectPhotonMap(); // need to allocate a butt load of memory
    // ~~~Gathering level~~~~~//
    void gatherIndirectPhotons(QList<Photon> indirect_photon_list);
    QList<Photon> shootIndirectPhotons();

    //~~~~~Photon Level~~~~//
    bool isCaustic(const Geometry* scene_obj);
    void placeIndirectPhoton(const Intersection &light_isx, QList<Photon> &photon_list);
    Intersection RandomSampleLight();
    glm::vec3 getPhotonEnergy(const Intersection & light_isx, Intersection &obj,
                                                Ray photon_ray, glm::vec3 & wi_ret ,float & pdf_bxdf);

    Photon createPhoton(const glm::vec3 point, const glm::vec3 &energy,
                                        float alpha, const glm::vec3 direction);

    // Housekeeping//
    void setMaxPhotons(int maxPhotons) {max_photons = maxPhotons;
                                       num_photons = (float)maxPhotons;}
    void setNumBounces(int max_bounce){num_bounces = max_bounce;}
    Scene* scene;
    IntersectionEngine* intersection_engine;
    void setLights(QList<Geometry*> L){lights =L;}
protected:
    QList<Geometry*> lights;
    std::mt19937 mersenne_generator;
    std::uniform_real_distribution<float> unif_distribution;
    int  max_photons;
    int num_bounces;
    float num_photons;

};

