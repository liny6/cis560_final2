#pragma once

#include <src/raytracing/intersection.h>
#include <src/la.h>
#include <scene/scene.h>
#include <src/raytracing/integrator.h>
#define  PHOTONAREA  10

typedef struct __attribute__((__packed__))
{
    glm::vec3 pos;
    unsigned char XYZ; // KD splitting plane
    glm::vec3 dir;
    unsigned char photon_type; // {0 - direct, 1- indirect, 2 - whatever this thing is}
    glm::vec3 energy;
    float alpha; // Power of light source / number of emitted photons
} Photon;


class PhotonNode
{
public:
    PhotonNode();
    PhotonNode* leftChild;
    PhotonNode* rightChild;
    Photon photo_area[PHOTONAREA];
    QString name;
//    Intersection GetIntersection(Ray r);
//    Intersection isx;
private:
};
class PhotonMap: public Integrator
{
public:
    PhotonMap();
    PhotonNode* createIndirectPhotonMap(); // need to allocate a butt load of memory
    PhotonNode* createDirectPhotonMap(PhotonNode* root, QList <Geometry*> &scene_geom, QList<Photon*> &photons );
    PhotonNode* createCausticPhotonMap();
    PhotonNode* createRadiancePhotonMap();
    PhotonNode* createVolumePhotonMap();
    void setMaxPhotons(int maxPhotons) {max_photons = maxPhotons; }
    void setNumBounces(int max_bounce){num_bounces = max_bounce;}
    Scene* scene;
    IntersectionEngine* intersection_engine;
    void setLights(QList<Geometry*> L){lights =L;}
private:
QList<Geometry*> lights;
int  max_photons;
int num_bounces;

};

