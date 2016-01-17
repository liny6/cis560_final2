#pragma once
#include<la.h>
#include<raytracing/intersection.h>

class Path
{
public:
    Path();
    void append(Intersection isx, glm::vec3 dir, glm::vec3 energy_accum);
private:
    QList<Intersection> isx_list;
    QList<glm::vec3*> dir_list;
    QList<glm::vec3*> energy_accum_list;
};
