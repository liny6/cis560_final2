#pragma once
#include <scene/materials/bxdfs/bxdf.h>

class FresnalBxDF : public BxDF
{
public:
//Constructors/Destructors
    FresnalBxDF() : FresnalBxDF(glm::vec3(0.5f), 1.0f)
    {}
    FresnalBxDF(const glm::vec3 &color) : FresnalBxDF(color, 1.0f)
    {}
    FresnalBxDF(const glm::vec3 &color, float exp) : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)), reflection_color(color), exponent(exp)
    {}
//Functions
    float Fresnel(float cosi, float &nt_ret, float &ni_ret) const;
    virtual glm::vec3 EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const;
    virtual glm::vec3 EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const;
    virtual float PDF(const glm::vec3 &wo, const glm::vec3 &wi) const;
    virtual glm::vec3 SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const;
//Member variables
    float refract_idx_in;
    float refract_idx_out;
    glm::vec3 reflection_color;
    float exponent;
};
