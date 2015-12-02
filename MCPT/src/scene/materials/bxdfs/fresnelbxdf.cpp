#include "fresnelbxdf.h"
#include <src/la.h>
float CosTheta(const glm::vec3 &wo)
{
    return -wo.z;
}
float AbsCosTheta(const glm::vec3 &wo)
{
    return glm::abs(wo.z);
}

float SinTheta2(const glm::vec3 &wo)
{
    return glm::max(0.0f, (1.0f - CosTheta(wo)*CosTheta(wo)));
}
float FresDieletric(float cos_i, float cos_t, float ni, float nt)
{

    float rl = (nt*cos_i - ni*cos_t) / (nt*cos_i + ni*cos_t);
    float rp = (ni*cos_i - nt*cos_t) / (ni*cos_i + nt*cos_t);

    float F = (rl*rl + rp*rp)/2.0f;
    return F;

}

float FresnalBxDF::Fresnel(float cosi, float &nt_ret, float &ni_ret) const
{
    //calculates the fresnel term
    glm::vec3 N(0.0f, 0.0f, 1.0f);
    float cos_i =glm::clamp( cosi, -1.0f, 1.0f);
//    float cos_o = glm::dot(wo, N);
     ni_ret  = refract_idx_in;
     nt_ret = refract_idx_out;
    bool entering = cos_i > 0.0f;
    if (!entering)
    {
        nt_ret = refract_idx_in;
        ni_ret = refract_idx_out;
    }
//    qDebug()<<"ni value is "<< ni_ret;
    float sin_t = ni_ret/nt_ret * glm::sqrt(glm::max(0.0f, 1.0f - cos_i*cos_i));
    if (sin_t >1.0f)
    {
        return 1.0f;
    }
    float cos_t = glm::sqrt(1 - sin_t*sin_t);
    return FresDieletric(glm::abs(cos_i), cos_t, ni_ret, nt_ret);
}

glm::vec3 FresnalBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi) const
{
//    //TODO
//    //find half vector
    glm::vec3 H(glm::normalize(wi+wo));
//    //define local normal
    glm::vec3 N(0.0f, 0.0f, 1.0f);
//    //find G
//    float G1 = 1.0f;
//    float G2 = 2.0f * glm::dot(N, H) * glm::dot(N, wo) / glm::dot(wo, H);
//    float G3 = 2.0f * glm::dot(N, H) * glm::dot(N, wi) / glm::dot(wi, H);
//    float G = glm::min(G3, glm::min(G1, G2));
//    //use blinn phong distribution
//    float D = (exponent+2.0f)/2.0f/PI*glm::pow(glm::dot(N, H), exponent);
//    //calculate fresnel term
  //    glm::vec3 result(reflection_color*D*G*F/4.0f/glm::dot(N, wi)/glm::dot(N, wo));

    float ni_ret, nt_ret;
    float F = Fresnel(CosTheta(wo), nt_ret, ni_ret);
    float T =1.0f;
    glm::vec3 result =reflection_color*(nt_ret*nt_ret )/(ni_ret*ni_ret) *( 1.0f-F)*T/ AbsCosTheta(wi);
    return result;
}
glm::vec3 FresnalBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
//    //only sample the brdf lobe
      glm::vec3 N(0.0f, 0.0f, 1.0f);
//    glm::vec3 energy(0);
//    glm::vec3 lobe_vec(0);
//    glm::vec3 H_sample;

//    //this is samping the lobe in lobe space
//    float cos_theta = glm::pow(rand1, 1.0f/(exponent+1.0f));
//    float sin_theta = glm::sqrt(1.0f - cos_theta*cos_theta);
//    float phi = 2.0f*PI*rand2;
//    //lobe is reflected about the normal
//    lobe_vec.x = -wo.x;
//    lobe_vec.y = -wo.y;
//    lobe_vec.z = wo.z;
//    //the incident ray need to be converted back into tangent space
//    H_sample.x = glm::cos(phi)*sin_theta;
//    H_sample.y = glm::sin(phi)*sin_theta;
//    H_sample.z = cos_theta;
//    //check for half vectors that goes into the page
//    if(glm::dot(H_sample, N) <= 0.0f) H_sample = -H_sample;
//    //get wi given wo and H_sample

//    wi_ret = -wo + 2.0f * glm::dot(wo, H_sample) * H_sample;


    float sini2 = SinTheta2(wo);
   float eta = refract_idx_in/refract_idx_out;
   float sint_2 = eta*eta*sini2;
   float cos_t = glm::sqrt(glm::max(0.0f,1.0f -sint_2));
//   float cos_t = glm::sqrt(1.0f - glm::max(0.0f,1.0f -sint_2));

   bool enter = CosTheta(wo) > 0.0f;
   if(enter)
   {
       cos_t = -cos_t;
   }
//     wi_ret =  glm::refract(wo,N,eta);
    wi_ret  = glm::vec3(eta* -wo.x, eta*-wo.y, cos_t);
    pdf_ret = 1.0f;


    if (sint_2 >= 1.0f)
    {
       return glm::vec3(0);
    }

    glm::vec3 energy = EvaluateScatteredEnergy(wo, wi_ret);

    return energy;
}

glm::vec3 FresnalBxDF::EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const
{
    //TODO
    return glm::vec3(0);
}

float FresnalBxDF::PDF(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    glm::vec3 H(glm::normalize(wi+wo));
    glm::vec3 N(0.0f, 0.0f, 1.0f);
    float costheta = glm::dot(H, N);

    if(glm::dot(wo, H) <= 0.0f) return 0.0f;

    float blinnpdf = (exponent + 1.0f) * glm::pow(costheta, exponent) / (2.0f*PI*4.0f*glm::dot(wo, H));

    return blinnpdf;
}
