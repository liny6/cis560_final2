#include <raytracing/integrator.h>


Integrator::Integrator():
    max_depth(5),
    mersenne_generator(rand()),
    unif_distribution(0.0f, 1.0f)
{
    scene = NULL;
    intersection_engine = NULL;
}

glm::vec3 ComponentMult(const glm::vec3 &a, const glm::vec3 &b)
{
    return glm::vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

//Basic ray trace
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{
    //TODO
    /*Terminate if too deep*/
    if(depth >= max_depth) return glm::vec3(0, 0, 0);
    /*---------------------*/



    /*Variables We may need*/

    //for the first intersection with the scene and results
    Intersection isx; // first intersection hit,
    glm::vec3 Emitted_Light(0,0,0); // stores the color value of emitted light from the point hit
    glm::vec3 Direct_Lighting(0,0,0); //stores the direct lighting
    float epsilon(0.001); //small distance;
    glm::vec3 Reflected_Light(0,0,0);// stores the color value of reflected light from the point hit
    unsigned int num_samples = 1; //number of samples per intersection

    //for BRDF PDF sampling
    //

    /*---------------------*/


    /*Light Source sampling*/

    //first find out who and where I hit
    isx = intersection_engine->GetIntersection(r);
    if(isx.object_hit == NULL)
    {
        return glm::vec3(0, 0, 0);
    }
    //traverse back the normal to prevent shadow acne
    isx.point = isx.point + epsilon*isx.normal;
    //pass in isx into EstimateDirectLighting
    //Direct_Lighting = EstimateDirectLighting(isx, num_samples, -r.direction);
    //if light source, add my light value to it
    if(isx.object_hit->material->is_light_source)
    {Emitted_Light = isx.object_hit->material->base_color;}
    else
    {Direct_Lighting = EstimateDirectLighting(isx, num_samples, -r.direction);}
    /*---------------------*/

    return Emitted_Light + Direct_Lighting;
}

glm::vec3 Integrator::EstimateDirectLighting(const Intersection &isx, unsigned int &samples_taken, const glm::vec3 &woW)
{
    //for Light source sampling
    //QList<glm::vec3> light_sample_pts;
    glm::vec3 color_temp(0,0,0);
    glm::vec3 color_final(0,0,0);
    Intersection light_sample_isx; //randomly sampled intersection on the light source's surface
    glm::vec3 wiW; // incoming ray in world frame
    Ray light_sample_ray;
    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);

    Intersection obstruction_test;

    //iterate through all the light sources
    for (int i = 0; i < scene->lights.count(); i++)
    {
        for(unsigned int j = 0; j < samples_taken; j++)
        {
            light_sample_isx = scene->lights[i]->GetRandISX(rand1, rand2, isx.normal); //take 1 sample point(intersection) on the light source for now
            wiW = light_sample_isx.point - isx.point; //ray direction going from world point to light source
            light_sample_isx.t = glm::length(wiW);
            wiW = glm::normalize(wiW);
            light_sample_ray = Ray(isx.point, wiW);//remember, the direction is from point in scene to light source
            obstruction_test = intersection_engine->GetIntersection(light_sample_ray);
            //update random point
            rand1 = unif_distribution(mersenne_generator);
            rand2 = unif_distribution(mersenne_generator);

            if (obstruction_test.object_hit == scene->lights[i])
            {
                color_temp = color_temp + CalculateEnergy(light_sample_isx, isx, light_sample_ray, woW);
            }
            else
            {
                //the ray contributes zero energy
            }
        }
        color_temp = color_temp/static_cast<float>(samples_taken); //divide by samples taken
        color_final = color_final + color_temp; // accumulate energy per high source
        color_temp = glm::vec3(0, 0, 0);// zero out color_temp for the next light source
    }

    return color_final;
}

glm::vec3 Integrator::CalculateEnergy(const Intersection &light_sample_isx, const Intersection &isx, const Ray &light_sample, const glm::vec3 &woW)
{
    glm::vec3 ray_color(0, 0, 0);
    Material* M = isx.object_hit->material; //material of point hit
    Geometry* L = light_sample_isx.object_hit; //light source
    float dummy;
     //Intersection isx_light = L->GetIntersection(light_sample);

    ray_color = ComponentMult(L->material->EvaluateScatteredEnergy(light_sample_isx, woW, -light_sample.direction, dummy), M->EvaluateScatteredEnergy(isx, woW, light_sample.direction, dummy)); // multiply the energy of the light with BRDF reflected energy
    ray_color = ray_color/L->RayPDF(light_sample_isx, light_sample)*glm::abs(glm::dot(isx.normal, light_sample.direction)); // and then do the solid angle PDF and the cosine
    ray_color = glm::clamp(ray_color, 0.0f, 1.0f);
    return ray_color;
}

void Integrator::SetDepth(unsigned int depth)
{
    max_depth = depth;
}


glm::vec3 DirectLightingIntegrator::TraceRay(Ray r, unsigned int depth)
{
    //TODO
    //Terminate if too deep//
    if(depth >= max_depth) return glm::vec3(0, 0, 0);
    //---------------------//



    //Variables We may need//

    //for the first intersection with the scene and results
    Intersection isx; // first intersection hit,
    glm::vec3 Emitted_Light(0,0,0); // stores the color value of emitted light from the point hit
    glm::vec3 Direct_Lighting(0,0,0); //stores the direct lighting
    float epsilon(0.001); //small distance;
    unsigned int n_light = 100; //number of light samples
    unsigned int n_brdf = 100; //number of brdf samples

    //for BRDF PDF sampling
    //

    //---------------------//


    //Light Source sampling//

    //first find out who and where I hit
    isx = intersection_engine->GetIntersection(r);

    if(isx.object_hit == NULL)
    {
        return glm::vec3(0, 0, 0);
    }
    if(isx.object_hit->name == "Glossy_Plane_exp50")
    {
        int breakhere = 0;
    }
    //traverse back the normal to prevent shadow acne
    isx.point = isx.point + epsilon*isx.normal;
    //pass in isx into EstimateDirectLighting
    //Direct_Lighting = EstimateDirectLighting(isx, num_samples, -r.direction);
    //if light source, add my light value to it
    if(isx.object_hit->material->is_light_source)
    {Emitted_Light = isx.object_hit->material->base_color;}
    else
    {Direct_Lighting = EstimateDirectLighting(isx, n_light, n_brdf, -r.direction);}
    //---------------------//

    return Emitted_Light + Direct_Lighting;
}

glm::vec3 DirectLightingIntegrator::EstimateDirectLighting(const Intersection &isx, unsigned int &n_light, unsigned int &n_brdf, const glm::vec3 &woW)
{
    //for Light source sampling
    //QList<glm::vec3> light_sample_pts;
    glm::vec3 brdf_sampling_final(0,0,0);
    glm::vec3 light_sampling_final(0,0,0);
    glm::vec3 color_final(0,0,0);
    Intersection light_sample_isx; //randomly sampled intersection on the light source's surface
    glm::vec3 wiW; // incoming ray in world frame
    Ray light_sample_ray;
    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    int light_choice = 0;


    Intersection obstruction_test;

    //for all light sources
    //for (int i = 0; i < scene->lights.count(); i++)
    //{
    //light source sampling
    for(unsigned int j = 0; j < n_light; j++)
    {
        light_choice = rand()%scene->lights.count();
        light_sample_isx = scene->lights[light_choice]->GetRandISX(rand1, rand2, isx.normal); //take 1 sample point(intersection) on the light source for now
        wiW = light_sample_isx.point - isx.point; //ray direction going from world point to light source
        light_sample_isx.t = glm::length(wiW);
        wiW = glm::normalize(wiW);
        light_sample_ray = Ray(isx.point, wiW);//remember, the direction is from point in scene to light source
        obstruction_test = intersection_engine->GetIntersection(light_sample_ray);
        //update random point
        rand1 = unif_distribution(mersenne_generator);
        rand2 = unif_distribution(mersenne_generator);

        if (obstruction_test.object_hit == scene->lights[light_choice])
        {
            light_sampling_final = light_sampling_final + LightPDFEnergy(obstruction_test, isx, light_sample_ray, woW, n_light, n_brdf);
        }
        else
        {
            //the ray contributes zero energy
        }
    }


    light_sampling_final = light_sampling_final/static_cast<float>(n_light); //divide by samples taken
    //light_sampling_final = light_sampling_final + light_sampling_temp; // accumulate energy per high source
    //light_sampling_temp = glm::vec3(0, 0, 0);// zero out color_temp for the next light source
    //}

    //brdf sampling
    for(unsigned int j = 0; j < n_brdf; j++)
    {
        brdf_sampling_final = brdf_sampling_final + BxDFPDFEnergy(isx, woW, n_light, n_brdf);
    }
    brdf_sampling_final = brdf_sampling_final/static_cast<float>(n_brdf);
    color_final = brdf_sampling_final + light_sampling_final;
    return color_final;
}

glm::vec3 DirectLightingIntegrator::LightPDFEnergy(const Intersection &light_sample_isx, const Intersection &isx, const Ray &light_sample, const glm::vec3 &woW, unsigned int n_light, unsigned int n_brdf)
{
    glm::vec3 ray_color(0, 0, 0);
    Material* M = isx.object_hit->material; //material of point hit
    Geometry* L = light_sample_isx.object_hit; //light source
    float lightPDF = L->RayPDF(light_sample_isx, light_sample);
    //if light pdf is less than zero, return no light
    if (lightPDF <= 0.0f)
    {
        return glm::vec3(0);
    }
    //get BRDFPDF and energy from the material
    float brdfPDF;
    float dummy;

    glm::vec3 M_energy(M->EvaluateScatteredEnergy(isx, woW, light_sample.direction, brdfPDF));
    //terminate early if brdf pdf is zero;
    if (brdfPDF <= 0.0f) return ray_color;

    glm::vec3 L_energy(L->material->EvaluateScatteredEnergy(light_sample_isx, woW, -light_sample.direction, dummy));
    float W = MIS(lightPDF, brdfPDF); //MIS power heuristic weighing function

    ray_color = ComponentMult(L_energy, M_energy); // multiply the energy of the light with BRDF reflected energy
    ray_color = ComponentMult(ComponentMult(ray_color, M->base_color), isx.texture_color);
    ray_color = ray_color*W/lightPDF*glm::abs(glm::dot(isx.normal, light_sample.direction)); // and then do the solid angle PDF and the cosine
    return ray_color;
}

glm::vec3 DirectLightingIntegrator::BxDFPDFEnergy(const Intersection &isx, const glm::vec3 &woW, unsigned int n_light, unsigned int n_brdf)
{
    glm::vec3 ray_color(0, 0, 0);
    glm::vec3 wiW(0, 0, 0);//this will be obtained by sampling BxDf
    Material* M = isx.object_hit->material; //material of point hit
    float brdfPDF;
    float lightPDF;
    float dummy;
    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    glm::vec3 M_energy(M->SampleAndEvaluateScatteredEnergy(isx, woW, wiW, brdfPDF, rand1, rand2));
    //use sampled wiW to check if I can hit the light
    Ray shadow_feeler(isx.point, wiW);
    Intersection light_isx = intersection_engine->GetIntersection(shadow_feeler);
    Geometry* L; //this holds the intersected light source


    //terminate early if brdf pdf is zero;
    if (brdfPDF <= 0.0f) return ray_color;
    if (light_isx.object_hit == NULL)//if ray didnt hit anything
        return ray_color;
    if (light_isx.object_hit->material->is_light_source)
    {
        L = light_isx.object_hit;
        lightPDF = L->RayPDF(light_isx, shadow_feeler);
        if (lightPDF <= 0)
        {
            return ray_color;
        }
        glm::vec3 L_energy(L->material->EvaluateScatteredEnergy(light_isx, woW, -shadow_feeler.direction, dummy));
        float W = MIS(brdfPDF, lightPDF);

        ray_color = ComponentMult(L_energy, M_energy);
        ray_color = ComponentMult(ComponentMult(ray_color, M->base_color), isx.texture_color);
        ray_color = ray_color*W/brdfPDF*glm::abs(glm::dot(isx.normal, shadow_feeler.direction));
        return ray_color;
    }
    else
    {
        return ray_color;
    }

}

float DirectLightingIntegrator::MIS(float f_PDF, float g_PDF)
{
    return glm::pow((1.0f*f_PDF), 2.0f)/(glm::pow(1.0f*f_PDF, 2.0f) + glm::pow(1.0f*g_PDF, 2.0f));
}

glm::vec3 AllLightingIntegrator::TraceRay(Ray r, unsigned int depth)
{
    //Terminate if too deep//
    if(depth >= max_depth) return glm::vec3(0, 0, 0);
    //---------------------//



    //Variables We may need//

    //for the first intersection with the scene and results
    Intersection isx; // first intersection hit,
    glm::vec3 Emitted_Light(0,0,0); // stores the color value of emitted light from the point hit
    glm::vec3 Direct_Lighting(0,0,0); //stores the direct lighting
    float epsilon(0.001); //small distance;
    glm::vec3 Indirect_Lighting(0,0,0);// stores the color value of reflected light from the point hit
    unsigned int n_light = 10; //number of light samples
    unsigned int n_brdf = 10; //number of brdf samples
    unsigned int n_indirect = 10; //number of sample splits

    //for BRDF PDF sampling
    //

    //---------------------//


    //Light Source sampling//

    //first find out who and where I hit
    isx = intersection_engine->GetIntersection(r);

    if(isx.object_hit == NULL)
    {
        return glm::vec3(0, 0, 0);
    }

    //traverse back the normal to prevent shadow acne
    isx.point = isx.point + epsilon*isx.normal;
    //pass in isx into EstimateDirectLighting
    //Direct_Lighting = EstimateDirectLighting(isx, num_samples, -r.direction);
    //if light source, add my light value to it
    if(isx.object_hit->material->is_light_source)
    {Emitted_Light = isx.object_hit->material->base_color;}
    else
    {
        Direct_Lighting = EstimateDirectLighting(isx, n_light, n_brdf, -r.direction);
        Indirect_Lighting = BiDirIndirectEnergy(isx, n_indirect, -r.direction);
    }
    //---------------------//

    return Emitted_Light + Direct_Lighting + Indirect_Lighting;
}

glm::vec3 AllLightingIntegrator::EstimateIndirectLighting(const Intersection &isx, const unsigned int &n_split, const glm::vec3 &woW)
{
    //wrapper that estimates indirect lighting
    glm::vec3 color_accum(0.0f, 0.0f, 0.0f); //accumulated color
    glm::vec3 color_temp_light(0.0f, 0.0f, 0.0f);
    glm::vec3 color_temp_bxdf(0.0f, 0.0f, 0.0f);

    //Light Source "objects in scene" Sampling

    for(int i = 0; i < n_split; i++)
    {
        color_temp_light = LightIndirectEnergy(isx, n_split, woW);
        color_temp_bxdf = BxDFIndirectEnergy(isx, n_split, woW);
        color_accum = color_accum + color_temp_light + color_temp_bxdf;
    }

    return color_accum;
}

glm::vec3 AllLightingIntegrator::BxDFIndirectEnergy(const Intersection &isx, unsigned int n_split, const glm::vec3 &woW)
{
    int depth = 0;
    Intersection isx_temp = isx;//reflected intersection
    Intersection isx_light; //sampled intersection with "light source"
    glm::vec3 color_accum(0.0f, 0.0f, 0.0f); //accumulated color
    glm::vec3 color_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wi_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wo_temp(woW);
    glm::vec3 brdf_energy_accum(1.0f, 1.0f, 1.0f);
    glm::vec3 brdf_energy_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 L_temp(0.0f, 0.0f, 0.0f);
    Material* M_temp = isx.object_hit->material;
    Ray sampler;
    float pdf_temp_brdf(0);
    float pdf_light_temp(0);
    float W(0); //for MIS


    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    float epsilon = 0.0001;
    float throughput = 1.000001f;
    float russian = unif_distribution(mersenne_generator);
    //default samples for direct lighting when estimating the irradiance of some other point in the scene
    unsigned int n_light = 10;
    unsigned int n_brdf = 10;

    while(depth < max_depth && (russian < throughput || depth < 2))
    {
        //sample random brdf starting at the input isx direction to get reflected ray to begin
        brdf_energy_temp = M_temp->SampleAndEvaluateScatteredEnergy(isx_temp, wo_temp, wi_temp, pdf_temp_brdf, rand1, rand2);
        //update accumulated brdf energy
        brdf_energy_accum = ComponentMult(brdf_energy_accum, brdf_energy_temp);
        //use the sampled incoming ray to find a reflected intersection
        sampler = Ray(isx_temp.point, wi_temp);
        isx_light = intersection_engine->GetIntersection(sampler);

        //this ray dies if I hit a real light source or nothing
        if(isx_light.object_hit == NULL) break;
        else if(isx_light.object_hit->material->is_light_source) break;
        //

        //to avoid shadow acne
        isx_light.point = isx_light.point + epsilon*isx_light.normal;
        //find the direct lighting irradiance of this point towards my original intersection
        wo_temp = -wi_temp; //now the old incoming ray is the outgoing ray for the new intersection
        L_temp = EstimateDirectLighting(isx_light, n_light, n_brdf, wo_temp); //the direct lighting towards the isx
        pdf_light_temp = isx_light.object_hit->RayPDF(isx_light, sampler);
        W = MIS(pdf_temp_brdf, pdf_light_temp);
        //this is BRDF sampling so use the illumination equation for BRDF sampling to accumulate color
        color_temp = ComponentMult(brdf_energy_accum, L_temp);
        color_temp = ComponentMult(ComponentMult(color_temp, M_temp->base_color), isx_temp.texture_color);
        color_temp = color_temp*W/pdf_temp_brdf*glm::abs(glm::dot(isx_temp.normal, wi_temp));
        color_accum = color_accum + color_temp/static_cast<float>(n_split);

        throughput = throughput * glm::max(glm::max(color_accum.r, color_accum.g), color_accum.b);
        //update the temporary material
        M_temp = isx_light.object_hit->material;
        isx_temp = isx_light;
        //update random number
        rand1 = unif_distribution(mersenne_generator);
        rand2 = unif_distribution(mersenne_generator);
        //update depth
        depth++;
        russian = unif_distribution(mersenne_generator);
    }

    return color_accum;
}

glm::vec3 AllLightingIntegrator::LightIndirectEnergy(const Intersection &isx, unsigned int n_split, const glm::vec3 &woW)
{
    int depth = 0;
    Intersection isx_temp = isx;//reflected intersection
    Intersection isx_light; //sampled intersection with "light source"
    glm::vec3 color_accum(0.0f, 0.0f, 0.0f); //accumulated color
    glm::vec3 color_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wi_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wo_temp(woW);
    glm::vec3 brdf_energy_accum(1.0f, 1.0f, 1.0f);
    glm::vec3 brdf_energy_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 L_temp(0.0f, 0.0f, 0.0f);
    Material* M_temp = isx.object_hit->material;
    Geometry* obj_temp; //stores the temporary sampled object
    Ray sampler;
    float pdf_temp_brdf(0);
    float pdf_light_temp(0);
    float W(0); //for MIS


    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    float epsilon = 0.0001f;
    float throughput = 1.000001f;
    float russian = unif_distribution(mersenne_generator);
    //default samples for direct lighting when estimating the irradiance of some other point in the scene
    unsigned int n_light = 10;
    unsigned int n_brdf = 10;

    int light_source_choice(0);

    while(depth < max_depth && (russian < throughput || depth < 2))
    {
        //sample a random point on a random object in the scene
        light_source_choice = rand()%scene->objects.count();
        obj_temp = scene->objects[light_source_choice];
        //if I hit a real light source, kill the ray
        if (scene->objects[light_source_choice]->material->is_light_source) break;

        isx_light = obj_temp->GetRandISX(rand1, rand2, isx_temp.normal);
        //update random numbers
        rand1 = unif_distribution(mersenne_generator);
        rand2 = unif_distribution(mersenne_generator);
        //make ray towards these points
        wi_temp = glm::normalize(isx_light.point - isx_temp.point);
        sampler = Ray(isx_temp.point, wi_temp);
        //update my light intersection
        isx_light = intersection_engine->GetIntersection(sampler);

        //this ray dies if it hit nothing or is blocked, kill the ray as well
        if(isx_light.object_hit == NULL) break;
        if(isx_light.object_hit != obj_temp) break;

        //to avoid shadow acne
        isx_light.point = isx_light.point + epsilon*isx_light.normal;

        //find out the pdf w/r/t light
        pdf_light_temp = obj_temp->RayPDF(isx_light, sampler);
        //if my pdf is negative, kill the ray as well
        if(pdf_light_temp <= 0) break;

        //update accumulated brdf energy
        brdf_energy_temp = M_temp->EvaluateScatteredEnergy(isx_temp, wo_temp, wi_temp, pdf_temp_brdf);
        brdf_energy_accum = ComponentMult(brdf_energy_accum, brdf_energy_temp);

        //find the direct lighting irradiance of this point towards my original intersection
        wo_temp = -wi_temp; //now the old incoming ray is the outgoing ray for the new intersection
        L_temp = EstimateDirectLighting(isx_light, n_light, n_brdf, wo_temp);
        W = MIS(pdf_light_temp, pdf_temp_brdf);
        //this is light source sampling so use the illumination equation for BRDF sampling to accumulate color
        color_temp = ComponentMult(brdf_energy_accum, L_temp);
        color_temp = ComponentMult(ComponentMult(color_temp, M_temp->base_color), isx_temp.texture_color);
        color_temp = color_temp*W/pdf_light_temp*glm::abs(glm::dot(isx_temp.normal, wi_temp));
        color_accum = color_accum + color_temp/static_cast<float>(n_split);

        throughput = throughput * glm::max(glm::max(color_accum.r, color_accum.g), color_accum.b);

        //update the temporary material
        M_temp = isx_light.object_hit->material;
        isx_temp = isx_light;
        //update random number
        //update depth
        depth++;
        russian = unif_distribution(mersenne_generator);
    }
    return color_accum;
}

glm::vec3 AllLightingIntegrator::BiDirIndirectEnergy(const Intersection &isx_camera_start, int splits, const glm::vec3 &woW)
{
    //first allocate 2 lists to store sampled intersections, ray directions and energy
    QList<Path> forward; //this starts from camera
    QList<Path> backward; //this starts from light
    glm::vec3 color_final(0);
    glm::vec3 color_temp(0);
    glm::vec3 wo_temp(0);
    glm::vec3 wi_temp(0);
    glm::vec3 brdf_light(0);
    glm::vec3 brdf_temp(0);
    glm::vec3 L_temp(0);
    glm::vec3 L_additional(0);//this is the radiance from backward sample
    Material* M_light;
    Material* M_temp;
    Intersection isx_light;
    Ray sampler;

    float epsilon(0.01f);
    float pdf_brdf_light(0);
    float pdf_brdf_temp(0);
    float pdf_light(0);
    unsigned int n_light = 10;
    unsigned int n_brdf = 10;
    float W(0);
    float total_paths = 0.0f;
    //call grow path to construct the lists
    grow_path(isx_camera_start, forward, backward, woW, max_depth/2, splits);
    /*
    //let's check the if the forward grow is right...
    if (forward.count() > 0)
    {
        return forward.last().energy_accum;
    }
    */



    total_paths  = static_cast<float>(forward.count()*backward.count());
    if (backward.count() == 0)
    {
        //if grow path fails...?
        std::cout<<"grow failed \n";
        return forward.last.energy_accum;
    }
    //for all permutations of the forward and backward paths
    for(int i = 0; i < forward.count(); i++)
    {
        for(int j = 0; j < backward.count(); j++)
        {
            //calculate the indirect lighting from backward list sample path to forward list sample path
            //remember, backward.dir is the wi direction for the backward list sample, and forward.dir is the wo for that sample
            //backward.energy is the luminance the backward sample receives from the light, and forward.energy is the irradiance the forward sample already has towards the camera along the path
            wo_temp = forward[i].dir;
            wi_temp = glm::normalize(backward[j].isx.point - forward[i].isx.point);
            M_light = backward[j].isx.object_hit->material;
            M_temp = forward[i].isx.object_hit->material;
            //let's connect the two paths, note, this is pretty much light source sampling
            sampler = Ray(forward[i].isx.point, wi_temp);
            //check for occlusion
            isx_light = intersection_engine->GetIntersection(sampler);
            if(isx_light.object_hit!=backward[j].isx.object_hit)
            {
                color_final = color_final + forward[i].energy_accum;
                //dont bother with other stuff
                continue;
            }
            else if(glm::distance(isx_light.point, backward[j].isx.point)>epsilon)
            {
                color_final = color_final + forward[i].energy_accum;
                //this means i hit the other side of the object
                continue;
            }
            else
            {
                //find the light pdf and check for negatives
                pdf_light = isx_light.object_hit->RayPDF(isx_light, sampler);
                if(pdf_light <= 0) continue;
                //evaluate the additional radiance and check brdf pdf for negatives
                brdf_light = M_light->EvaluateScatteredEnergy(backward[j].isx, -wi_temp, backward[j].dir, pdf_brdf_light);
                if(pdf_brdf_light <= 0) continue;
                L_additional = ComponentMult(backward[j].energy_accum, brdf_light);
                //evaluate the brdf at the forward sample
                brdf_temp = M_temp->EvaluateScatteredEnergy(forward[i].isx, wo_temp, wi_temp, pdf_brdf_temp);
                brdf_temp = ComponentMult(brdf_temp, forward[i].brdf_accum); //this gets me to the final brdf to the camera direction
                L_temp = EstimateDirectLighting(isx_light, n_light, n_brdf, -wi_temp) + L_additional;
                //connecting the lines is akin to light source sampling
                W = MIS(pdf_light, pdf_brdf_temp);
                color_temp = ComponentMult(L_temp, brdf_temp);
                color_temp = ComponentMult(ComponentMult(color_temp, M_temp->base_color), forward[i].isx.texture_color);
                color_temp = color_temp*W/pdf_light*glm::abs(glm::dot(forward[i].isx.normal, wi_temp));

                //sum shit up
                color_final = color_final + forward[i].energy_accum + color_temp;
            }
        }
    }
    color_final = color_final/total_paths;
    return color_final;

}

void AllLightingIntegrator::grow_path(const Intersection &isx_camera_start, QList<Path> &forward_ret, QList<Path> &backward_ret, const glm::vec3 &woW, int half_max_depth, int splits)
{
    //do number of split times of light source sampling and number of split times of brdf sampling

    //--------grow path using two different kinds of sampling
    for(int i = 0; i < splits; i++)
    {
        grow_Bxdf_sampling(isx_camera_start, forward_ret, woW, half_max_depth, splits);
        grow_light_sampling(isx_camera_start, forward_ret, woW, half_max_depth, splits);
        grow_from_light_source(backward_ret, half_max_depth);
    }
    return;
}

void AllLightingIntegrator::grow_Bxdf_sampling(const Intersection &isx_camera_start, QList<Path> &forward_ret, const glm::vec3 &woW, int half_max_depth, int splits)
{
    int depth = 0;
    Intersection isx_temp = isx_camera_start;//reflected intersection
    Intersection isx_light; //sampled intersection with "light source"
    glm::vec3 color_accum(0.0f, 0.0f, 0.0f); //accumulated color
    glm::vec3 color_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wi_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wo_temp(woW);
    glm::vec3 brdf_energy_accum(1.0f, 1.0f, 1.0f);
    glm::vec3 brdf_energy_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 L_temp(0.0f, 0.0f, 0.0f);
    Material* M_temp = isx_camera_start.object_hit->material;
    Ray sampler;
    Path path_sample;
    float pdf_temp_brdf(0);
    float pdf_light_temp(0);
    float W(0); //for MIS


    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    float epsilon = 0.0001;

    //default 10 samples for direct lighting when estimating the irradiance of some other point in the scene
    unsigned int n_light = 10;
    unsigned int n_brdf = 10;

    while(depth < half_max_depth)
    {

        //sample random brdf starting at the input isx direction to get reflected ray to begin
        brdf_energy_temp = M_temp->SampleAndEvaluateScatteredEnergy(isx_temp, wo_temp, wi_temp, pdf_temp_brdf, rand1, rand2);
        //update accumulated brdf energy
        brdf_energy_accum = ComponentMult(brdf_energy_accum, brdf_energy_temp);
        //now with everything updated, we can construct a sample path
        path_sample.isx = isx_temp;
        path_sample.brdf_accum = brdf_energy_accum;
        path_sample.dir = wo_temp;

        //use the sampled incoming ray to find a reflected intersection
        sampler = Ray(isx_temp.point, wi_temp);
        isx_light = intersection_engine->GetIntersection(sampler);
        //this ray dies if I hit a real light source or nothing
        if(isx_light.object_hit == NULL) break;
        else if(isx_light.object_hit->material->is_light_source) break;

        //to avoid shadow acne
        isx_light.point = isx_light.point + epsilon*isx_light.normal;
        //find the direct lighting irradiance of this point towards my original intersection
        wo_temp = -wi_temp; //now the old incoming ray is the outgoing ray for the new intersection
        L_temp = EstimateDirectLighting(isx_light, n_light, n_brdf, wo_temp); //the direct lighting towards the isx
        pdf_light_temp = isx_light.object_hit->RayPDF(isx_light, sampler);
        //if the light pdf is not possible, kill the ray
        if(pdf_light_temp <= 0) break;

        W = MIS(pdf_temp_brdf, pdf_light_temp);
        //this is BRDF sampling so use the illumination equation for BRDF sampling to accumulate color
        color_temp = ComponentMult(brdf_energy_accum, L_temp);
        color_temp = ComponentMult(ComponentMult(color_temp, M_temp->base_color), isx_temp.texture_color);
        color_temp = color_temp*W/pdf_temp_brdf*glm::abs(glm::dot(isx_temp.normal, wi_temp));
        color_accum = color_accum + color_temp;
        path_sample.energy_accum = color_accum;

        //append to forward list
        forward_ret.append(path_sample);

        //update the temporary material
        M_temp = isx_light.object_hit->material;
        isx_temp = isx_light;
        //update random number
        rand1 = unif_distribution(mersenne_generator);
        rand2 = unif_distribution(mersenne_generator);
        //update depth
        depth++;
    }
    return;
}

void AllLightingIntegrator::grow_light_sampling(const Intersection &isx_camera_start, QList<Path> &forward_ret, const glm::vec3 &woW, int half_max_depth, int splits)
{
    int depth = 0;
    Intersection isx_temp = isx_camera_start;//reflected intersection
    Intersection isx_light; //sampled intersection with "light source"
    glm::vec3 color_accum(0.0f, 0.0f, 0.0f); //accumulated color
    glm::vec3 color_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wi_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wo_temp(woW);
    glm::vec3 brdf_energy_accum(1.0f, 1.0f, 1.0f);
    glm::vec3 brdf_energy_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 L_temp(0.0f, 0.0f, 0.0f);
    Material* M_temp = isx_camera_start.object_hit->material;
    Geometry* obj_temp; //stores the temporary sampled object
    Ray sampler;
    Path path_sample;
    float pdf_temp_brdf(0);
    float pdf_light_temp(0);
    float W(0); //for MIS

    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    float epsilon = 0.0001f;

    //default 10 samples for direct lighting when estimating the irradiance of some other point in the scene
    unsigned int n_light = 10;
    unsigned int n_brdf = 10;

    int light_source_choice(0);

    while(depth < max_depth)
    {
        //sample a random point on a random object in the scene
        light_source_choice = rand()%scene->objects.count();
        obj_temp = scene->objects[light_source_choice];
        //if I hit a real light source, kill the ray
        if (scene->objects[light_source_choice]->material->is_light_source) break;

        isx_light = obj_temp->GetRandISX(rand1, rand2, isx_temp.normal);
        //update random numbers
        rand1 = unif_distribution(mersenne_generator);
        rand2 = unif_distribution(mersenne_generator);
        //make ray towards these points
        wi_temp = glm::normalize(isx_light.point - isx_temp.point);
        sampler = Ray(isx_temp.point, wi_temp);
        //update my light intersection
        isx_light = intersection_engine->GetIntersection(sampler);

        //this ray dies if it hit nothing or is blocked
        if(isx_light.object_hit == NULL) break;
        if(isx_light.object_hit != obj_temp) break;

        //to avoid shadow acne
        isx_light.point = isx_light.point + epsilon*isx_light.normal;

        //find out the pdf w/r/t light
        pdf_light_temp = obj_temp->RayPDF(isx_light, sampler);
        //if my pdf is negative, kill the ray as well
        if(pdf_light_temp <= 0) break;


        //update accumulated brdf energy
        brdf_energy_temp = M_temp->EvaluateScatteredEnergy(isx_temp, wo_temp, wi_temp, pdf_temp_brdf);
        brdf_energy_accum = ComponentMult(brdf_energy_accum, brdf_energy_temp);
        //if brdf pdf doesn't work, kill the ray
        if(pdf_temp_brdf <= 0) break;

        //update sampled path
        path_sample.isx = isx_temp;
        path_sample.dir = wo_temp;
        path_sample.brdf_accum = brdf_energy_accum;

        //find the direct lighting irradiance of this point towards my original intersection
        wo_temp = -wi_temp; //now the old incoming ray is the outgoing ray for the new intersection
        L_temp = EstimateDirectLighting(isx_light, n_light, n_brdf, wo_temp);
        W = MIS(pdf_light_temp, pdf_temp_brdf);
        //this is light source sampling so use the illumination equation for BRDF sampling to accumulate color
        color_temp = ComponentMult(brdf_energy_accum, L_temp);
        color_temp = ComponentMult(ComponentMult(color_temp, M_temp->base_color), isx_temp.texture_color);
        color_temp = color_temp*W/pdf_light_temp*glm::abs(glm::dot(isx_temp.normal, wi_temp));
        color_accum = color_accum;

        path_sample.energy_accum = color_accum;
        //append to the list
        forward_ret.append(path_sample);

        //update the temporary material
        M_temp = isx_light.object_hit->material;
        isx_temp = isx_light;
        //update random number
        //update depth
        depth++;
    }
    return;
}

void AllLightingIntegrator::grow_from_light_source(QList<Path> &backward_ret, int half_max_depth)
{
    int depth = 0;
    Intersection isx_object1;//
    Intersection isx_light; //sampled intersection on "light source"
    glm::vec3 color_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wi_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 wo_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 brdf_energy_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 L_temp(0.0f, 0.0f, 0.0f);
    glm::vec3 L_last(0.0f, 0.0f, 0.0f);
    Material* M_temp;
    Ray sampler;
    Path path_sample;
    float pdf_temp_brdf(0);
    float pdf_light_temp(0);
    float W(0); //for MIS

    float rand1 = unif_distribution(mersenne_generator);
    float rand2 = unif_distribution(mersenne_generator);
    float epsilon = 0.0001;

    int light_source_choice;
    unsigned int n_light = 10;
    unsigned int n_brdf = 10;


    //sample a random point on one of the light source
    light_source_choice = rand()%scene->lights.count();
    isx_light = scene->lights[light_source_choice]->GetRandISX(rand1, rand2, glm::vec3(0));
    //to avoid self intersection due to floating point error
    isx_light.point = isx_light.point + epsilon*isx_light.normal;
    //renew random numbers
    rand1 = unif_distribution(mersenne_generator);
    rand2 = unif_distribution(mersenne_generator);
    //sample some other object in the scene
    //randomlly shoot out rays from the light source from a hemisphere, this will be naive monte carlo sampling so don't bother with the pdf
    sampler = isx_light.object_hit->GetRandRay(rand1, rand2, isx_light);
    isx_object1 = intersection_engine->GetIntersection(sampler);
    if(isx_object1.object_hit == NULL) return;
    if(isx_object1.object_hit->material->is_light_source) return;

    M_temp = isx_object1.object_hit->material;

    //set up wi
    sampler.direction = -sampler.direction;
    wi_temp = sampler.direction;

    //update path and append to the list
    path_sample.isx = isx_object1;
    path_sample.energy_accum = L_temp;
    path_sample.dir = wi_temp;
    backward_ret.append(path_sample);


    while(depth < half_max_depth)
    {
        //renew random numbers
        rand1 = unif_distribution(mersenne_generator);
        rand2 = unif_distribution(mersenne_generator);
        //now the first intersection is set up, we can start traversing and caculate the energy a point receives from some direction
        //let's stick to brdf sampling for determining wos, since light sampling will likely be occluded
        //note that the wi and wo are reversed since i am doing things backwards, however, it's okay since brdf is the same when i swap the two directions
        brdf_energy_temp = M_temp->SampleAndEvaluateScatteredEnergy(isx_object1, wi_temp, wo_temp, pdf_temp_brdf, rand1, rand2);
        //find solid angle light pdf
        pdf_light_temp = isx_light.object_hit->RayPDF(isx_light, sampler);
        if(pdf_light_temp<=0)
        {
            break;
        }
        W = MIS(pdf_temp_brdf, pdf_light_temp);
        //calculate the energy object reflects towards the wo direction
        L_temp = EstimateDirectLighting(isx_object1, n_light, n_brdf, wo_temp) + L_last;
        L_temp = ComponentMult(brdf_energy_temp, L_temp);
        L_temp = ComponentMult(ComponentMult(color_temp, M_temp->base_color), isx_object1.texture_color);
        L_temp = color_temp*W/pdf_temp_brdf*glm::abs(glm::dot(isx_object1.normal, wi_temp));
        //update sampler ray
        sampler = Ray(isx_object1.point, wo_temp);
        //find the next intersection and update accordingly
        isx_object1 = intersection_engine->GetIntersection(sampler);
        M_temp = isx_object1.object_hit->material;
        //update wi
        wi_temp = -wo_temp;

        //update path and append to the list
        path_sample.isx = isx_object1;
        path_sample.energy_accum = L_temp;
        path_sample.dir = wi_temp;
        backward_ret.append(path_sample);

        L_last = L_temp;

        depth++;
    }
    return;
}
