#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <math.h>
#include <limits>
#include <algorithm>
#include <time.h>

typedef unsigned char RGB[3];

using namespace std;
using namespace parser;

int maxDepth;

float determinant(Vec3f a, Vec3f b, Vec3f c){
    float det = 0.0;
    det = a.x * (b.y*c.z - c.y*b.z) + a.y * (b.z*c.x - c.z*b.x) + a.z * (b.x*c.y - c.x*b.y);

    return det;
}

Vec3i convertV3ftoV3i(Vec3f v) {
    Vec3i i;
    i.x = int(v.x);
    i.y = int(v.y);
    i.z = int(v.z);

    return i;
}

Vec3f convertV3itoV3f(Vec3i i) {
    Vec3f v;
    v.x = float(i.x);
    v.y = float(i.y);
    v.z = float(i.z);

    return v;
}

Vec3f addV3f(Vec3f a, Vec3f b){
    Vec3f c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;

    return c;
}

Vec3f subV3f(Vec3f a, Vec3f b){
    Vec3f c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;

    return c;
}

Vec3f scalarMultV3f(float a, Vec3f b){
    Vec3f c;
    c.x = a * b.x;
    c.y = a * b.y;
    c.z = a * b.z;

    return c;
}

Vec3f vectorMultV3f(Vec3f a, Vec3f b) {
    Vec3f c;
    c.x = a.x * b.x;
    c.y = a.y * b.y;
    c.z = a.z * b.z;

    return c;   
}

Vec3f crossProdV3f(Vec3f a, Vec3f b){
    Vec3f c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;

    return c;
}

float dotProdV3f(Vec3f a, Vec3f b){
    float c = a.x * b.x + a.y * b.y + a.z * b.z;
    
    return c;
}

Vec3f normalizeVec3f(Vec3f v) {
    Vec3f normalized;
    float length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    normalized.x = v.x / length;
    normalized.y = v.y / length;
    normalized.z = v.z / length;

    return normalized;
}

float lengthVec3f(Vec3f v) {
    return sqrt(pow(v.x,2) + pow(v.y,2) + pow(v.z,2));
}

Vec3f reflectVec3f(Vec3f D, Vec3f N) {
    return subV3f(D, scalarMultV3f(2 * dotProdV3f(D, N), N));
}

unsigned char clampColor(int color) {
    unsigned char val;

    if (color > 255) {
        val = 255;
    }
    else if (color < 0) {
        val = 0;
    }
    else {
        val = (unsigned char) color;
    }

    return val;
}

Vec3f calcIrradiance(Vec3f &point, PointLight p) {
    Vec3f irr;

    Vec3f rad = subV3f(p.position, point);
    float rad2 = dotProdV3f(rad,rad);

    if (rad2 != 0.0) {
        irr = scalarMultV3f(1/rad2, p.intensity);
    }

    return irr;
}

Vec3i applyShading(Ray r, HitRecord hr, Scene &s) {
    Vec3f color = {0.0,0.0,0.0};

    Material m;
    if(hr.typeObject == "SPHERE"){
        m = s.materials[s.spheres[hr.objectID].material_id-1];
    }
    else if (hr.typeObject == "TRIANGLE"){
        m = s.materials[s.triangles[hr.objectID].material_id-1];
    }
    else if (hr.typeObject == "MESH"){
        m = s.materials[s.meshes[hr.objectID].material_id-1];
    }

    color = addV3f(color, vectorMultV3f(s.ambient_light, m.ambient));
    
    for (auto & light : s.point_lights) {

        Vec3f w_i = normalizeVec3f(subV3f(light.position, hr.point));

        Ray shadowRay;
        shadowRay.o = addV3f(hr.point, scalarMultV3f(s.shadow_ray_epsilon, w_i));
        shadowRay.d = w_i;

        HitRecord shadow_hr;

        float light_t = subV3f(light.position, shadowRay.o).x/shadowRay.d.x;

        if (!(closestHit(shadowRay, s, shadow_hr) && (shadow_hr.t < light_t && shadow_hr.t >= 0))) {
            
            Vec3f irr = calcIrradiance(hr.point,light);
            
            float dotProd = dotProdV3f(hr.normal, w_i);
            if (dotProd > 0) {
                Vec3f diffuseTerm = vectorMultV3f(scalarMultV3f(dotProd, m.diffuse), irr);
                color = addV3f(color, diffuseTerm);
            }

            Vec3f h = normalizeVec3f(subV3f(w_i, r.d));

            float dotProdSpec = dotProdV3f(hr.normal, h);
            
            if(dotProdSpec > 0) {
                Vec3f specularTerm = vectorMultV3f(scalarMultV3f(pow(dotProdSpec,m.phong_exponent), m.specular), irr);
                color = addV3f(color, specularTerm);
            }
        }
    }

    return convertV3ftoV3i(color);
}

void raySphereCalc(Sphere &sphere, int &sphereIndex, Ray &r, HitRecord &hr, Scene &s) {
    // (x – cx)^2 + (y – cy)^2 + (z – cz)^2 – R2=0
    float A = dotProdV3f(r.d, r.d);
    float B = 2*dotProdV3f(r.d, subV3f(r.o, s.vertex_data[sphere.center_vertex_id-1]));
    float C = dotProdV3f(subV3f(r.o, s.vertex_data[sphere.center_vertex_id-1]), subV3f(r.o, s.vertex_data[sphere.center_vertex_id-1])) - pow(sphere.radius,2);
    float discr = pow(B,2) - 4*A*C;

    if (discr < 0) {
        hr.hit = false;
    }
    else {
        float t1 = (-1 * B + sqrt(discr)) / 2*A;
        float t2 = (-1 * B - sqrt(discr)) / 2*A;
        if (t1 < 0 && t2 < 0) {
            hr.hit = false;
            hr.t = __FLT_MAX__;
        }
        else if (t1 < 0 && t2 >0){
            if (hr.hit == false) {
                hr.hit = true;
                hr.t = t2;
                hr.objectID = sphereIndex;
                hr.material_id = sphere.material_id;
                hr.typeObject = "SPHERE";
            }
            else {
                if (t2 < hr.t) {
                    hr.t = t2;
                    hr.objectID = sphereIndex;
                    hr.material_id = sphere.material_id;
                    hr.typeObject = "SPHERE";
                }
            }
        }
        else if (t1 >0 && t2 <0){
            if (hr.hit == false) {
                hr.hit = true;
                hr.t = t1;
                hr.objectID = sphereIndex;
                hr.material_id = sphere.material_id;
                hr.typeObject = "SPHERE";
            }
            else {
                if (t1 < hr.t) {
                    hr.t = t1;
                    hr.objectID = sphereIndex;
                    hr.material_id = sphere.material_id;
                    hr.typeObject = "SPHERE";
                }
            }
        }
        else if (t1 >0 && t2 >0){
            float mint=min(t1,t2);
            if (hr.hit == false) {
                hr.hit = true;
                hr.t = mint;
                hr.objectID = sphereIndex;
                hr.material_id = sphere.material_id;
                hr.typeObject = "SPHERE";
            }
            else {
                if (mint < hr.t) {
                    hr.t = mint;
                    hr.objectID = sphereIndex;
                    hr.material_id = sphere.material_id;
                    hr.typeObject = "SPHERE";
                }
            }

        }
        else {
            hr.hit = false;
        }
    }
}

float rayTriangleCalc(Face &f, Ray &r, Scene &s){
    float detA = determinant(subV3f(s.vertex_data[f.v0_id-1],s.vertex_data[f.v1_id-1]),
                        subV3f(s.vertex_data[f.v0_id-1],s.vertex_data[f.v2_id-1]), r.d);
    if(detA == 0){
        return __FLT_MAX__;
    }

    float t = determinant(subV3f(s.vertex_data[f.v0_id-1],s.vertex_data[f.v1_id-1]),
                        subV3f(s.vertex_data[f.v0_id-1],s.vertex_data[f.v2_id-1]),
                        subV3f(s.vertex_data[f.v0_id-1],r.o)) / detA;
    if(t <= 0){
        return __FLT_MAX__;
    }
    
    float gamma = determinant(subV3f(s.vertex_data[f.v0_id-1],s.vertex_data[f.v1_id-1]),
                        subV3f(s.vertex_data[f.v0_id-1],r.o), r.d) / detA;
    if(gamma < 0 || gamma > 1){
        return __FLT_MAX__;
    }

    float beta = determinant(subV3f(s.vertex_data[f.v0_id-1],r.o),
                        subV3f(s.vertex_data[f.v0_id-1],s.vertex_data[f.v2_id-1]), r.d) / detA;
    if(beta < 0 || beta > 1 - gamma){
        return __FLT_MAX__;
    }
    
    return t;
}

bool closestHit(Ray &r, Scene &s, HitRecord &hr) {
    hr.hit = false;
    hr.material_id = -1;
    hr.t = __FLT_MAX__;
    hr.typeObject = "none";
    hr.objectID = -1;
    vector<HitRecord> hits;

    for (int sphereIndex = 0; sphereIndex < s.spheres.size(); sphereIndex++){
        raySphereCalc(s.spheres[sphereIndex], sphereIndex, r, hr, s);
        if (hr.hit == true && hr.t >= 0 && hr.t < __FLT_MAX__) {
            hr.objectID = sphereIndex;
            hr.point = addV3f(r.o, scalarMultV3f(hr.t,r.d));
            hr.normal = scalarMultV3f(1/s.spheres[hr.objectID].radius, subV3f(hr.point, s.vertex_data[s.spheres[hr.objectID].center_vertex_id-1]));
            hits.push_back(hr);
        }
    }

    // iterate over triangle, function to check if point is inside triangle
    for (int triangleIndex = 0; triangleIndex < s.triangles.size(); triangleIndex++) {
        float t = rayTriangleCalc(s.triangles[triangleIndex].indices, r, s);

        if (t == __FLT_MAX__ ) {
            continue;
        }
        else {
            if (hr.hit == false ) {
                hr.hit = true;
                hr.t = t;
                hr.objectID = triangleIndex;
                hr.typeObject = "TRIANGLE";
                hr.material_id = s.triangles[hr.objectID].material_id;
                hr.normal = s.triangles[hr.objectID].indices.normal;
                hr.point = addV3f(r.o, scalarMultV3f(hr.t,r.d));
                hits.push_back(hr);
            }
            else {
                if (t < hr.t) {
                    hr.t = t;
                    hr.objectID = triangleIndex;
                    hr.typeObject = "TRIANGLE";
                    hr.normal = s.triangles[hr.objectID].indices.normal;
                    hr.material_id = s.triangles[hr.objectID].material_id;
                    hr.point = addV3f(r.o, scalarMultV3f(hr.t,r.d));
                    hits.push_back(hr);
                }
            }
        }
    }
    
    // iterate over mesh,use trangle function to check faces
    for (int meshIndex = 0; meshIndex < s.meshes.size(); meshIndex++) {
        for (int faceIndex = 0 ; faceIndex < s.meshes[meshIndex].faces.size(); faceIndex++) {
            float t = rayTriangleCalc(s.meshes[meshIndex].faces[faceIndex], r, s);
            if (t == __FLT_MAX__){
                continue;
            }
            else {
                if (hr.hit == false) {
                    hr.hit = true;
                    hr.t = t;
                    hr.objectID = meshIndex;
                    hr.typeObject = "MESH";
                    hr.faceID = faceIndex;
                    hr.material_id = s.meshes[hr.objectID].material_id;
                    hr.normal = s.meshes[hr.objectID].faces[hr.faceID].normal;
                    hr.point = addV3f(r.o, scalarMultV3f(hr.t,r.d));
                    hits.push_back(hr);
                }
                else {
                    if (t < hr.t) {
                        hr.t = t;
                        hr.objectID = meshIndex;
                        hr.typeObject = "MESH";
                        hr.faceID = faceIndex;
                        hr.material_id = s.meshes[hr.objectID].material_id;
                        hr.normal = s.meshes[hr.objectID].faces[hr.faceID].normal;
                        hr.point = addV3f(r.o, scalarMultV3f(hr.t,r.d));
                        hits.push_back(hr);
                    }
                }
            }
        }
    }
    //find closest hit in hits
    if (hits.size() > 0) {
        float minT = hits[0].t;
        int minIndex = 0;
        for (int i = 1; i < hits.size(); i++) {
            if (hits[i].t < minT) {
                minT = hits[i].t;
                minIndex = i;
            }
        }
        hr = hits[minIndex];
    }
    return hr.hit;
}

Vec3i computeColor(Ray &r, Scene &s) {
    Vec3f color;
    Vec3i def;
    def.x = 0;
    def.y = 0;
    def.z = 0;

    if (r.depth > maxDepth) {
        return def;
    }

    HitRecord hr;
    if (closestHit(r, s, hr)) {
        Vec3i localColor = applyShading(r, hr, s);
        color = convertV3itoV3f(localColor);

        Material m = s.materials[hr.material_id - 1];
        bool mrr = (lengthVec3f(m.mirror) > 0);

        if (r.depth < maxDepth && mrr) {
            Vec3f reflectionDir = reflectVec3f(r.d, hr.normal); 
            reflectionDir = normalizeVec3f(reflectionDir);
            Vec3f reflectionOrigin = addV3f(hr.point, scalarMultV3f(s.shadow_ray_epsilon, reflectionDir));

            Ray reflectionRay;
            reflectionRay.o = reflectionOrigin;
            reflectionRay.d = reflectionDir;
            reflectionRay.depth = r.depth + 1;

            Vec3i reflectedColor = computeColor(reflectionRay, s);
            Vec3f reflectionContribution = vectorMultV3f(convertV3itoV3f(reflectedColor), m.mirror);

            color = addV3f(color, reflectionContribution);
        }

        return convertV3ftoV3i(color);
    } else if (r.depth == 0) {
        return s.background_color; 
    } else {
        return def; 
    }
}



void calculateTriangleNormals(Scene &s) {
    for (int i = 0; i < s.triangles.size(); i++){  
        Vec3f v0 = s.vertex_data[s.triangles[i].indices.v0_id - 1];
        Vec3f v1 = s.vertex_data[s.triangles[i].indices.v1_id - 1];
        Vec3f v2 = s.vertex_data[s.triangles[i].indices.v2_id - 1];

        Vec3f edge1 = subV3f(v1, v0);
        Vec3f edge2 = subV3f(v2, v0);

        Vec3f normal = normalizeVec3f(crossProdV3f(edge1, edge2));

        s.triangles[i].indices.normal = normal;
    }
}

void calculateMeshNormals(Scene &s) {
    for (int i = 0; i < s.meshes.size(); i++){
        for (int j = 0; j < s.meshes[i].faces.size(); j++){ 
            Vec3f v0 = s.vertex_data[s.meshes[i].faces[j].v0_id - 1];
            Vec3f v1 = s.vertex_data[s.meshes[i].faces[j].v1_id - 1];
            Vec3f v2 = s.vertex_data[s.meshes[i].faces[j].v2_id - 1];

            Vec3f edge1 = subV3f(v1, v0);
            Vec3f edge2 = subV3f(v2, v0);

            Vec3f normal = normalizeVec3f(crossProdV3f(edge1, edge2));

            s.meshes[i].faces[j].normal = normal;
        }
    }
}

int main(int argc, char* argv[]) { 
    // Sample usage for reading an XML scene file
    Scene scene;
    scene.loadFromXml(argv[1]);
    
    calculateTriangleNormals(scene);
    calculateMeshNormals(scene);
    maxDepth = scene.max_recursion_depth;

    for (auto & camera : scene.cameras){       

        unsigned char* image = new unsigned char [camera.image_width * camera.image_height * 3];

        int bf = 0;
        Vec3f normalizedGaze = normalizeVec3f(camera.gaze);
        Vec3f m = addV3f(camera.position, scalarMultV3f(camera.near_distance, normalizedGaze));

        Vec3f u = crossProdV3f(normalizedGaze,camera.up);
        Vec3f normalizedU = normalizeVec3f(u);

        Vec3f v = crossProdV3f(normalizedU, normalizedGaze);

        Vec3f q = addV3f(m, scalarMultV3f(camera.near_plane.x, normalizedU));
        q = addV3f(scalarMultV3f(camera.near_plane.w, camera.up), q);
        

        for (int i = 0; i < camera.image_height; i++) {
            float s_v = (i + 0.5) * ((camera.near_plane.w - camera.near_plane.z) / (camera.image_height));

            for (int j = 0; j < camera.image_width; j++){
                float s_u = (j + 0.5) * ((camera.near_plane.y - camera.near_plane.x) / (camera.image_width));
                Vec3f s;
                s = addV3f(q, subV3f(scalarMultV3f(s_u, normalizedU), scalarMultV3f(s_v, v)));

                Ray r;  
                r.o = camera.position;
                r.d = subV3f(s, camera.position); 
                r.d = normalizeVec3f(r.d);

                Vec3i color = computeColor(r, scene);

                image[bf++] = clampColor(color.x);
                image[bf++] = clampColor(color.y);
                image[bf++] = clampColor(color.z);
            }
        }
        write_ppm(camera.image_name.c_str(), image, camera.image_width, camera.image_height);

    
    }
}