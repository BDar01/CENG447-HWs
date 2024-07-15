#version 120
#extension GL_EXT_gpu_shader4 : enable
vec3 lightPos = vec3(5, 5, 5);
vec3 eyePos = vec3(0, 0, 0);

vec3 I = vec3(2, 2, 2);
vec3 Iamb = vec3(0.8, 0.8, 0.8);

vec3 kd = vec3(0.2, 0, 0.7);
vec3 ka = vec3(0.1, 0.1, 0.1);
vec3 ks = vec3(0.8, 0.8, 0.8);

varying vec4 fragPos;
varying vec3 N;

// Time uniform
uniform float currentOffset;

bool checkerboard(vec3 pos, float offset, float scale) {
    // Use time to modulate the position of the checkerboard pattern
    float movingOffset = offset - currentOffset;

    bool x = int((pos.x + offset) * scale) % 2 == 1;
    bool y = int((pos.y + offset) * scale) % 2 == 1;
    bool z = int((pos.z + movingOffset) * scale) % 2 == 1;
    bool xOrY = x != y;
    
    if (xOrY != z) {
        return true;
    } else {
        return false;
    }
}

void main(void)
{
    vec3 L = normalize(lightPos - vec3(fragPos));
    vec3 V = normalize(eyePos - vec3(fragPos));
    vec3 H = normalize(L + V);
    float NdotL = dot(N, L);
    float NdotH = dot(N, H);

    // Use the checkerboard pattern to determine color
    if (checkerboard(fragPos.xyz, -9.5, 0.205)) {
        gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0); // Black
    } else {
        vec3 diffuseColor = I * kd * max(0, NdotL);
        vec3 ambientColor = Iamb * ka;
        vec3 specularColor = I * ks * pow(max(0, NdotH), 20);

        gl_FragColor = vec4(diffuseColor + ambientColor + specularColor, 1);
    }
}
