#version 120

vec3 lightPos = vec3(5, 5, 5);
vec3 eyePos = vec3(0, 0, 0);

vec3 I = vec3(2, 2, 2);
vec3 Iamb = vec3(0.8, 0.8, 0.8);

vec3 kd = vec3(1.0, 1.0, 0.0); // Full intensity red and green, no blue
vec3 ka = vec3(0.1, 0.1, 0.0); // Similarly adjust the ambient color
// You can also adjust ks if you want to affect the specular color
vec3 ks = vec3(0.8, 0.8, 0.0); // Reducing the blue component in specular color

varying vec4 fragPos;
varying vec3 N;

void main(void)
{
	vec3 L = normalize(lightPos - vec3(fragPos));
	vec3 V = normalize(eyePos - vec3(fragPos));
	vec3 H = normalize(L + V);
	float NdotL = dot(N, L);
	float NdotH = dot(N, H);

	vec3 diffuseColor = I * kd * max(0, NdotL);
	vec3 ambientColor = Iamb * ka;
	vec3 specularColor = I * ks * pow(max(0, NdotH), 20);

    gl_FragColor = vec4(diffuseColor + ambientColor + specularColor, 1);
}
