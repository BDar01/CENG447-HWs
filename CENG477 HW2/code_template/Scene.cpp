#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	//command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	//system(command.c_str());
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera &camera)
{
	for(auto &mesh : meshes) {
        
		//Perform modelling transformations
        Matrix4 M = applyModellingTransformations(mesh);

		//Perform backface culling
		if (cullingEnabled) {
			calculateMeshNormals(mesh);
			performBackfaceCulling(mesh, camera);
		}

		// Perform camera transformation
		Matrix4 C = performCameraTransformation(camera);


        // Perform projection transformation
		Matrix4 P = performProjectionTransformation(camera);

		
		Matrix4 PC = multiplyMatrixWithMatrix(P, C);
	
		Matrix4 MVP = multiplyMatrixWithMatrix(PC, M);

		// Perform viewport transformation
		Matrix4 V = performViewportTransformation(mesh, camera);
		mesh->V = V;
		

		applyTransformations(mesh, camera);
    }
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/
Matrix4 Scene::applyModellingTransformations(Mesh *mesh) {
    // Apply translation transformations to modelling matrix M
	Matrix4 M = getIdentityMatrix();

    for (int i = 0; i < mesh->numberOfTransformations; i++) {
		// Check each transformation in mesh and multiply with composite matrix M
		if(mesh->transformationTypes[i] == 't'){
			Translation* translation;
			for (int j = 0; j < translations.size(); j++) {
				if (mesh->transformationIds[i] == translations[j]->translationId) {
					translation = translations[j];
					break;
				}
			}
			M = applyTranslations(&M, translation);
		}

		else if(mesh->transformationTypes[i] == 's'){
			Scaling* scaling;
			for (int j = 0; j < scalings.size(); j++) {
				if (mesh->transformationIds[i] == scalings[j]->scalingId) {
					scaling = scalings[j];
					break;
				}
			}
			M = applyScalings(&M, scaling);
		}

		else if(mesh->transformationTypes[i] == 'r'){
			Rotation* rotation;
			for (int j = 0; j < scalings.size(); j++) {
				if (mesh->transformationIds[i] == rotations[j]->rotationId) {
					rotation = rotations[j];
					break;
				}
			}
			M = applyRotations(&M, rotation);
		}
    }

	mesh->M = M;

	return M;
}

void Scene::applyTransformations(Mesh* mesh, Camera& camera) {
	for (int n = 0; n < mesh->triangles.size(); n++) {
		if (!cullingEnabled || mesh->triangles[n].isVisible) {
			Vec4 triangle[3];
			for (int i = 0; i < 3; i++) {
				Vec3* v3 = vertices[mesh->triangles[n].vertexIds[i]-1];
				Vec4 v4 = Vec4(v3->x, v3->y, v3->z, 1.0);
				Vec4 tmp = multiplyMatrixWithVec4(mesh->MVP, v4);


				triangle[i] = Vec4(tmp.x, tmp.y, tmp.z, tmp.t, vertices[mesh->triangles[n].vertexIds[i]-1]->colorId);
			}
			performRasterization(triangle, mesh, camera);
		}
	}
}

Matrix4 Scene::applyTranslations(Matrix4* transform, Translation* translation){
	Matrix4 T = getIdentityMatrix();
	T.values[0][3] = translation->tx;
	T.values[1][3] = translation->ty;
	T.values[2][3] = translation->tz;

	return multiplyMatrixWithMatrix(T, *transform);
}

Matrix4 Scene::applyScalings(Matrix4* transform, Scaling* scaling) {
	Matrix4 S = getIdentityMatrix();
	S.values[0][0] = scaling->sx;
	S.values[1][1] = scaling->sy;
	S.values[2][2] = scaling->sz;

	return multiplyMatrixWithMatrix(S, *transform);
}

Matrix4 Scene::applyRotations(Matrix4* transform, Rotation* rotation) {
	Vec3 u = normalizeVec3(Vec3(rotation->ux,rotation->uy,rotation->uz));
	Matrix4 R = rotationMatrix(u, rotation->angle);

	return multiplyMatrixWithMatrix(R, *transform);
}

Matrix4 Scene::rotationMatrix(Vec3 u, double theta) {
	Matrix4 R = getIdentityMatrix();
	double t = (theta*M_PI)/180.0;
	double cosA = cos(t);
	double sinA = sin(t);

	R.values[0][0] = cosA + u.x*u.x*(1-cosA);
	R.values[0][1] = u.x*u.y*(1.0-cosA) - u.z*sinA;
	R.values[0][2] = u.x*u.z*(1.0-cosA) + u.y*sinA;
	R.values[0][3] = 0.0;

	R.values[1][0] = u.x*u.y*(1.0-cosA) + u.z*sinA;
	R.values[1][1] = cosA + u.y*u.y*(1.0-cosA);
	R.values[1][2] = u.z*u.y*(1.0-cosA) - u.x*sinA;
	R.values[1][3] = 0.0;

	R.values[2][0] = u.x*u.z*(1.0-cosA) - u.y*sinA;
	R.values[2][1] = u.z*u.y*(1.0-cosA) + u.x*sinA;
	R.values[2][2] = cosA + u.z*u.z*(1.0-cosA);
	R.values[2][3] = 0.0;

	R.values[3][0] = 0.0;
	R.values[3][1] = 0.0;
	R.values[3][2] = 0.0;
	R.values[3][3] = 1.0;

	return R;
}

void Scene::performBackfaceCulling(Mesh *mesh, Camera& camera) {
	for (auto & triangle : mesh->triangles) {
		Vec4 vertex = Vec4(vertices[triangle.vertexIds[0]-1]->x, vertices[triangle.vertexIds[0]-1]->y, vertices[triangle.vertexIds[0]-1]->z, 1.0);

		Vec4 v0 = multiplyMatrixWithVec4(mesh->M, vertex);
		Vec4 n = normalizeVec4(multiplyMatrixWithVec4(inverseTransposeMatrix(mesh->M), triangle.normal));

		Vec3 normal = Vec3(n.x, n.y, n.z);

		Vec3 v = normalizeVec3(subtractVec3(Vec3(v0.x, v0.y, v0.z), camera.position));

		if (dotProductVec3(v, normal) < 0.0) {
			triangle.isVisible = false;
		}
		else{
			triangle.isVisible = true;
		}
	}
}

Matrix4 Scene::performCameraTransformation(Camera &camera) {
	Matrix4 Cam;

	Cam.values[0][0] = camera.u.x;
	Cam.values[0][1] = camera.u.y;
	Cam.values[0][2] = camera.u.z;
	Cam.values[0][3] = -1.0*(camera.position.x*camera.u.x + camera.position.y*camera.u.y + camera.position.z*camera.u.z);

	Cam.values[1][0] = camera.v.x;
	Cam.values[1][1] = camera.v.y;
	Cam.values[1][2] = camera.v.z;
	Cam.values[1][3] = -1.0*(camera.position.x*camera.v.x + camera.position.y*camera.v.y + camera.position.z*camera.v.z);

	Cam.values[2][0] = camera.w.x;
	Cam.values[2][1] = camera.w.y;
	Cam.values[2][2] = camera.w.z;
	Cam.values[2][3] = -1.0*(camera.position.x*camera.w.x + camera.position.y*camera.w.y + camera.position.z*camera.w.z);

	Cam.values[3][0] = 0.0;
	Cam.values[3][1] = 0.0;
	Cam.values[3][2] = 0.0;
	Cam.values[3][3] = 1.0;

	return Cam;
}

Matrix4 Scene::performProjectionTransformation(Camera &camera) {
    Matrix4 Proj;

	if (camera.projectionType == 0) {
		Proj.values[0][0] = 2.0/(camera.right-camera.left);
		Proj.values[0][1] = 0.0;
		Proj.values[0][2] = 0.0;
		Proj.values[0][3] = -1.0*(camera.right+camera.left)/(camera.right-camera.left);

		Proj.values[1][0] = 0.0;
		Proj.values[1][1] = 2.0/(camera.top-camera.bottom);
		Proj.values[1][2] = 0.0;
		Proj.values[1][3] = -1.0*(camera.top+camera.bottom)/(camera.top-camera.bottom);

		Proj.values[2][0] = 0.0;
		Proj.values[2][1] = 0.0;
		Proj.values[2][2] = -2.0/(camera.far-camera.near);
		Proj.values[2][3] = -1.0*(camera.far+camera.near)/(camera.far-camera.near);

		Proj.values[3][0] = 0.0;
		Proj.values[3][1] = 0.0;
		Proj.values[3][2] = 0.0;
		Proj.values[3][3] = 1.0;
	}
	else if (camera.projectionType == 1) {

		Proj.values[0][0] = 2.0*camera.near/(camera.right-camera.left);
		Proj.values[0][1] = 0.0;
		Proj.values[0][2] = (camera.right+camera.left)/(camera.right-camera.left);
		Proj.values[0][3] = 0.0;

		Proj.values[1][0] = 0.0;
		Proj.values[1][1] = 2.0*camera.near/(camera.top-camera.bottom);
		Proj.values[1][2] = (camera.top+camera.bottom)/(camera.top-camera.bottom);
		Proj.values[1][3] = 0.0;

		Proj.values[2][0] = 0.0;
		Proj.values[2][1] = 0.0;
		Proj.values[2][2] = (-1.0)*(camera.far+camera.near)/(camera.far-camera.near);
		Proj.values[2][3] = (-2.0)*(camera.far*camera.near)/(camera.far-camera.near);

		Proj.values[3][0] = 0.0;
		Proj.values[3][1] = 0.0;
		Proj.values[3][2] = -1.0;
		Proj.values[3][3] = 0.0;
	}
    
	return Proj;

}	

Matrix4 Scene::performViewportTransformation(Mesh *mesh, Camera &camera) {
    Matrix4 VP;

	VP.values[0][0] = camera.horRes/2.0;
	VP.values[0][1] = 0.0;
	VP.values[0][2] = 0.0;
	VP.values[0][3] = (camera.horRes - 1.0)/2.0;

	VP.values[1][0] = 0.0;
	VP.values[1][1] = camera.verRes/2.0;
	VP.values[1][2] = 0.0;
	VP.values[1][3] = (camera.verRes - 1.0)/2.0;

	VP.values[2][0] = 0.0;
	VP.values[2][1] = 0.0;
	VP.values[2][2] = 0.5;
	VP.values[2][3] = 0.5;

	return VP;
}

void Scene::calculateMeshNormals(Mesh *mesh) {
	for (int j = 0; j < mesh->triangles.size(); j++){ 
		Vec3 v0 = *vertices[mesh->triangles[j].vertexIds[0] - 1];
		Vec3 v1 = *vertices[mesh->triangles[j].vertexIds[1] - 1];
		Vec3 v2 = *vertices[mesh->triangles[j].vertexIds[2] - 1];

		Vec3 edge1 = subtractVec3(v1, v0);
		Vec3 edge2 = subtractVec3(v2, v0);

		Vec3 n = normalizeVec3(crossProductVec3(edge1, edge2));
		Vec4 normal = Vec4(n.x, n.y, n.z, 0.0);

		mesh->triangles[j].normal = normal;
	}
}

Matrix4 Scene::multiplyMatrixWithMatrix(Matrix4 &m1, Matrix4 &m2) {
    Matrix4 result{};

    result.values[0][0] = m1.values[0][0] * m2.values[0][0] + m1.values[0][1] * m2.values[1][0] + m1.values[0][2] * m2.values[2][0] + m1.values[0][3] * m2.values[3][0];
    result.values[0][1] = m1.values[0][0] * m2.values[0][1] + m1.values[0][1] * m2.values[1][1] + m1.values[0][2] * m2.values[2][1] + m1.values[0][3] * m2.values[3][1];
    result.values[0][2] = m1.values[0][0] * m2.values[0][2] + m1.values[0][1] * m2.values[1][2] + m1.values[0][2] * m2.values[2][2] + m1.values[0][3] * m2.values[3][2];
    result.values[0][3] = m1.values[0][0] * m2.values[0][3] + m1.values[0][1] * m2.values[1][3] + m1.values[0][2] * m2.values[2][3] + m1.values[0][3] * m2.values[3][3];

    result.values[1][0] = m1.values[1][0] * m2.values[0][0] + m1.values[1][1] * m2.values[1][0] + m1.values[1][2] * m2.values[2][0] + m1.values[1][3] * m2.values[3][0];
    result.values[1][1] = m1.values[1][0] * m2.values[0][1] + m1.values[1][1] * m2.values[1][1] + m1.values[1][2] * m2.values[2][1] + m1.values[1][3] * m2.values[3][1];
    result.values[1][2] = m1.values[1][0] * m2.values[0][2] + m1.values[1][1] * m2.values[1][2] + m1.values[1][2] * m2.values[2][2] + m1.values[1][3] * m2.values[3][2];
    result.values[1][3] = m1.values[1][0] * m2.values[0][3] + m1.values[1][1] * m2.values[1][3] + m1.values[1][2] * m2.values[2][3] + m1.values[1][3] * m2.values[3][3];

    result.values[2][0] = m1.values[2][0] * m2.values[0][0] + m1.values[2][1] * m2.values[1][0] + m1.values[2][2] * m2.values[2][0] + m1.values[2][3] * m2.values[3][0];
    result.values[2][1] = m1.values[2][0] * m2.values[0][1] + m1.values[2][1] * m2.values[1][1] + m1.values[2][2] * m2.values[2][1] + m1.values[2][3] * m2.values[3][1];
    result.values[2][2] = m1.values[2][0] * m2.values[0][2] + m1.values[2][1] * m2.values[1][2] + m1.values[2][2] * m2.values[2][2] + m1.values[2][3] * m2.values[3][2];
    result.values[2][3] = m1.values[2][0] * m2.values[0][3] + m1.values[2][1] * m2.values[1][3] + m1.values[2][2] * m2.values[2][3] + m1.values[2][3] * m2.values[3][3];

    result.values[3][0] = m1.values[3][0] * m2.values[0][0] + m1.values[3][1] * m2.values[1][0] + m1.values[3][2] * m2.values[2][0] + m1.values[3][3] * m2.values[3][0];
    result.values[3][1] = m1.values[3][0] * m2.values[0][1] + m1.values[3][1] * m2.values[1][1] + m1.values[3][2] * m2.values[2][1] + m1.values[3][3] * m2.values[3][1];
    result.values[3][2] = m1.values[3][0] * m2.values[0][2] + m1.values[3][1] * m2.values[1][2] + m1.values[3][2] * m2.values[2][2] + m1.values[3][3] * m2.values[3][2];
    result.values[3][3] = m1.values[3][0] * m2.values[0][3] + m1.values[3][1] * m2.values[1][3] + m1.values[3][2] * m2.values[2][3] + m1.values[3][3] * m2.values[3][3];

    return result;
}

Matrix4 Scene::inverseTransposeMatrix(Matrix4& tform) {
    double det = tform.values[0][0] * (tform.values[1][1] * tform.values[2][2] - tform.values[1][2] * tform.values[2][1])
               - tform.values[0][1] * (tform.values[1][0] * tform.values[2][2] - tform.values[1][2] * tform.values[2][0])
               + tform.values[0][2] * (tform.values[1][0] * tform.values[2][1] - tform.values[1][1] * tform.values[2][0]);

    if (det == 0) {
        throw std::runtime_error("Matrix is not invertible");
    }

    double invDet = 1.0 / det;

    Matrix4 inv;
    inv.values[0][0] = invDet * (tform.values[1][1] * tform.values[2][2] - tform.values[1][2] * tform.values[2][1]);
    inv.values[0][1] = (-1.0f*invDet) * (tform.values[1][0] * tform.values[2][2] - tform.values[0][2] * tform.values[1][2]);
    inv.values[0][2] = invDet * (tform.values[1][0] * tform.values[2][1] - tform.values[2][0] * tform.values[1][1]);
    inv.values[1][0] = (-1.0f*invDet)  * (tform.values[0][1] * tform.values[2][2] - tform.values[2][1] * tform.values[0][2]);
    inv.values[1][1] = invDet * (tform.values[0][0] * tform.values[2][2] - tform.values[0][2] * tform.values[2][0]);
    inv.values[1][2] = (-1.0f*invDet) * (tform.values[0][0] * tform.values[2][1] - tform.values[2][0] * tform.values[0][1]);
    inv.values[2][0] = invDet * (tform.values[0][1] * tform.values[1][2] - tform.values[1][1] * tform.values[0][2]);
    inv.values[2][1] = (-1.0f*invDet) * (tform.values[0][0] * tform.values[1][2] - tform.values[1][0] * tform.values[0][2]);
    inv.values[2][2] = invDet * (tform.values[0][0] * tform.values[1][1] - tform.values[0][1] * tform.values[1][0]);

    return inv;
}

bool Scene::isInsideFrustum(const Vec4& v1, const Vec4& v2) {
	if (v1.t < 0.0 && v2.t < 0.0) {
		return false;
    }
	else {
		return true;
	}
}

bool Scene::isOutsideSamePlane(const Vec4& v1, const Vec4& v2) {
    if (v1.t < 0.0 && v2.t > 0.0){
		return true;
	}
	else {
		return false;
	}
}

bool Scene::vis(double num, double denom, double &te, double &to) {
    if (denom == 0) {
        return num >= 0;
    }

    double t = num / denom;

    if (denom < 0) {
        if (t > to) return false;
        te = std::max(t, te);
    }
    else if (denom > 0) {
        if (t < te) return false;
        to = std::min(t, to);
    }
    return true;
}

bool Scene::clipLinePers(Vec4& v1, Vec4& v2, Color& c1, Color& c2) {
    if(!isInsideFrustum(v1, v2)) {
		return false;
	}
	
	bool outside1 = isOutsideSamePlane(v1, v2);
	bool outside2 = isOutsideSamePlane(v2, v1);

	if (outside1 || outside2) {
		if (outside1) {
			v1 = v1 + ((v2 - v1) * ((0.00001 - v1.t) / (v2.t - v1.t)));
			c1 = *colorsOfVertices[v1.colorId];
			c1 += (c2 - c1) * ((0.00001 - v1.t) / (v2.t - v1.t));
		}
		else {
			c2 = *colorsOfVertices[v2.colorId];
			v2 = v1 + ((v2 - v1) * ((0.00001 - v1.t) / (v2.t - v1.t)));
			c2 += (c2 - c1) * ((0.00001 - v1.t) / (v2.t - v1.t));
		}
	}

	double te = 0.0;
	double to = 1.0;
	Vec4 df = v2 - v1;
	if (vis(v1.t - v1.x, df.x - df.t, te,to)) {
		if (vis(v1.t + v1.x, -df.x - df.t, te,to)) {
			if (vis(v1.t - v1.y, df.y - df.t, te,to)) {
				if (vis(v1.t + v1.y, -df.t - df.y, te,to)) {
					if (vis(v1.t - v1.z, df.z - df.t, te,to)) {
						if (vis(v1.t + v1.z, -df.t - df.z, te,to)) {
							if (to < 1) {
								v2 = v1 + (v2-v1) * to;
								c2 = c1 + (c2-c1) * to;
							}
							if (te > 0) {
								v1 = v1 + (v2-v1) * te;
								c1 = c1 + (c2-c1) * te;
							}
							
							// Handle edge case where w component becomes zero after interpolation
							if (v1.t == 0) {
								v1.t = 0.00001;
							}
							if (v2.t == 0) {
								v2.t = 0.00001;
							}
							
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}

bool Scene::clipLineOr(Vec4& v1, Vec4& v2, Color& c1, Color& c2){
	double te = 0.0;
	double to = 1.0;
	double dx = v2.x-v1.x;
	double dy = v2.y-v1.y;
	double dz = v2.z-v1.z;
	
	if(vis(v1.x+1, -dx, te, to)){
		if(vis(1-v1.x, dx, te, to)){
			if(vis(v1.y+1, -dy, te, to)){
				if(vis(1-v1.y, dy, te, to)){
					if(vis(v1.z+1, -dz, te, to)){
						if(vis(1-v1.z, dz, te, to)){
							c1 = *colorsOfVertices[v1.colorId];
							c2 = *colorsOfVertices[v2.colorId];

							Vec4 df = Vec4(dx, dy, dz, 0, -1);
							Color c_df = c2-c1;
							if(to < 1){
								v2 = v1 + (df * to);
								c2 = c1 + (c_df * to);
							}
							if(te > 0){
								v1 = v1 + (df * te);
								c1 = c1 + (c_df * te);
							}
							return true;
						}
					}
				}
			}
		}
	}

	return false;	

}

void Scene::rasterizeWireframeTriangle(Mesh *mesh, Vec4* triangle, Camera &camera) {
	Vec4 v0 = triangle[0];
	Vec4 v1 = triangle[1];
	Vec4 v2 = triangle[2];

	// Clipping
	Color color1, color2;

	Vec4 clipped1_V0 = v0;
	Vec4 clipped1_V1 = v1;

	if (clipLinePers(clipped1_V0, clipped1_V1, color1, color2)) {
		clipped1_V0.x = clipped1_V0.x/clipped1_V0.t;
		clipped1_V0.y = clipped1_V0.y/clipped1_V0.t;
		clipped1_V1.x = clipped1_V1.x/clipped1_V1.t;
		clipped1_V1.y = clipped1_V1.y/clipped1_V1.t;

		clipped1_V0 = multiplyMatrixWithVec4(mesh->V, clipped1_V0);
		clipped1_V1 = multiplyMatrixWithVec4(mesh->V, clipped1_V1);
		
		int pixel_x0 = max(0, min((int)(clipped1_V0.x+0.5), camera.horRes-1));
		int pixel_y0 = max(0, min((int)(clipped1_V0.y+0.5), camera.verRes-1));
		int pixel_x1 = max(0, min((int)(clipped1_V1.x+0.5), camera.horRes-1));
		int pixel_y1 = max(0, min((int)(clipped1_V1.y+0.5), camera.verRes-1));

		drawLine(pixel_x0, pixel_y0, pixel_x1, pixel_y1, &color1, &color2);
	}

	Vec4 clipped2_V1 = v1;
	Vec4 clipped2_V2 = v2;

	if (clipLinePers(clipped2_V1, clipped2_V2, color1, color2)) {
		clipped2_V1.x = clipped2_V1.x/clipped2_V1.t;
		clipped2_V1.y = clipped2_V1.y/clipped2_V1.t;
		clipped2_V2.x = clipped2_V2.x/clipped2_V2.t;
		clipped2_V2.y = clipped2_V2.y/clipped2_V2.t;

		clipped2_V1 = multiplyMatrixWithVec4(mesh->V, clipped2_V1);
		clipped2_V2 = multiplyMatrixWithVec4(mesh->V, clipped2_V2);
		
		int pixel_x1 = max(0, min((int)(clipped2_V1.x+0.5), camera.horRes-1));
		int pixel_y1 = max(0, min((int)(clipped2_V1.y+0.5), camera.verRes-1));
		int pixel_x2 = max(0, min((int)(clipped2_V2.x+0.5), camera.horRes-1));
		int pixel_y2 = max(0, min((int)(clipped2_V2.y+0.5), camera.verRes-1));
		
		drawLine(pixel_x1, pixel_y1, pixel_x2, pixel_y2, &color1, &color2);
	}

	Vec4 clipped3_V0 = v0;
	Vec4 clipped3_V2 = v2;

	if (clipLinePers(clipped3_V2, clipped3_V0, color1, color2)) {
		clipped3_V0.x = clipped3_V0.x/clipped3_V0.t;
		clipped3_V0.y = clipped3_V0.y/clipped3_V0.t;
		clipped3_V2.x = clipped3_V2.x/clipped3_V2.t;
		clipped3_V2.y = clipped3_V2.y/clipped3_V2.t;

		clipped3_V0 = multiplyMatrixWithVec4(mesh->V, clipped3_V0);
		clipped3_V2 = multiplyMatrixWithVec4(mesh->V, clipped3_V2);
		
		int pixel_x0 = max(0, min((int)(clipped3_V0.x+0.5), camera.horRes-1));
		int pixel_y0 = max(0, min((int)(clipped3_V0.y+0.5), camera.verRes-1));
		int pixel_x2 = max(0, min((int)(clipped3_V2.x+0.5), camera.horRes-1));
		int pixel_y2 = max(0, min((int)(clipped3_V2.y+0.5), camera.verRes-1));
		
		drawLine(pixel_x2, pixel_y2, pixel_x0, pixel_y0, &color1, &color2);
	}
}

void Scene::rasterizeWireframeTriangleOrtho(Mesh *mesh, Vec4* triangle, Camera &camera) {
	Vec4 v0 = triangle[0];
	Vec4 v1 = triangle[1];
	Vec4 v2 = triangle[2];

	// Clipping
	Color color1, color2;

	Vec4 clipped1_V0 = v0;
	Vec4 clipped1_V1 = v1;

	if (clipLineOr(clipped1_V0, clipped1_V1, color1, color2)) {

		clipped1_V0 = multiplyMatrixWithVec4(mesh->V, clipped1_V0);
		clipped1_V1 = multiplyMatrixWithVec4(mesh->V, clipped1_V1);
		
		int pixel_x0 = max(0, min((int)(clipped1_V0.x+0.5), camera.horRes-1));
		int pixel_y0 = max(0, min((int)(clipped1_V0.y+0.5), camera.verRes-1));
		int pixel_x1 = max(0, min((int)(clipped1_V1.x+0.5), camera.horRes-1));
		int pixel_y1 = max(0, min((int)(clipped1_V1.y+0.5), camera.verRes-1));

		drawLine(pixel_x0, pixel_y0, pixel_x1, pixel_y1, &color1, &color2);
	}

	Vec4 clipped2_V1 = v1;
	Vec4 clipped2_V2 = v2;

	if (clipLineOr(clipped2_V1, clipped2_V2, color1, color2)) {

		clipped2_V1 = multiplyMatrixWithVec4(mesh->V, clipped2_V1);
		clipped2_V2 = multiplyMatrixWithVec4(mesh->V, clipped2_V2);
		
		int pixel_x1 = max(0, min((int)(clipped2_V1.x+0.5), camera.horRes-1));
		int pixel_y1 = max(0, min((int)(clipped2_V1.y+0.5), camera.verRes-1));
		int pixel_x2 = max(0, min((int)(clipped2_V2.x+0.5), camera.horRes-1));
		int pixel_y2 = max(0, min((int)(clipped2_V2.y+0.5), camera.verRes-1));
		
		drawLine(pixel_x1, pixel_y1, pixel_x2, pixel_y2, &color1, &color2);
	}

	Vec4 clipped3_V0 = v0;
	Vec4 clipped3_V2 = v2;

	if (clipLineOr(clipped3_V2, clipped3_V0, color1, color2)) {

		clipped3_V0 = multiplyMatrixWithVec4(mesh->V, clipped3_V0);
		clipped3_V2 = multiplyMatrixWithVec4(mesh->V, clipped3_V2);
		
		int pixel_x0 = max(0, min((int)(clipped3_V0.x+0.5), camera.horRes-1));
		int pixel_y0 = max(0, min((int)(clipped3_V0.y+0.5), camera.verRes-1));
		int pixel_x2 = max(0, min((int)(clipped3_V2.x+0.5), camera.horRes-1));
		int pixel_y2 = max(0, min((int)(clipped3_V2.y+0.5), camera.verRes-1));
		
		drawLine(pixel_x2, pixel_y2, pixel_x0, pixel_y0, &color1, &color2);
	}
}
 
void Scene::drawLine(int x0, int y0, int x1, int y1, Color* c1, Color* c2){
    bool notSteep = abs(y1 - y0) < abs(x1 - x0);
	bool flag = false;

    if (!notSteep) {
        std::swap(x0, y0);
        std::swap(x1, y1);
    }

    if (x0 > x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
        std::swap(c1, c2);
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
	int yStep;

	if (dy < 0) {
		yStep = 1;
		dy = -dy;
	}
	else {
		yStep = 1;
	}

    int y_coord = y0;
	int error = 2*dy - dx;

    Color colorIncrement = (*c2 - *c1);
	colorIncrement.b /= dx;
	colorIncrement.g /= dx;
	colorIncrement.r /= dx;

    Color currentColor = *c1;

    for (int x_coord = x0; x_coord <= x1; x_coord++) {
        if (!notSteep) {
            image[y_coord][x_coord] = currentColor;  
        } else {
            image[x_coord][y_coord] = currentColor;
        }

        if (error > 0) {
            y_coord += yStep;
            error += 2*(dy-dx);
        }
		else {
			error += 2*dy;
		}

        currentColor += colorIncrement;
    }
}

void Scene::rasterizeSolidTriangle(Mesh *mesh, Vec4* triangle, Camera &camera){
	Vec4 v0 = triangle[0];
	Vec4 v1 = triangle[1];
	Vec4 v2 = triangle[2];

	v0.x /= v0.t;
	v0.y /= v0.t;
	v1.x /= v1.t;
	v1.y /= v1.t;
	v2.x /= v2.t;
	v2.y /= v2.t;

	v0 = multiplyMatrixWithVec4(mesh->V, v0);
	v1 = multiplyMatrixWithVec4(mesh->V, v1);
	v2 = multiplyMatrixWithVec4(mesh->V, v2); 

	int pixel_x0 = max(0, min(camera.horRes-1, (int)round(v0.x)));
	int pixel_y0 = max(0, min(camera.verRes-1, (int)round(v0.y)));
	int pixel_x1 = max(0, min(camera.horRes-1, (int)round(v1.x)));
	int pixel_y1 = max(0, min(camera.verRes-1, (int)round(v1.y)));
	int pixel_x2 = max(0, min(camera.horRes-1, (int)round(v2.x)));
	int pixel_y2 = max(0, min(camera.verRes-1, (int)round(v2.y)));
	
	double l1 = (v0.x+0.5)*(v1.y+0.5) - (v1.x+0.5)*(v0.y+0.5);
	double l2 = (v1.x+0.5)*(v2.y+0.5) - (v2.x+0.5)*(v1.y+0.5);
	double l3 = (v2.x+0.5)*(v0.y+0.5) - (v0.x+0.5)*(v2.y+0.5);

	double alpha = (v0.x+0.5)* (v1.y-v2.y) + (v0.y+0.5)*(v2.x-v1.x) + l2;
	double beta = (v1.x+0.5)* (v2.y-v0.y) + (v1.y+0.5)*(v0.x-v2.x) + l3;
	double gamma = (v2.x+0.5)* (v0.y-v1.y) + (v2.y+0.5)*(v1.x-v0.x) + l1;

	drawSolidLine(colorsOfVertices[v0.colorId], colorsOfVertices[v1.colorId], 
	colorsOfVertices[v2.colorId],min(pixel_x0, min(pixel_x1, pixel_x2)), 
	min(pixel_y0, min(pixel_y1, pixel_y2)), max(pixel_x0, max(pixel_x1, pixel_x2)),
	max(pixel_y0, max(pixel_y1, pixel_y2)), v0.y-v1.y, v1.y-v2.y, v2.y-v0.y, v1.x-v0.x, v2.x-v1.x, v0.x-v2.x, l1, l2, l3, alpha, beta, gamma);
}

void Scene::rasterizeSolidTriangleOrtho(Mesh *mesh, Vec4* triangle, Camera &camera){
	Vec4 v0 = triangle[0];
	Vec4 v1 = triangle[1];
	Vec4 v2 = triangle[2];

	v0 = multiplyMatrixWithVec4(mesh->V, v0);
	v1 = multiplyMatrixWithVec4(mesh->V, v1);
	v2 = multiplyMatrixWithVec4(mesh->V, v2); 

	int pixel_x0 = max(0, min(camera.horRes-1, (int)round(v0.x)));
	int pixel_y0 = max(0, min(camera.verRes-1, (int)round(v0.y)));
	int pixel_x1 = max(0, min(camera.horRes-1, (int)round(v1.x)));
	int pixel_y1 = max(0, min(camera.verRes-1, (int)round(v1.y)));
	int pixel_x2 = max(0, min(camera.horRes-1, (int)round(v2.x)));
	int pixel_y2 = max(0, min(camera.verRes-1, (int)round(v2.y)));
	
	double l1 = (v0.x+0.5)*(v1.y+0.5) - (v1.x+0.5)*(v0.y+0.5);
	double l2 = (v1.x+0.5)*(v2.y+0.5) - (v2.x+0.5)*(v1.y+0.5);
	double l3 = (v2.x+0.5)*(v0.y+0.5) - (v0.x+0.5)*(v2.y+0.5);

	double alpha = (v0.x+0.5)* (v1.y-v2.y) + (v0.y+0.5)*(v2.x-v1.x) + l2;
	double beta = (v1.x+0.5)* (v2.y-v0.y) + (v1.y+0.5)*(v0.x-v2.x) + l3;
	double gamma = (v2.x+0.5)* (v0.y-v1.y) + (v2.y+0.5)*(v1.x-v0.x) + l1;

	drawSolidLine(colorsOfVertices[v0.colorId], colorsOfVertices[v1.colorId], 
	colorsOfVertices[v2.colorId],min(pixel_x0, min(pixel_x1, pixel_x2)), 
	min(pixel_y0, min(pixel_y1, pixel_y2)), max(pixel_x0, max(pixel_x1, pixel_x2)),
	max(pixel_y0, max(pixel_y1, pixel_y2)), v0.y-v1.y, v1.y-v2.y, v2.y-v0.y, v1.x-v0.x, v2.x-v1.x, v0.x-v2.x, l1, l2, l3, alpha, beta, gamma);
}

void Scene::drawSolidLine(Color* c0, Color *c1, Color *c2, int min_x, int min_y,int max_x,int max_y,double a, double b, double c, double v,double w, double z,double l1, double l2, double l3,double alpha_d, double beta_d, double gamma_d){
	for (int i = min_y; i <= max_y; i++) {
		for (int j = min_x; j <= max_x; j++) {
			double alpha = ((i+0.5)* (b) + (j+0.5)*(w) + l2)/alpha_d;
			double beta = ((i+0.5)* (c) + (j+0.5)*(z) + l3)/beta_d;
			double gamma = ((i+0.5)* (a) + (j+0.5)*(v) + l1)/gamma_d;
			if (alpha >= 0 && beta >= 0 && gamma >= 0) {
				Color color;
				color.r = alpha * c0->r + beta * c1->r + gamma * c2->r;
				color.g = alpha * c0->g + beta * c1->g + gamma * c2->g;
				color.b = alpha * c0->b + beta * c1->b + gamma * c2->b;
				image[i][j] = color;
			}
		}
	}
}

 
void Scene::performRasterization(Vec4* triangle, Mesh* mesh, Camera& camera) {
	if (camera.projectionType == 1){
		if(mesh->type == 0){
			rasterizeWireframeTriangle(mesh, triangle, camera); //change loop inside of it
		}
		else if(mesh->type == 1){
			rasterizeSolidTriangle(mesh, triangle, camera);
		}
	}
	else {
		if(mesh->type == 0){
			rasterizeWireframeTriangleOrtho(mesh, triangle, camera);
		}
		else if(mesh->type == 1){
			rasterizeSolidTriangleOrtho(mesh, triangle, camera);
		}
	}
}