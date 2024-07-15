#ifndef _SCENE_H_
#define _SCENE_H_
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include "Vec3.h"
#include "Vec4.h"
#include "Matrix4.h"
#include "Color.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Camera.h"
#include "Mesh.h"

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	std::vector<std::vector<Color>> image;
	std::vector<std::vector<double>> depth;
	std::vector<Camera *> cameras;
	std::vector<Vec3 *> vertices;
	std::vector<Color *> colorsOfVertices;
	std::vector<Scaling *> scalings;
	std::vector<Rotation *> rotations;
	std::vector<Translation *> translations;
	std::vector<Mesh *> meshes;

	Scene(const char *xmlPath);

	void assignColorToPixel(int i, int j, Color c);
	void initializeImage(Camera *camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera *camera);
	void convertPPMToPNG(std::string ppmFileName);
	void forwardRenderingPipeline(Camera &camera);
	
	Matrix4 applyModellingTransformations(Mesh *mesh);
	void applyTransformations(Mesh* mesh, Camera& camera);

	Matrix4 applyTranslations(Matrix4* transform, Translation* translation);
	Matrix4 applyRotations(Matrix4* transform, Rotation* rotation);
	Matrix4 applyScalings(Matrix4* transform, Scaling* scaling);
	Matrix4 rotationMatrix(Vec3 u, double theta);
	Matrix4 inverseTransposeMatrix(Matrix4& M);
	Matrix4 multiplyMatrixWithMatrix(Matrix4& M1, Matrix4& M2);
	void calculateMeshNormals(Mesh *mesh);
	void performBackfaceCulling(Mesh *mesh, Camera &camera);
	Matrix4 performCameraTransformation(Camera& camera);
	Matrix4 performProjectionTransformation(Camera &camera);
	Matrix4 performViewportTransformation(Mesh *mesh, Camera &camera);
	void performRasterization(Vec4* triangle, Mesh *mesh, Camera &camera);
	void rasterizeWireframeTriangle(Mesh *mesh, Vec4* triangle, Camera &camera);
	void rasterizeWireframeTriangleOrtho(Mesh *mesh, Vec4* triangle, Camera &camera);

	bool clipLinePers(Vec4& v1, Vec4& v2, Color& c1, Color& c2);
	bool clipLineOr(Vec4& v1, Vec4& v2, Color& c1, Color& c2);
	bool isInsideFrustum(const Vec4& v1, const Vec4& v2);
	bool calculateIntersection(const Vec4& v1, const Vec4& v2, Vec4& intersection);
	Color interpolateColor(const Vec4& v1, const Vec4& v2, const Color& c1, const Color& c2);
	bool isOutsideSamePlane(const Vec4& v1, const Vec4& v2);
	void drawLine(int x0, int y0, int x1, int y1, Color* c_0, Color* c_1);
	bool vis(double num, double denom, double &te, double &to);
	void rasterizeSolidTriangle(Mesh *mesh, Vec4* triangle, Camera &camera);
	void rasterizeSolidTriangleOrtho(Mesh *mesh, Vec4* triangle, Camera &camera);
	void drawSolidLine(Color* c0, Color *c1, Color *c2, int min_x, int min_y,int max_x,int max_y,double a, double b, double c, double v,double w, double z,double l1, double l2, double l3,double alpha, double beta, double gamma);
};

#endif
