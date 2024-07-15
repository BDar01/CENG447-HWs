#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <GL/glew.h>   // The GL Header File
#include <GL/gl.h>   // The GL Header File
#include <GLFW/glfw3.h> // The GLFW header
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <ft2build.h>
#include FT_FREETYPE_H

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

using namespace std;
static float gBunnyRotAngleY = 0.0f; // Global Y rotation angle ıf gBunnyFlag is true
static float gBunnyRotAngleX = 0.0f; // Global X rotation angle ıf gReturnFlag is true
float gBunyPosX = 0.0f; // Global position in X
float gBunyPosY = -5.0f; // Global position in Y
bool gBunnyUpFlag = false;
bool gBunnyFlag = false;
bool gReturnFlag = false;
bool gRestartFlag = false;
bool gCollisionFlag = false;
float gCubeObjectPosZ = 0.0f; // Global cube position in Z
float gGroundObjectPosZ = -80.0f; // Global ground position in Z
float moveSpeed = 0.5f;   // Movement speed
static int score = 0;
float gFixObjectPosZ = -100.0f; // Fixed position in Z

GLuint gProgram[5];
GLint gIntensityLoc;
float gIntensity = 1000;
int gWidth = 640, gHeight = 480;
int randomNum = rand() % 3;

struct Vertex
{
    Vertex(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Texture
{
    Texture(GLfloat inU, GLfloat inV) : u(inU), v(inV) { }
    GLfloat u, v;
};

struct Normal
{
    Normal(GLfloat inX, GLfloat inY, GLfloat inZ) : x(inX), y(inY), z(inZ) { }
    GLfloat x, y, z;
};

struct Face
{
	Face(int v[], int t[], int n[]) {
		vIndex[0] = v[0];
		vIndex[1] = v[1];
		vIndex[2] = v[2];
		tIndex[0] = t[0];
		tIndex[1] = t[1];
		tIndex[2] = t[2];
		nIndex[0] = n[0];
		nIndex[1] = n[1];
		nIndex[2] = n[2];
	}
    GLuint vIndex[3], tIndex[3], nIndex[3];
};

vector<Vertex> gVertices;
vector<Texture> gTextures;
vector<Normal> gNormals;
vector<Face> gFaces;

vector<Vertex> gQuadVertices;
vector<Texture> gQuadTextures;
vector<Normal> gQuadNormals;
vector<Face> gQuadFaces;

vector<Vertex> gCubeVertices;
vector<Texture> gCubeTextures;
vector<Normal> gCubeNormals;
vector<Face> gCubeFaces;
//for bunny
GLuint gVertexAttribBuffer, gTextVBO, gIndexBuffer;
GLint gInVertexLoc, gInNormalLoc;
int gVertexDataSizeInBytes, gNormalDataSizeInBytes;
//for quad
GLuint quadVertexAttribBuffer, quadIndexBuffer;
GLint quadInVertexLoc, quadInNormalLoc;
int quadVertexDataSizeInBytes, quadNormalDataSizeInBytes;
//for cube
GLuint cubeVertexAttribBuffer, cubeIndexBuffer;
GLint cubeInVertexLoc, cubeInNormalLoc;
int cubeVertexDataSizeInBytes, cubeNormalDataSizeInBytes;


/// Holds all state information relevant to a character as loaded using FreeType
struct Character {
    GLuint TextureID;   // ID handle of the glyph texture
    glm::ivec2 Size;    // Size of glyph
    glm::ivec2 Bearing;  // Offset from baseline to left/top of glyph
    GLuint Advance;    // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;


bool ParseObj(const string& fileName)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == '#') // comment
                {
                    continue;
                }
                else if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        gTextures.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
					char c;
					int vIndex[3],  nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0]; 
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1]; 
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2]; 

					assert(vIndex[0] == nIndex[0] &&
						   vIndex[1] == nIndex[1] &&
						   vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

                    gFaces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            //data += curLine;
            if (!myfile.eof())
            {
                //data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

	/*
	for (int i = 0; i < gVertices.size(); ++i)
	{
		Vector3 n;

		for (int j = 0; j < gFaces.size(); ++j)
		{
			for (int k = 0; k < 3; ++k)
			{
				if (gFaces[j].vIndex[k] == i)
				{
					// face j contains vertex i
					Vector3 a(gVertices[gFaces[j].vIndex[0]].x, 
							  gVertices[gFaces[j].vIndex[0]].y,
							  gVertices[gFaces[j].vIndex[0]].z);

					Vector3 b(gVertices[gFaces[j].vIndex[1]].x, 
							  gVertices[gFaces[j].vIndex[1]].y,
							  gVertices[gFaces[j].vIndex[1]].z);

					Vector3 c(gVertices[gFaces[j].vIndex[2]].x, 
							  gVertices[gFaces[j].vIndex[2]].y,
							  gVertices[gFaces[j].vIndex[2]].z);

					Vector3 ab = b - a;
					Vector3 ac = c - a;
					Vector3 normalFromThisFace = (ab.cross(ac)).getNormalized();
					n += normalFromThisFace;
				}

			}
		}

		n.normalize();

		gNormals.push_back(Normal(n.x, n.y, n.z));
	}
	*/

	assert(gVertices.size() == gNormals.size());

    return true;
}
bool ParseQuadObj(const string& fileName)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == '#') // comment
                {
                    continue;
                }
                else if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        gQuadTextures.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gQuadNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gQuadVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
					char c;
					int vIndex[3],  nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0]; 
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1]; 
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2]; 

					assert(vIndex[0] == nIndex[0] &&
						   vIndex[1] == nIndex[1] &&
						   vIndex[2] == nIndex[2]); // a limitation for now
					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

                    gQuadFaces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            //data += curLine;
            if (!myfile.eof())
            {
                //data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }


	assert(gVertices.size() == gNormals.size());

    return true;
}
bool ParseCubeObj(const string& fileName)
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine; 

        while (getline(myfile, curLine))
        {
            stringstream str(curLine);
            GLfloat c1, c2, c3;
            GLuint index[9];
            string tmp;

            if (curLine.length() >= 2)
            {
                if (curLine[0] == '#') // comment
                {
                    continue;
                }
                else if (curLine[0] == 'v')
                {
                    if (curLine[1] == 't') // texture
                    {
                        str >> tmp; // consume "vt"
                        str >> c1 >> c2;
                        gCubeTextures.push_back(Texture(c1, c2));
                    }
                    else if (curLine[1] == 'n') // normal
                    {
                        str >> tmp; // consume "vn"
                        str >> c1 >> c2 >> c3;
                        gCubeNormals.push_back(Normal(c1, c2, c3));
                    }
                    else // vertex
                    {
                        str >> tmp; // consume "v"
                        str >> c1 >> c2 >> c3;
                        gCubeVertices.push_back(Vertex(c1, c2, c3));
                    }
                }
                else if (curLine[0] == 'f') // face
                {
                    str >> tmp; // consume "f"
					char c;
					int vIndex[3],  nIndex[3], tIndex[3];
					str >> vIndex[0]; str >> c >> c; // consume "//"
					str >> nIndex[0]; 
					str >> vIndex[1]; str >> c >> c; // consume "//"
					str >> nIndex[1]; 
					str >> vIndex[2]; str >> c >> c; // consume "//"
					str >> nIndex[2]; 

					assert(vIndex[0] == nIndex[0] &&
						   vIndex[1] == nIndex[1] &&
						   vIndex[2] == nIndex[2]); // a limitation for now

					// make indices start from 0
					for (int c = 0; c < 3; ++c)
					{
						vIndex[c] -= 1;
						nIndex[c] -= 1;
						tIndex[c] -= 1;
					}

                    gCubeFaces.push_back(Face(vIndex, tIndex, nIndex));
                }
                else
                {
                    cout << "Ignoring unidentified line in obj file: " << curLine << endl;
                }
            }

            //data += curLine;
            if (!myfile.eof())
            {
                //data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

	assert(gVertices.size() == gNormals.size());

    return true;
}
bool ReadDataFromFile(
    const string& fileName, ///< [in]  Name of the shader file
    string&       data)     ///< [out] The contents of the file
{
    fstream myfile;

    // Open the input 
    myfile.open(fileName.c_str(), std::ios::in);

    if (myfile.is_open())
    {
        string curLine;

        while (getline(myfile, curLine))
        {
            data += curLine;
            if (!myfile.eof())
            {
                data += "\n";
            }
        }

        myfile.close();
    }
    else
    {
        return false;
    }

    return true;
}

void createVS(GLuint& program, const string& filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &shader, &length);
    glCompileShader(vs);

    char output[1024] = {0};
    glGetShaderInfoLog(vs, 1024, &length, output);
    printf("VS compile log: %s\n", output);

    glAttachShader(program, vs);
}

void createFS(GLuint& program, const string& filename)
{
    string shaderSource;

    if (!ReadDataFromFile(filename, shaderSource))
    {
        cout << "Cannot find file name: " + filename << endl;
        exit(-1);
    }

    GLint length = shaderSource.length();
    const GLchar* shader = (const GLchar*) shaderSource.c_str();

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &shader, &length);
    glCompileShader(fs);

    char output[1024] = {0};
    glGetShaderInfoLog(fs, 1024, &length, output);
    printf("FS compile log: %s\n", output);

    glAttachShader(program, fs);
}

void initShaders()
{
    gProgram[0] = glCreateProgram();
    gProgram[1] = glCreateProgram();
    gProgram[2] = glCreateProgram();
    gProgram[3] = glCreateProgram();
    gProgram[4] = glCreateProgram();

    createVS(gProgram[0], "vert0.glsl");
    createFS(gProgram[0], "frag0.glsl");

    createVS(gProgram[1], "vert1.glsl");
    createFS(gProgram[1], "frag1.glsl");
    
    createVS(gProgram[2], "vert_text.glsl");
    createFS(gProgram[2], "frag_text.glsl");

    createVS(gProgram[3], "vert2.glsl");
    createFS(gProgram[3], "frag2.glsl");

    createVS(gProgram[4], "vert2.glsl");
    createFS(gProgram[4], "frag3.glsl");

    glBindAttribLocation(gProgram[0], 0, "inVertex");
    glBindAttribLocation(gProgram[0], 1, "inNormal");

    glBindAttribLocation(gProgram[1], 0, "inVertex");
    glBindAttribLocation(gProgram[1], 1, "inNormal");

    glBindAttribLocation(gProgram[3], 0, "inVertex");
    glBindAttribLocation(gProgram[3], 1, "inNormal");

    glBindAttribLocation(gProgram[4], 0, "inVertex");
    glBindAttribLocation(gProgram[4], 1, "inNormal");

    glBindAttribLocation(gProgram[2], 2, "vertex");

    glLinkProgram(gProgram[0]);
    glLinkProgram(gProgram[1]);
    glLinkProgram(gProgram[2]);
    glLinkProgram(gProgram[3]);
    glLinkProgram(gProgram[4]);
    glUseProgram(gProgram[0]);

    gIntensityLoc = glGetUniformLocation(gProgram[0], "intensity");
    cout << "gIntensityLoc = " << gIntensityLoc << endl;
    glUniform1f(gIntensityLoc, gIntensity);
}

void initVBO()
{
    glEnableVertexAttribArray(0); // Vertex positions
    glEnableVertexAttribArray(1); // Vertex normals
    assert(glGetError() == GL_NONE);

    glGenBuffers(1, &gVertexAttribBuffer);
    glGenBuffers(1, &gIndexBuffer);

    assert(gVertexAttribBuffer > 0 && gIndexBuffer > 0);

    glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

    // Calculate the sizes of the vertex and normal data
    gVertexDataSizeInBytes = gVertices.size() * sizeof(Vertex);
    gNormalDataSizeInBytes = gNormals.size() * sizeof(Normal);
    int indexDataSizeInBytes = gFaces.size() * 3 * sizeof(GLuint);

    // Create arrays to hold the vertex and normal data
    GLfloat* vertexData = new GLfloat[gVertices.size() * 3];
    GLfloat* normalData = new GLfloat[gNormals.size() * 3];
    GLuint* indexData = new GLuint[gFaces.size() * 3];

    // Copy vertex and normal data into the arrays
    for (size_t i = 0; i < gVertices.size(); ++i) {
        vertexData[3 * i] = gVertices[i].x;
        vertexData[3 * i + 1] = gVertices[i].y;
        vertexData[3 * i + 2] = gVertices[i].z;
    }

    for (size_t i = 0; i < gNormals.size(); ++i) {
        normalData[3 * i] = gNormals[i].x;
        normalData[3 * i + 1] = gNormals[i].y;
        normalData[3 * i + 2] = gNormals[i].z;
    }

    // Copy index data into the array
    for (size_t i = 0; i < gFaces.size(); ++i) {
        indexData[3 * i] = gFaces[i].vIndex[0];
        indexData[3 * i + 1] = gFaces[i].vIndex[1];
        indexData[3 * i + 2] = gFaces[i].vIndex[2];
    }

    // Upload the vertex and normal data to the GPU
    glBufferData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes + gNormalDataSizeInBytes, NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, gVertexDataSizeInBytes, vertexData);
    glBufferSubData(GL_ARRAY_BUFFER, gVertexDataSizeInBytes, gNormalDataSizeInBytes, normalData);

    // Upload the index data to the GPU
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexDataSizeInBytes, indexData, GL_STATIC_DRAW);

    // Set up the vertex attribute pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

    // Clean up
    delete[] vertexData;
    delete[] normalData;
    delete[] indexData;
}
void initQuadVBO()
{
    glEnableVertexAttribArray(0); // Vertex positions
    glEnableVertexAttribArray(1); // Vertex normals
    assert(glGetError() == GL_NONE);

    glGenBuffers(1, &quadVertexAttribBuffer);
    glGenBuffers(1, &quadIndexBuffer);

    assert(quadVertexAttribBuffer > 0 && quadIndexBuffer > 0);

    glBindBuffer(GL_ARRAY_BUFFER, quadVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, quadIndexBuffer);

    // Calculate the sizes of the vertex and normal data
    quadVertexDataSizeInBytes = gQuadVertices.size() * sizeof(Vertex);
    quadNormalDataSizeInBytes = gQuadNormals.size() * sizeof(Normal);
    int quadIndexDataSizeInBytes = gQuadFaces.size() * 3 * sizeof(GLuint);

    // Create arrays to hold the vertex and normal data
    GLfloat* vertexData = new GLfloat[gQuadVertices.size() * 3];
    GLfloat* normalData = new GLfloat[gQuadNormals.size() * 3];
    GLuint* indexData = new GLuint[gQuadFaces.size() * 3];

    // Copy vertex and normal data into the arrays
    for (size_t i = 0; i < gQuadVertices.size(); ++i) {
        vertexData[3 * i] = gQuadVertices[i].x;
        vertexData[3 * i + 1] = gQuadVertices[i].y;
        vertexData[3 * i + 2] = gQuadVertices[i].z;
    }

    for (size_t i = 0; i < gQuadNormals.size(); ++i) {
        normalData[3 * i] = gQuadNormals[i].x;
        normalData[3 * i + 1] = gQuadNormals[i].y;
        normalData[3 * i + 2] = gQuadNormals[i].z;
    }

    // Copy index data into the array
    for (size_t i = 0; i < gQuadFaces.size(); ++i) {
        indexData[3 * i] = gQuadFaces[i].vIndex[0];
        indexData[3 * i + 1] = gQuadFaces[i].vIndex[1];
        indexData[3 * i + 2] = gQuadFaces[i].vIndex[2];
    }

    // Upload the vertex and normal data to the GPU
    glBufferData(GL_ARRAY_BUFFER, quadVertexDataSizeInBytes + quadNormalDataSizeInBytes, NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, quadVertexDataSizeInBytes, vertexData);
    glBufferSubData(GL_ARRAY_BUFFER, quadVertexDataSizeInBytes, quadNormalDataSizeInBytes, normalData);

    // Upload the index data to the GPU
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, quadIndexDataSizeInBytes, indexData, GL_STATIC_DRAW);

    // Set up the vertex attribute pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(quadVertexDataSizeInBytes));

    // Clean up
    delete[] vertexData;
    delete[] normalData;
    delete[] indexData;
}

void initCubeVBO()
{
    glEnableVertexAttribArray(0); // Vertex positions
    glEnableVertexAttribArray(1); // Vertex normals
    assert(glGetError() == GL_NONE);

    glGenBuffers(1, &cubeVertexAttribBuffer);
    glGenBuffers(1, &cubeIndexBuffer);

    assert(cubeVertexAttribBuffer > 0 && cubeIndexBuffer > 0);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeIndexBuffer);

    // Calculate the sizes of the vertex and normal data
    cubeVertexDataSizeInBytes = gCubeVertices.size() * sizeof(Vertex);
    cubeNormalDataSizeInBytes = gCubeNormals.size() * sizeof(Normal);
    int cubeIndexDataSizeInBytes = gCubeFaces.size() * 3 * sizeof(GLuint);

    // Create arrays to hold the vertex and normal data
    GLfloat* vertexData = new GLfloat[gCubeVertices.size() * 3];
    GLfloat* normalData = new GLfloat[gCubeNormals.size() * 3];
    GLuint* indexData = new GLuint[gCubeFaces.size() * 3];

    // Copy vertex and normal data into the arrays
    for (size_t i = 0; i < gCubeVertices.size(); ++i) {
        vertexData[3 * i] = gCubeVertices[i].x;
        vertexData[3 * i + 1] = gCubeVertices[i].y;
        vertexData[3 * i + 2] = gCubeVertices[i].z;
    }

    for (size_t i = 0; i < gCubeNormals.size(); ++i) {
        normalData[3 * i] = gCubeNormals[i].x;
        normalData[3 * i + 1] = gCubeNormals[i].y;
        normalData[3 * i + 2] = gCubeNormals[i].z;
    }

    // Copy index data into the array
    for (size_t i = 0; i < gCubeFaces.size(); ++i) {
        indexData[3 * i] = gCubeFaces[i].vIndex[0];
        indexData[3 * i + 1] = gCubeFaces[i].vIndex[1];
        indexData[3 * i + 2] = gCubeFaces[i].vIndex[2];
    }

    // Upload the vertex and normal data to the GPU
    glBufferData(GL_ARRAY_BUFFER, cubeVertexDataSizeInBytes + cubeNormalDataSizeInBytes, NULL, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, cubeVertexDataSizeInBytes, vertexData);
    glBufferSubData(GL_ARRAY_BUFFER, cubeVertexDataSizeInBytes, cubeNormalDataSizeInBytes, normalData);

    // Upload the index data to the GPU GLFW_REPEA
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, cubeIndexDataSizeInBytes, indexData, GL_STATIC_DRAW);

    // Set up the vertex attribute pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(cubeVertexDataSizeInBytes));

    // Clean up
    delete[] vertexData;
    delete[] normalData;
    delete[] indexData;
}

void initFonts(int windowWidth, int windowHeight)
{
    // Set OpenGL options
    //glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(windowWidth), 0.0f, static_cast<GLfloat>(windowHeight));
    glUseProgram(gProgram[2]);
    glUniformMatrix4fv(glGetUniformLocation(gProgram[2], "projection"), 1, GL_FALSE, glm::value_ptr(projection));

    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
    }

    // Load font as face
    FT_Face face;
    if (FT_New_Face(ft, "C:\\Windows\\Fonts\\Arial.ttf", 0, &face)) //Change text for Linux/Windows
    {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
    }

    // Set size to load glyphs as
    FT_Set_Pixel_Sizes(face, 0, 48);

    // Disable byte-alignment restriction
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); 

    // Load first 128 characters of ASCII set
    for (GLubyte c = 0; c < 128; c++)
    {
        // Load character glyph 
        if (FT_Load_Char(face, c, FT_LOAD_RENDER))
        {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
            continue;
        }
        // Generate texture
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                face->glyph->bitmap.width,
                face->glyph->bitmap.rows,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                face->glyph->bitmap.buffer
                );
        // Set texture options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // Now store character for later use
        Character character = {
            texture,
            glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
            glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
            face->glyph->advance.x
        };
        Characters.insert(std::pair<GLchar, Character>(c, character));
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    // Destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);

    //
    // Configure VBO for texture quads
    //
    glGenBuffers(1, &gTextVBO);
    glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 6 * 4, NULL, GL_DYNAMIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), 0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void init() 
{
	//ParseObj("armadillo.obj");
	ParseObj("bunny.obj");
    ParseQuadObj("quad.obj");
    ParseCubeObj("cube.obj");

    glEnable(GL_DEPTH_TEST);
    initShaders();
    initFonts(gWidth, gHeight);
    initVBO();
    initQuadVBO();
    initCubeVBO();
}

void drawModel()
{
	glBindBuffer(GL_ARRAY_BUFFER, gVertexAttribBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gIndexBuffer);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(gVertexDataSizeInBytes));

	glDrawElements(GL_TRIANGLES, gFaces.size() * 3, GL_UNSIGNED_INT, 0);
    // Unbind VBO and EBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void drawQuad()
{
    glBindBuffer(GL_ARRAY_BUFFER, quadVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, quadIndexBuffer);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(quadVertexDataSizeInBytes));

    glDrawElements(GL_TRIANGLES, gQuadFaces.size() * 3, GL_UNSIGNED_INT, 0);
    // Unbind VBO and EBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void drawCube()
{
    glBindBuffer(GL_ARRAY_BUFFER, cubeVertexAttribBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeIndexBuffer); 

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(cubeVertexDataSizeInBytes));

    glDrawElements(GL_TRIANGLES, gCubeFaces.size() * 3, GL_UNSIGNED_INT, 0);
    // Unbind VBO and EBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void renderText(const std::string& text, GLfloat x, GLfloat y, GLfloat scale, glm::vec3 color)
{
    // Activate corresponding render state	
    glUseProgram(gProgram[2]);
    glUniform3f(glGetUniformLocation(gProgram[2], "textColor"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);

    // Iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) 
    {
        Character ch = Characters[*c];

        GLfloat xpos = x + ch.Bearing.x * scale;
        GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

        GLfloat w = ch.Size.x * scale;
        GLfloat h = ch.Size.y * scale;

        // Update VBO for each character
        GLfloat vertices[6][4] = {
            { xpos,     ypos + h,   0.0, 0.0 },            
            { xpos,     ypos,       0.0, 1.0 },
            { xpos + w, ypos,       1.0, 1.0 },

            { xpos,     ypos + h,   0.0, 0.0 },
            { xpos + w, ypos,       1.0, 1.0 },
            { xpos + w, ypos + h,   1.0, 0.0 }           
        };

        // Render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);

        // Update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, gTextVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

        //glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)

        x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 translate= 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }

    glBindTexture(GL_TEXTURE_2D, 0);
}


void display()
{
    //random number generator
    glClearColor(0, 0, 0, 1);
    glClearDepth(1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    if (gReturnFlag) {
        gBunyPosY = -5.0f;
    }
    else {
        if (gBunyPosY <= -5.0f) {
            gBunnyUpFlag = true;
        }
        else if (gBunyPosY >= -3.0f) {
            gBunnyUpFlag = false;
        }
        if (gBunnyUpFlag) {
            gBunyPosY += (moveSpeed / 2.f);
        }
        else {
            gBunyPosY -= (moveSpeed / 2.f);
        }
    }

    if(gCubeObjectPosZ >= 0.0f)
    {
        randomNum = rand() % 3;
        printf("randomNum = %d\n", randomNum);
        gCubeObjectPosZ = gFixObjectPosZ;
        gGroundObjectPosZ = -80.0f;
        gCollisionFlag = false;
    }

    // Collision - Bunny and Cube
    //printf("gCubeObjectPosZ = %f\n", gCubeObjectPosZ );
    if (!gCollisionFlag && gCubeObjectPosZ > -6.3f){
        if (gBunyPosX <= -5.0f){
            gCollisionFlag = true;
            printf("hitted left most cube\n");
            if (randomNum == 0) {                
                gBunnyFlag = true;
                score += 1000;
                gBunnyRotAngleY = 0.0f;
                printf("hitted yellow\n");
            }
            else {
                gBunnyFlag = false;
                gReturnFlag = true;
                printf("hitted red\n");
            }
        }
        if (gBunyPosX >= -1.2f && gBunyPosX <= 1.2f){
            gCollisionFlag = true;
            printf("hitted middle cube\n");
            if (randomNum == 1) {       
                gBunnyFlag = true;
                score += 1000;
                gBunnyRotAngleY = 0.0f;
                printf("hitted yellow\n");
            }
            else {
                gBunnyFlag = false;
                gReturnFlag = true;
                printf("hitted red\n");
            }
        }
        if (gBunyPosX >= 5.f){
            gCollisionFlag = true;
            score += 1000;
            printf("hitted right most cube\n");
            if (randomNum == 2) {
                gBunnyFlag = true;
                gBunnyRotAngleY = 0.0f;
                printf("hitted yellow\n");
            }
            else {
                gBunnyFlag = false;
                gReturnFlag = true;
                printf("hitted red\n");
            }
        }
    }

    static float bunnyAngle = -100;

    if (gBunnyFlag && gBunnyRotAngleY < 360.0f){
        gBunnyRotAngleY += 10*moveSpeed;
    }

    if (gReturnFlag && gBunnyRotAngleX > -70.0f) {
        gBunnyRotAngleX -= 3.f;
    }

    // Object's transformation matrix
    glm::mat4 T;
    glm::mat4 R;
    if (gReturnFlag) {
        glm::mat4 R1 = glm::rotate(glm::mat4(1.f), glm::radians(bunnyAngle), glm::vec3(0, 1, 0));
        glm::mat4 R2 = glm::rotate(glm::mat4(1.f), glm::radians(gBunnyRotAngleX), glm::vec3(0, 0, 1));
        T = glm::translate(glm::mat4(1.f), glm::vec3(gBunyPosX, gBunyPosY-(0.01*gBunnyRotAngleX), -6.f));

        R = R2 * R1;
    }
    else {
        R = glm::rotate(glm::mat4(1.f), glm::radians(bunnyAngle + gBunnyRotAngleY), glm::vec3(0, 1, 0));
        T = glm::translate(glm::mat4(1.f), glm::vec3(gBunyPosX, gBunyPosY, -6.f));
    }
    glm::mat4 modelMat = T * R;

    glm::vec3 objectPos = glm::vec3(modelMat[3]);

    float cameraHeight = 10.0f; // Height of the camera above the object
    float cameraDistanceZ = 10.0f; // Fixed distance from the object in Z-direction

    // Calculate the camera position to maintain a fixed distance in Z-direction
    glm::vec3 cameraPos = glm::vec3(objectPos.x, objectPos.y + cameraHeight, objectPos.z + cameraDistanceZ);

    // Create the view matrix to look at the object
    glm::mat4 viewMat = glm::lookAt(
                        cameraPos,                       // Camera's position
                        objectPos,                       // Point to look at (object's position)
                        glm::vec3(0.0f, 1.0f, 0.0f));
    
    glm::mat4 perspMat = glm::perspective(glm::radians(90.0f), (float)gWidth / (float)gHeight, 0.1f, 200.0f);

    glUseProgram(gProgram[0]);
    glm::mat4 modelMatInv = glm::transpose(glm::inverse(modelMat));

    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[0], "viewMat"), 1, GL_FALSE, glm::value_ptr(viewMat));
    
    drawModel();
    glUseProgram(gProgram[1]);
    static float quadAngle = -90;
    T = glm::translate(glm::mat4(1.f), glm::vec3(0.f, -5.f, gGroundObjectPosZ));
    printf("Ground z: %f\n", gGroundObjectPosZ);
    glm::mat4 S = glm::scale(glm::mat4(1.0), glm::vec3(10.0, 1.0, 400.0));
    R = glm::rotate(glm::mat4(1.f), glm::radians(quadAngle), glm::vec3(1, 0, 0));
    modelMat = T * S * R;
    modelMatInv = glm::transpose(glm::inverse(modelMat));

    float currentOffset = gCubeObjectPosZ;
    if (gReturnFlag) {
        currentOffset = 0.0;
    }
    GLfloat currentOffValue = static_cast<GLfloat>(currentOffset);

    GLint timeUniformLocation = glGetUniformLocation(gProgram[1], "currentOffset");
    glUniform1f(timeUniformLocation, currentOffValue);
    glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[1], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    drawQuad();
    
    if(randomNum == 0){ 
        glUseProgram(gProgram[3]);
    }
    else {
        glUseProgram(gProgram[4]);
    }
    T = glm::translate(glm::mat4(1.f), glm::vec3(-7.f, -3.f, gCubeObjectPosZ));
    S = glm::scale(glm::mat4(1.0), glm::vec3(1.0, 2.0, 1.0));
    modelMat = T * S;
    modelMatInv = glm::transpose(glm::inverse(modelMat));

    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    if (!gCollisionFlag) {
        drawCube();
    }

    if(randomNum == 1){ 
        glUseProgram(gProgram[3]);
    }
    else {
        glUseProgram(gProgram[4]);
    }
    T = glm::translate(glm::mat4(1.f), glm::vec3(0.f, -3.f, gCubeObjectPosZ));
    S = glm::scale(glm::mat4(1.0), glm::vec3(1.0, 2.0, 1.0));
    modelMat = T * S;
    modelMatInv = glm::transpose(glm::inverse(modelMat));

    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    if (!gCollisionFlag) {
        drawCube();
    }

    if(randomNum == 2){ 
        glUseProgram(gProgram[3]);
    }
    else {
        glUseProgram(gProgram[4]);
    }
    T = glm::translate(glm::mat4(1.f), glm::vec3(7.f, -3.f, gCubeObjectPosZ));
    S = glm::scale(glm::mat4(1.0), glm::vec3(1.0, 2.0, 1.0));
    modelMat = T * S;
    modelMatInv = glm::transpose(glm::inverse(modelMat));

    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "modelingMat"), 1, GL_FALSE, glm::value_ptr(modelMat));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "modelingMatInvTr"), 1, GL_FALSE, glm::value_ptr(modelMatInv));
    glUniformMatrix4fv(glGetUniformLocation(gProgram[3], "perspectiveMat"), 1, GL_FALSE, glm::value_ptr(perspMat));

    if (!gCollisionFlag) {
        drawCube();
    }

    if (gReturnFlag) {
        renderText("Score: " + std::to_string(score), 0, 440, 1, glm::vec3(1, 0, 0));
        return;
    }
    else {
        renderText("Score: " + std::to_string(score), 0, 440, 1, glm::vec3(1, 1, 0));
        gCubeObjectPosZ += moveSpeed;
        gGroundObjectPosZ += moveSpeed;
        score += ceil(moveSpeed);
        moveSpeed += 0.0001f;
    }
}



void reshape(GLFWwindow* window, int w, int h)
{
    w = w < 1 ? 1 : w;
    h = h < 1 ? 1 : h;

    gWidth = w;
    gHeight = h;

    glViewport(0, 0, w, h);
}

void keyboard(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS || action == GLFW_REPEAT)
    {
        switch (key)
        {
            case GLFW_KEY_UP:
                //gObjectPosZ -= moveSpeed;
                break;
            case GLFW_KEY_DOWN:
                //gObjectPosZ += moveSpeed;
                break;
            case GLFW_KEY_LEFT:
                if(gReturnFlag || gBunyPosX <= -6.0f)
                    break; 
                gBunyPosX -= 1.5*moveSpeed;
                
                break;
            case GLFW_KEY_RIGHT:
                if(gReturnFlag || gBunyPosX >= 6.0f)
                    break;
                gBunyPosX += 1.5*moveSpeed;
                break;
            case GLFW_KEY_R:
                gRestartFlag = true;
        }
    }
}

void mainLoop(GLFWwindow* window)
{
    while (!glfwWindowShouldClose(window))
    {
        if (gRestartFlag) {
            gBunnyRotAngleY = 0.0f; // Global Y rotation angle ıf gBunnyFlag is true
            gBunnyRotAngleX = 0.0f; // Global X rotation angle ıf gReturnFlag is true
            gBunyPosX = 0.0f; // Global position in X
            gBunyPosY = -5.0f; // Global position in Y
            gBunnyUpFlag = false;
            gBunnyFlag = false;
            gReturnFlag = false;
            gRestartFlag = false;
            gCollisionFlag = false;
            gCubeObjectPosZ = 0.0f; // Global cube position in Z
            gGroundObjectPosZ = -20.0f; // Global ground position in Z
            moveSpeed = 0.5f;   // Movement speed
            score = 0;
            gFixObjectPosZ = -100.0f; // Fixed position in Z
            randomNum = rand() % 3;
            
            gRestartFlag = false;
        }
        display();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

int main(int argc, char** argv)   // Create Main Function For Bringing It All Together
{
    GLFWwindow* window;
    if (!glfwInit())
    {
        exit(-1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(gWidth, gHeight, "Simple Example", NULL, NULL);

    if (!window)
    {
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize GLEW to setup the OpenGL Function pointers
    if (GLEW_OK != glewInit())
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return EXIT_FAILURE;
    }

    char rendererInfo[512] = {0};
    strcpy(rendererInfo, (const char*) glGetString(GL_RENDERER));
    strcat(rendererInfo, " - ");
    strcat(rendererInfo, (const char*) glGetString(GL_VERSION));
    glfwSetWindowTitle(window, rendererInfo);

    init();

    glfwSetKeyCallback(window, keyboard);
    glfwSetWindowSizeCallback(window, reshape);

    reshape(window, gWidth, gHeight); // need to call this once ourselves
    mainLoop(window); // this does not return unless the window is closed

    glfwDestroyWindow(window);
    glfwTerminate();
       
    std::cin.get();
    return 0;
}

