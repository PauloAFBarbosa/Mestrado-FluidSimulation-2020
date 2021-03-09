#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>


#include <stdio.h>

using namespace std;


#endif

#define _USE_MATH_DEFINES
#include <math.h>


#define NOQUERY 0
#define QUERYV1 1

GLuint program;
GLuint index_buffer;
int *ret = new int [1000];

bool mortonCode = true;

float RESTITUTION=0.5;
float TIMESTEP=0.01;
float SURFACETENSION=0.0728;
float THRESHOLD = 7.065;
float H = 0.0457; // Com o H a 10 ja apanha algumas particulas vizinhas
float STIFF = 3.0;
float DECLIVE = 0.0;
float MASS = 0.02;
float RESTDENSITY = 998.29;
float VISCOSITY = 3.5;

float GRAVITYVALUE = -9.8;
float GRAVITY[3] = { 0, GRAVITYVALUE, 0 };
//Valores max e min da simulação - determina o voluma da sim
float XMIN = 0, XMAX = 200;
float YMIN = 0, YMAX = 200;
float ZMIN = 0, ZMAX = 200;

//usado para a geração das particulas
const float fluidVolume = 1000 * MASS / RESTDENSITY;
const float particleDiameter = powf(fluidVolume, 1.0f / 3.0f) / 10;
const float particleRadius = particleDiameter / 2;
int const hashSize = 1000; //o cubo tem 0.4 de lado, e o h tem 0.04.. isto da 8 divisões por cada lado, arredondo para 10, por isso fica 10*10*10
//------------

//usado para determinar quando começar a simular
bool start = false;
bool simulateFrame = false; //usado para simular apenas um frame
int framesSimulated = 0;

//Particulas por lado do cubo, ou seja, particulas total é 7^3
int pontos_lado = 15;
int countPoints = 0;
//usado para inserir multiplos batches de pontos
int pontos_inseridos = 0;


float size = 2.5;
//numero maximo de pontos em cada celula da octree
int max_points = 100;
int max_depth = 5;

int frame = 0;
float fps = 0;
float center[3] = { 0,0,0 };
float size_query = H * 2;
//std::vector<Point*> points;

int mytime;
int timebase;
int maxTime = 10000;
float stepT = 100.0f / maxTime;

float alfa = 0.0f, beta = 0.5f, radius = 5.0f;
float camX = (XMIN+XMAX)/2, camY = 0, camZ = (ZMIN + ZMAX) / 2;
float rot = 0;
int winid = 0;




void processKeys(unsigned char c, int xx, int yy) {
	switch (c)
	{
	default:
		break;
	}
	glutPostRedisplay();
}

void processSpecialKeys(int key, int xx, int yy) {
	switch (key) {

	case GLUT_KEY_RIGHT:
		alfa -= 0.1; break;

	case GLUT_KEY_LEFT:
		alfa += 0.1; break;

	case GLUT_KEY_UP:
		beta += 0.1f;
		if (beta > 1.5f)
			beta = 1.5f;
		break;

	case GLUT_KEY_DOWN:
		beta -= 0.1f;
		if (beta < -1.5f)
			beta = -1.5f;
		break;

	case GLUT_KEY_PAGE_DOWN: radius -= 1.0f;
		if (radius < 1.0f)
			radius = 1.0f;
		break;

	case GLUT_KEY_PAGE_UP: radius += 1.0f; break;
	}
	
	glutPostRedisplay();
}

void changeSize(int w, int h) {
	// Prevent a divide by zero, when window is too short
	// (you cant make a window with zero width).
	if (h == 0)
		h = 1;

	// compute window's aspect ratio 
	float ratio = w * 1.0 / h;

	// Set the projection matrix as current
	glMatrixMode(GL_PROJECTION);
	// Load Identity Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set perspective
	gluPerspective(45.0f, ratio, 1.0f, 1000.0f);
	
	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}



void renderScene(void) {
	char fpss[200];

	static float t = 0;

	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set the camera
	glLoadIdentity();
	gluLookAt(camX, camY, camZ,
		(XMIN+XMAX)/2, (YMIN + YMAX) / 2, (ZMIN + ZMAX) / 2,
		0.0f, 1.0f, 0.0f);
	
	
	glUseProgram(program); // Compute shader program.
	//buffer 
	glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, index_buffer);
	//--------
	glDispatchCompute(1, 1, 1);
	
	glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, 1000 * 32, ret);
	
	//Fps

	frame++;
	mytime = glutGet(GLUT_ELAPSED_TIME);
	if (mytime - timebase > 1000) {
		fps = frame * 1000.0 / (mytime - timebase);
		timebase = mytime;
		frame = 0;
	}
	sprintf(fpss, "%f", fps);
	glutSetWindowTitle(fpss);
	

	glEnable(GL_LIGHTING);

	// End of frame
	glutSwapBuffers();
}



void setLighting(void) {
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	float lightPos[] = { 0.0f, 0.4f, 1.0f, 0.0f };
	float lightAmb[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float lightDif[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	float lightSpc[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDif);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpc);

	float matAmb[] = { 0.7f, 0.7f, 0.9f, 1.0f };
	float matDif[] = { 0.7f, 0.7f, 0.9f, 1.0f };
	float matSpc[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	float matShi[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_AMBIENT, matAmb);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, matDif);
	glMaterialfv(GL_FRONT, GL_SPECULAR, matSpc);
	glMaterialfv(GL_FRONT, GL_SHININESS, matShi);
}



int main(int argc, char** argv) {
	// put GLUT init here
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(799, 599);
	winid = glutCreateWindow("Octree");
	// put callback registration here

	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);

	glutKeyboardFunc(processKeys);
	glutSpecialFunc(processSpecialKeys);
	

	// OpenGL settings 
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	setLighting();
	// enter GLUT's main loop
	
	//compute shader
	const char* csSrc[] = {
		"#version 430\n",
		"layout (local_size_x = 1, local_size_y = 1) in;\
		layout(std430, binding = 0) buffer Index\
		{\
			int index_buffer[];\
		};\
         void main() {\
			for(int i =0; i < 1000 ; i ++){\
				index_buffer[i]=i;\
		}\
         }"
	};
	GLint GlewInitResult = glewInit();
	if (GLEW_OK != GlewInitResult)
	{
		printf("ERROR: %s ",glewGetErrorString(GlewInitResult));
		//exit(EXIT_FAILURE);
	}

	program = glCreateProgram();
	GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(shader,2, csSrc,NULL);
	glCompileShader(shader);

	//check errors
	GLint isCompiled = 0;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &isCompiled);
	if (isCompiled == GL_FALSE)
	{
		fprintf(stderr, "Error in compiling the compute shader\n");
		GLchar log[10240];
		GLsizei length;
		glGetShaderInfoLog(shader, 10239, &length, log);
		fprintf(stderr, "Compiler log:\n%s\n", log);
		//exit(40);
		
	}

	glAttachShader(program, shader);
	glLinkProgram(program);
	int rvalue;
	glGetProgramiv(program, GL_LINK_STATUS, &rvalue);
	if (!rvalue) {
		fprintf(stderr, "Error in linking compute shader program\n");
		GLchar log[10240];
		GLsizei length;
		glGetProgramInfoLog(program, 10239, &length, log);
		fprintf(stderr, "Linker log:\n%s\n", log);
		
	}

	
	

	//buffers
	
	glGenBuffers(1, &index_buffer);
	printf("buffer %d\n", index_buffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, index_buffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, 32*1000,NULL, GL_DYNAMIC_COPY);
	glutMainLoop();

	return 1;
}