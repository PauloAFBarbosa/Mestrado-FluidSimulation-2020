#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <iostream>
#include <chrono>
#include "Octree.h"
#include<GL/glut.h>
#include<vector>
#include "../glm/gtx/norm.hpp"
#include "../AntTweakBar/include/AntTweakBar.h"

#endif

#define _USE_MATH_DEFINES
#include <math.h>



double RESTITUTION=0.5;
double TIMESTEP=0.01f;
double SURFACETENSION=0.0728;
double THRESHOLD = 7.065;
double H = 0.0457; // Com o H a 10 ja apanha algumas particulas vizinhas
double STIFF = 3.0;
double DECLIVE = 0.0;
double MASS = 0.02;
double RESTDENSITY = 998.29;
double VISCOSITY = 3.5;

double GRAVITYVALUE = -9.8;
glm::vec3 GRAVITY = glm::vec3(0, GRAVITYVALUE, 0);
//Valores max e min da simulação - determina o voluma da sim
double XMIN = 0, XMAX = 200;
double YMIN = 0, YMAX = 200;
double ZMIN = 0, ZMAX = 200;

//usado para a geração das particulas
const float fluidVolume = 1000 * MASS / RESTDENSITY;
const float particleDiameter = powf(fluidVolume, 1.0f / 3.0f) / 10;
const float particleRadius = particleDiameter / 2;
//------------

//usado para determinar quando começar a simular
bool start = false;

//Diz quantos pontos sao inseridos (points_p*points_p*points_p)
int points_p = 9;
//usado para inserir multiplos batches de pontos
int pontos_inseridos = 0;

Octree* o;
double size = 2.5;
//numero maximo de pontos em cada celula da octree
int max_points = 100;
int max_depth = 5;

int frame = 0;
double fps = 0;
glm::vec3 center = glm::vec3(0, 0, 0);
double size_query = H * 2;
std::vector<Point*> points;

int mytime;
int timebase;
int maxTime = 10000;
double stepT = 100.0f / maxTime;

double alfa = 0.0f, beta = 0.5f, radius = 5.0f;
double camX = (XMIN+XMAX)/2, camY = 0, camZ = (ZMIN + ZMAX) / 2;
double rot = 0;
int winid = 0;

//Da maneira que fiz agora ele verifica 1734 pontos. Em vez de 10 000 pontos.



/**
 * @brief calculates the cam values (used in GluLookAt function) from the alteration of alfa, beta and radius
 *
 */
void spherical2Cartesian() {
	camX = radius * cos(beta) * sin(alfa) + (XMIN+XMAX)/2;
	camY = radius * sin(beta) + (YMIN + YMAX) / 2;
	camZ = radius * cos(beta) * cos(alfa) + (ZMIN + ZMAX) / 2;
}

void restart() {
	points.clear();
	pontos_inseridos = 0;
	int countPoints = 0;
	std::chrono::steady_clock::time_point begin, end;
	
	
	Point* p;
	
	
	for (float x = -particleRadius * points_p; x <= particleRadius * points_p; x += particleDiameter) {
		for (float y = -particleRadius * points_p; y <= particleRadius * points_p; y += particleDiameter) {
			for (float z = -particleRadius * points_p; z <= particleRadius * points_p; z += particleDiameter) {
				p = new Point(x, y, z);
				points.push_back(p);
				countPoints++;
			}
		}
	}

	//separei a parte de inserir os pontos porque quero medir apenas a parte de inserir.
	//no programa ele vai gerar os pontos 1 vez e vai inserir numa nova arvore cada frame
	begin = std::chrono::steady_clock::now();
	for (int i = 0; i < points_p; i++) {
		o->insertPoint(points.at(pontos_inseridos + i));
	}

	end = std::chrono::steady_clock::now();
	pontos_inseridos += points_p;
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[micro seconds]" << std::endl;
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[nano seconds]" << std::endl;
}

void processKeys(unsigned char c, int xx, int yy) {
	double step = 1.0;
	

	int countPoints = 0;
	
	std::chrono::steady_clock::time_point begin, end;
	Point* p;

	double sizex = XMAX - XMIN;
	if (sizex < 0)
		sizex = sizex * (-1);

	double sizez = YMAX - YMIN;
	if (sizez < 0)
		sizez = sizez * (-1);

	switch (c)
	{
	case 'o':
		start = !start;
		break;
	case 'w':
		center.y += step;
		break;
	case 'd':
		center.x += step;
		break;
	case 's':center.y -= step;
		break;
	case 'a':
		center.x -= step;
		break;
	case 'q':
		center.z += step;
		break;
	case 'e':
		center.z -= step;
		break;
	case 'k':
		size_query += step;
		break;
	case 'j':
		size_query -= step;
		break;
	case 'r':
		restart();
		break;
	case 'p':

		for (float x = -particleRadius * points_p; x <= particleRadius * points_p; x += particleDiameter) {
			for (float y = -particleRadius * points_p; y <= particleRadius * points_p; y += particleDiameter) {
				for (float z = -particleRadius * points_p; z <= particleRadius * points_p; z += particleDiameter) {
					p = new Point(x, y, z);
					points.push_back(p);
					countPoints++;
				}
			}
		}

		//separei a parte de inserir os pontos porque quero medir apenas a parte de inserir.
		//no programa ele vai gerar os pontos 1 vez e vai inserir numa nova arvore cada frame
		begin = std::chrono::steady_clock::now();
		for (int i = 0; i < points_p; i++) {
			o->insertPoint(points.at(pontos_inseridos + i));
		}

		end = std::chrono::steady_clock::now();
		pontos_inseridos += points_p;
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[micro s]" << std::endl;
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[nano s]" << std::endl;

		break;
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
	spherical2Cartesian();
	glutPostRedisplay();
}

void changeSize(int w, int h) {
	// Prevent a divide by zero, when window is too short
	// (you cant make a window with zero width).
	if (h == 0)
		h = 1;

	// compute window's aspect ratio 
	double ratio = w * 1.0 / h;

	// Set the projection matrix as current
	glMatrixMode(GL_PROJECTION);
	// Load Identity Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set perspective
	gluPerspective(45.0f, ratio, 1.0f, 1000.0f);
	TwWindowSize(w, h);
	// return to the model view matrix mode
	glMatrixMode(GL_MODELVIEW);
}

/**
 * @brief draw an x,y,z axis for better visualization of scenario
 *
 */
void drawAxis() {
	// red x
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(-4.0, 0.0f, 0.0f);
	glVertex3f(4.0, 0.0f, 0.0f);

	// arrow
	glVertex3f(4.0, 0.0f, 0.0f);
	glVertex3f(3.0, 1.0f, 0.0f);

	glVertex3f(4.0, 0.0f, 0.0f);
	glVertex3f(3.0, -1.0f, 0.0f);

	glEnd();
	glFlush();

	// green y
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0.0, -4.0f, 0.0f);
	glVertex3f(0.0, 4.0f, 0.0f);

	// arrow
	glVertex3f(0.0, 4.0f, 0.0f);
	glVertex3f(1.0, 3.0f, 0.0f);

	glVertex3f(0.0, 4.0f, 0.0f);
	glVertex3f(-1.0, 3.0f, 0.0f);

	glEnd();
	glFlush();

	// blue z
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(0.0, 0.0f, -4.0f);
	glVertex3f(0.0, 0.0f, 4.0f);

	// arrow
	glVertex3f(0.0, 0.0f, 4.0f);
	glVertex3f(0.0, 1.0f, 3.0f);

	glVertex3f(0.0, 0.0f, 4.0f);
	glVertex3f(0.0, -1.0f, 3.0f);

	glEnd();
	glFlush();
}

void drawQueryVolume() {
	glPushMatrix();

	glColor3f(1.0f, 0.0f, 0.0f);
	glTranslatef(center.x, center.y, center.z);
	glutWireCube(H);

	glPopMatrix();
}

void updatePoints() {
	double x, y, z;
	x = o->pos.x;
	y = o->pos.y;
	z = o->pos.z;
	delete o;
	o = new Octree(x, y, z, size, max_points, max_depth);
	for (size_t i = 0; i < points.size(); i++)
	{
		//points.at(i)->force.x = 0;
		//points.at(i)->force.y = 0;
		//points.at(i)->force.z = 0;
		o->insertPoint(points.at(i));

	}
}

//-------------------V 2
//Kernels
double useDefaultKernel(glm::vec3 distVector, double supportRadius) {
	double dist = glm::length(distVector);
	
	if (dist > supportRadius) {
		
		return 0.0;
	}
	else {
		return (315 / (64 * M_PI * powf(supportRadius, 9.0f))) * powf(supportRadius * supportRadius - dist * dist, 3.0f);
	}
}

void useDefaultKernel_gradient(glm::vec3 distVector, double supportRadius,glm::vec3 * ret) {
	double dist = glm::length(distVector);
	if (dist > supportRadius) {
		ret->x = 0;
		ret->y = 0;
		ret->z = 0;
	}
	else {
		ret->x = -(distVector.x * (945 / (32 * M_PI * powf(supportRadius, 9.0f))) * powf(supportRadius * supportRadius - dist * dist, 2.0f));
		ret->y = -(distVector.y * (945 / (32 * M_PI * powf(supportRadius, 9.0f))) * powf(supportRadius * supportRadius - dist * dist, 2.0f));
		ret->z = -(distVector.z * (945 / (32 * M_PI * powf(supportRadius, 9.0f))) * powf(supportRadius * supportRadius - dist * dist, 2.0f));
	}
}

double useDefaultKernel_laplacian(glm::vec3 distVector, double supportRadius) {
	double dist = glm::length(distVector);
	if (dist > supportRadius)
		return 0.0f;
	else
		return -(945 / (32 * M_PI * powf(supportRadius, 9.0f))) * (supportRadius * supportRadius - dist * dist) * (3 * supportRadius * supportRadius - 7 * dist * dist);
}

void usePressureKernel_gradient(glm::vec3 distVector, double supportRadius, glm::vec3 * ret) {
	double dist = glm::length(distVector);
	if (dist > supportRadius) {
		
		ret->x = 0;
		ret->y = 0;
		ret->z = 0;
	}
	else if (dist < 10e-5) // If ||r|| -> 0+
	{
		
		glm::vec3 normalized = glm::normalize(glm::vec3(1.0f, 1.0f, 1.0f));
		ret->x= -(normalized.x * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
		ret->y = -(normalized.y * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
		ret->z = -(normalized.z * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
	}
	else
	{
		
		glm::vec3 normalized = glm::normalize(distVector);
		ret->x= -(normalized.x * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
		ret->y = -(normalized.y * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
		ret->z = -(normalized.z * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
	}
}

double useViscosityKernel_laplacian(glm::vec3 distVector, double supportRadius) {
	double dist = glm::length(distVector);
	if (dist > supportRadius)
		return 0.0f;
	else
		return (45 / (M_PI * powf(supportRadius, 6.0f))) * (supportRadius - dist);
}

double computeDensity(glm::vec3 position) {
	double sum = 0;
	//ainda nao usa as query. Vai buscar todas as particulas
	for each (Point * p in points)
	{
		//O meu kernel da 0 em todos
		//printf("Position %f %f %f Ponto2 %f %f %f\n", position.x, position.y, position.z, p->pos.x, p->pos.y, p->pos.z);
		//printf("Mass %f default kernel %f H %f\n", MASS, useDefaultKernel(position - p->pos, H), H);
		sum += MASS * useDefaultKernel(position - p->pos, H);
	}
	
	return sum;
}

double computePressure(double density) {
	return STIFF * (density - RESTDENSITY);
}

void computeForce( double density, double pressure, glm::vec3 position, glm::vec3 * ret) {
	glm::vec3 sum=glm::vec3(0.0f, 0.0f, 0.0f);

	for each (Point * p in points)
	{
		if (position.x == p->pos.x && position.y == p->pos.y && position.z == p->pos.z)
			continue;
		glm::vec3 ret = glm::vec3(0);
		usePressureKernel_gradient(position - p->pos, H, &ret);

		
		
		sum.x += ret.x * (pressure / (density * density) + p->pressure / (p->density * p->density)) * MASS;
		sum.y += ret.y * (pressure / (density * density) + p->pressure / (p->density * p->density)) * MASS;
		sum.z += ret.z * (pressure / (density * density) + p->pressure / (p->density * p->density)) * MASS;
	}
	
	ret->x = -(sum.x * density);
	ret->y = -(sum.y * density);
	ret->z = -(sum.z * density);
	
}

void computeViscosity( glm::vec3 velocity, glm::vec3 position,glm::vec3 * ret) {
	glm::vec3 sum= glm::vec3(0.0f, 0.0f, 0.0f);

	for each (Point *p in points)
	{
		if (position.x == p->pos.x && position.y == p->pos.y && position.z == p->pos.z)
			continue;
		sum.x += (p->velocity.x - velocity.x) * (MASS / p->density) * useViscosityKernel_laplacian(position - p->pos, H);
		sum.y += (p->velocity.y - velocity.y) * (MASS / p->density) * useViscosityKernel_laplacian(position - p->pos, H);
		sum.z += (p->velocity.z - velocity.z) * (MASS / p->density) * useViscosityKernel_laplacian(position - p->pos, H);
	}
	ret->x = sum.x*VISCOSITY;
	ret->y = sum.y * VISCOSITY;
	ret->z = sum.z * VISCOSITY;
}

void computeGravity(double density,glm::vec3 * ret) {
	ret->x= GRAVITY.x * density;
	ret->y = GRAVITY.y * density;
	ret->z = GRAVITY.z * density;
}

void computeSurfaceNormal(glm::vec3 position,glm::vec3 * ret) {
	glm::vec3 sum= glm::vec3(0.0f, 0.0f, 0.0f);

	for each (Point * p in points)
	{
		glm::vec3 ret2 = glm::vec3(0);
		useDefaultKernel_gradient(position - p->pos, H, &ret2);
		sum.x +=  ret2.x * (MASS / p->density);
		sum.y += ret2.y * (MASS / p->density);
		sum.z += ret2.z * (MASS / p->density);
	}
	ret->x = sum.x;
	ret->y = sum.y;
	ret->z = sum.z;
}

void computeSurfaceTension(glm::vec3 surfaceNormal, glm::vec3 position, glm::vec3 * ret) {
	double sum = 0.0f;

	for each (Point *p in points)
	{
		sum += (MASS / p->density) * useDefaultKernel_laplacian(position - p->pos, H);
	}

	glm::vec3 surfaceNormalNormalized = glm::normalize(surfaceNormal);
	ret->x = -(surfaceNormalNormalized.x * SURFACETENSION * sum);
	ret->y = -(surfaceNormalNormalized.y * SURFACETENSION * sum);
	ret->z = -(surfaceNormalNormalized.z * SURFACETENSION * sum);
	
}

bool detectCollision(Point * particle, glm::vec3 * contactPoint, glm::vec3 * unitSurfaceNormal) {
	if (particle->pos.x <= XMAX && particle->pos.x >= XMIN && particle->pos.y <= YMAX && particle->pos.y >= YMIN && particle->pos.z <= ZMAX && particle->pos.z >= ZMIN)
		return false;

	char maxComponent = 'x';
	float maxValue = abs(particle->pos.x);
	if (maxValue < abs(particle->pos.y)) {
		maxComponent = 'y';
		maxValue = abs(particle->pos.y);
	}
	if (maxValue < abs(particle->pos.z)) {
		maxComponent = 'z';
		maxValue = abs(particle->pos.z);
	}
	// 'unitSurfaceNormal' is based on the current position component with the largest absolute value
	switch (maxComponent) {
	case 'x':
		if (particle->pos.x < XMIN) {
			contactPoint->x = XMIN;
			contactPoint->y = particle->pos.y; 
			contactPoint->z = particle->pos.z;

			if (particle->pos.y < YMIN)     contactPoint->y = YMIN;
			else if (particle->pos.y > YMAX) contactPoint->y = YMAX;
			if (particle->pos.z < ZMIN)     contactPoint->z = ZMIN;
			else if (particle->pos.z > ZMAX) contactPoint->z = ZMAX;
			
			unitSurfaceNormal->x = 1;
			unitSurfaceNormal->y = 0;
			unitSurfaceNormal->z = 0;
		}
		else if (particle->pos.x > XMAX) {
			contactPoint->x = XMAX;
			contactPoint->y = particle->pos.y;
			contactPoint->z = particle->pos.z;

			if (particle->pos.y < YMIN)     contactPoint->y = YMIN;
			else if (particle->pos.y > YMAX) contactPoint->y = YMAX;
			if (particle->pos.z < ZMIN)     contactPoint->z = ZMIN;
			else if (particle->pos.z > ZMAX) contactPoint->z = ZMAX;
			
			unitSurfaceNormal->x = -1;
			unitSurfaceNormal->y = 0;
			unitSurfaceNormal->z = 0;

		}
		break;
	case 'y':
		if (particle->pos.y < YMIN) {
			contactPoint->x = particle->pos.x;
			contactPoint->y = YMIN;
			contactPoint->z = particle->pos.z;

			if (particle->pos.x < XMIN)     contactPoint->x = XMIN;
			else if (particle->pos.x > XMAX) contactPoint->x = XMAX;
			if (particle->pos.z < ZMIN)     contactPoint->z = ZMIN;
			else if (particle->pos.z > ZMAX) contactPoint->z = ZMAX;

			unitSurfaceNormal->x = 0;
			unitSurfaceNormal->y = 1;
			unitSurfaceNormal->z = 0;
		}
		else if (particle->pos.y > YMAX) {
			contactPoint->x = particle->pos.x;
			contactPoint->y = YMAX;
			contactPoint->z = particle->pos.z;

			if (particle->pos.x < XMIN)     contactPoint->x = XMIN;
			else if (particle->pos.x > XMAX) contactPoint->x = XMAX;
			if (particle->pos.z < ZMIN)     contactPoint->z = ZMIN;
			else if (particle->pos.z > ZMAX) contactPoint->z = ZMAX;

			unitSurfaceNormal->x = 0;
			unitSurfaceNormal->y = -1;
			unitSurfaceNormal->z = 0;
		}
		break;
	case 'z':
		if (particle->pos.z < ZMIN) {
			contactPoint->x = particle->pos.x;
			contactPoint->y = particle->pos.y;
			contactPoint->z = ZMIN;

			if (particle->pos.x < XMIN)     contactPoint->x = XMIN;
			else if (particle->pos.x > XMAX) contactPoint->x = XMAX;
			if (particle->pos.y < YMIN)     contactPoint->y = YMIN;
			else if (particle->pos.y > YMAX) contactPoint->y = YMAX;
			unitSurfaceNormal->x = 0;
			unitSurfaceNormal->y = 0;
			unitSurfaceNormal->z = 1;
		}
		else if (particle->pos.z > ZMAX) {
			contactPoint->x = particle->pos.x;
			contactPoint->y = particle->pos.y;
			contactPoint->z = ZMAX;

			if (particle->pos.x < XMIN)     contactPoint->x = XMIN;
			else if (particle->pos.x > XMAX) contactPoint->x = XMAX;
			if (particle->pos.y < YMIN)     contactPoint->y = YMIN;
			else if (particle->pos.y > YMAX) contactPoint->y = YMAX;
			unitSurfaceNormal->x = 0;
			unitSurfaceNormal->y = 0;
			unitSurfaceNormal->z = -1;
		}
		break;
	}
	
	return true;
}

void updateVelocity(glm::vec3 velocity, glm::vec3 unitSurfaceNormal, float penetrationDepth, glm::vec3 * ret) {
	//ret = velocity - unitSurfaceNormal * (1 + RESTITUTION * penetrationDepth / (TIMESTEP * glm::length(velocity))) * glm::dot(velocity, unitSurfaceNormal);
	
	ret->x = velocity.x - unitSurfaceNormal.x * (1 + RESTITUTION * penetrationDepth / (TIMESTEP * glm::length(velocity))) * glm::dot(unitSurfaceNormal, velocity);
	ret->y = velocity.y - unitSurfaceNormal.y * (1 + RESTITUTION * penetrationDepth / (TIMESTEP * glm::length(velocity))) * glm::dot(unitSurfaceNormal, velocity);
	ret->z = velocity.z - unitSurfaceNormal.z * (1 + RESTITUTION * penetrationDepth / (TIMESTEP * glm::length(velocity))) * glm::dot(unitSurfaceNormal, velocity);
}

void simulate() {
	
	//density and pressure
	for each (Point * p in points) {
		p->density = computeDensity(p->pos);
		p->pressure = computePressure(p->density);
	}
	//Internal forces
	for each (Point * p in points) {
		glm::vec3 ret = glm::vec3(0);
		computeForce(p->density,p->pressure,p->pos,&ret);
		//printf("ret 1 %f %f %f\n", ret.x, ret.y, ret.z);
		p->force.x = ret.x;
		p->force.y = ret.y;
		p->force.z = ret.z;
		//printf("Force %f %f %f\n", p->force.x, p->force.y, p->force.z);

		computeViscosity( p->velocity, p->pos,&ret);
		//printf("ret 2 %f %f %f\n", ret.x, ret.y, ret.z);
		p->viscosity.x = ret.x;
		p->viscosity.y = ret.y;
		p->viscosity.z = ret.z;
		//printf("Viscosity %f %f %f\n", p->viscosity.x, p->viscosity.y, p->viscosity.z);
	}
	
	//external forces
	for each (Point * p in points)
	{
		glm::vec3 ret = glm::vec3(0);
		computeGravity(p->density,&ret);
		p->gravity.x = ret.x;
		p->gravity.y = ret.y;
		p->gravity.z = ret.z;

		computeSurfaceNormal(p->pos, &ret);
		p->surfaceNormal.x = ret.x;
		p->surfaceNormal.y = ret.y;
		p->surfaceNormal.z = ret.z;
		if (glm::length(p->surfaceNormal) >= THRESHOLD) {
			computeSurfaceTension(p->surfaceNormal, p->pos, &ret);
			p->surfaceTension.x = ret.x;
			p->surfaceTension.y = ret.y;
			p->surfaceTension.z = ret.z;
		}
		else
			p->surfaceTension = glm::vec3(0);
	}

	// Time integration and collision handling
	static double time = 0.0f;
	time += TIMESTEP;
	glm::vec3 totalForce=glm::vec3(0);

	for each (Point * p in points)
	{
		//printf("PressureForce %f %f %f Viscosity %f %f %f Gravity %f %f %f surface tension %f %f %f\n", mParticles[i].mPressureForce.x, mParticles[i].mPressureForce.y, mParticles[i].mPressureForce.z, mParticles[i].mViscosityForce.x, mParticles[i].mViscosityForce.y, mParticles[i].mViscosityForce.z, mParticles[i].mGravitationalForce.x, mParticles[i].mGravitationalForce.y, mParticles[i].mGravitationalForce.z, mParticles[i].mSurfaceTensionForce.x, mParticles[i].mSurfaceTensionForce.y, mParticles[i].mSurfaceTensionForce.z);
		totalForce.x = p->force.x + p->viscosity.x + p->gravity.x + p->surfaceTension.x;
		totalForce.y = p->force.y + p->viscosity.y + p->gravity.y + p->surfaceTension.y;
		totalForce.z = p->force.z + p->viscosity.z + p->gravity.z + p->surfaceTension.z;
		//employEulerIntegrator
		//printf("TotalFrce %f %f %f density %f\n", totalForce.x, totalForce.y, totalForce.z, p->density);
		p->acceleration.x = totalForce.x / p->density;
		p->acceleration.y = totalForce.y / p->density;
		p->acceleration.z = totalForce.z / p->density;

		//printf("Acceleration %f %f %f\n", p->acceleration.x, p->acceleration.y, p->acceleration.z);
		p->velocity.x = p->velocity.x + p->acceleration.x * TIMESTEP;
		p->velocity.y = p->velocity.y + p->acceleration.y * TIMESTEP;
		p->velocity.z = p->velocity.z + p->acceleration.z * TIMESTEP;
		//printf("Antes Position %f %f %f velocity %f %f %f\n", p->pos.x, p->pos.y, p->pos.z, p->velocity.x, p->velocity.y, p->velocity.z);
		p->pos.x = p->pos.x + p->velocity.x * TIMESTEP;
		p->pos.y = p->pos.y + p->velocity.y * TIMESTEP;
		p->pos.z = p->pos.z + p->velocity.z * TIMESTEP;

		//printf("Position %f %f %f\n", p->pos.x, p->pos.y, p->pos.z);
		
		glm::vec3 contactPoint = glm::vec3(0);
		glm::vec3 unitSurfaceNormal=glm::vec3(0);
		if (detectCollision(p, &contactPoint, &unitSurfaceNormal)) {
			//printf("unitSurfaceNormal %f %f %f\n", unitSurfaceNormal.x, unitSurfaceNormal.y, unitSurfaceNormal.z);
			glm::vec3 ret = glm::vec3(0);
			updateVelocity(p->velocity, unitSurfaceNormal, glm::length((p->pos - contactPoint)),&ret);
			//printf("Antes --- Velocity %f %f %f length %f\n", p->velocity.x, p->velocity.y, p->velocity.z, glm::length(p->velocity));
			p->velocity.x = ret.x;
			p->velocity.y = ret.y;
			p->velocity.z = ret.z;
			//printf("Depois --- Velocity %f %f %f length %f\n", p->velocity.x, p->velocity.y, p->velocity.z, glm::length(p->velocity));
			p->pos = contactPoint;
		}

	}

}

void drawBox() {
	// Draw bottom surface edges of the box
	glBegin(GL_LINE_LOOP);
	glVertex3f(XMIN, YMIN, ZMIN);
	glVertex3f(XMAX, YMIN, ZMIN);
	glVertex3f(XMAX, YMIN, ZMAX);
	glVertex3f(XMIN, YMIN, ZMAX);
	glEnd();

	// Draw top surface edges of the box
	glBegin(GL_LINE_LOOP);
	glVertex3f(XMIN, YMAX, ZMIN);
	glVertex3f(XMAX, YMAX, ZMIN);
	glVertex3f(XMAX, YMAX, ZMAX);
	glVertex3f(XMIN, YMAX, ZMAX);
	glEnd();

	// Draw left surface edges of the box
	glBegin(GL_LINE_LOOP);
	glVertex3f(XMIN, YMAX, ZMIN);
	glVertex3f(XMIN, YMAX, ZMAX);
	glVertex3f(XMIN, YMIN, ZMAX);
	glVertex3f(XMIN, YMIN, ZMIN);
	glEnd();

	// Draw right surface edges of the box
	glBegin(GL_LINE_LOOP);
	glVertex3f(XMAX, YMAX, ZMIN);
	glVertex3f(XMAX, YMAX, ZMAX);
	glVertex3f(XMAX, YMIN, ZMAX);
	glVertex3f(XMAX, YMIN, ZMIN);
	glEnd();
}

void renderScene(void) {
	char fpss[200];

	static double t = 0;

	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set the camera
	glLoadIdentity();
	gluLookAt(camX, camY, camZ,
		(XMIN+XMAX)/2, (YMIN + YMAX) / 2, (ZMIN + ZMAX) / 2,
		0.0f, 1.0f, 0.0f);
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
	
	if (start)
		simulate();

	//Draw stuff
	//drawAxis();

	//Desenha o volume da octree (apenas o maior cubo)
	o->drawOut();

	//Calcula o raio das esferas
	double sphereRadius = powf((3 * MASS) / (4 * M_PI * RESTDENSITY), 1.0f / 3.0f);
	
	for (size_t i = 0; i < points.size(); i++)
	{
		points.at(i)->draw(sphereRadius);
	}

	drawBox();

	glColor3f(1.0f, 0.0f, 0.0f);

	glPointSize(5.0f);
	glBegin(GL_POINTS);
	
	glVertex3f(0,0,0);
	glEnd();

	TwDraw();

	glEnable(GL_LIGHTING);

	// End of frame
	glutSwapBuffers();
}

void SetupNormal() {
	//definir o volume da simulação
	XMIN = -0.2;
	XMAX = 0.2;
	YMIN = -0.2;
	YMAX = 0.2;
	ZMIN = -0.2;
	ZMAX = 0.2;

	GRAVITY = glm::vec3(0.0f, -9.8f, 0.0f);

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
	o = new Octree(0, 0, 0, size, max_points, max_depth);

	//SetupDualWave();
	SetupNormal();


	TwBar* bar;         // Pointer to a tweak bar

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

	// Initialize AntTweakBar
	TwInit(TW_OPENGL, NULL);

	glutKeyboardFunc(processKeys);
	glutSpecialFunc(processSpecialKeys);
	//glutKeyboardFunc((GLUTkeyboardfun)TwEventKeyboardGLUT);
	//glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);


	glutMouseFunc((GLUTmousebuttonfun)TwEventMouseButtonGLUT);
	glutMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT);
	glutPassiveMotionFunc((GLUTmousemotionfun)TwEventMouseMotionGLUT); // same as MouseMotion

	TwGLUTModifiersFunc(glutGetModifiers);

	//Bar
// Create a tweak bar
	bar = TwNewBar("TweakBar");
	TwDefine(" TweakBar size='200 400' color='96 216 224' "); // change default tweak bar size and color

	// Add 'wire' to 'bar': it is a modifable variable of type TW_TYPE_BOOL32 (32 bits boolean). Its key shortcut is [w].

	TwAddVarRW(bar, "TimeStep", TW_TYPE_DOUBLE, &TIMESTEP,
		" label='TimeStep' min=0.001 max=2 step=0.001 keyIncr=s keyDecr=S help='Time step used in the simulation' ");

	TwAddVarRW(bar, "Particles", TW_TYPE_INT32, &points_p,
		" label='Particles' min=1 max=10000000 step=1000 key=w help='Number of particles.' ");

	TwAddVarRW(bar, "Declive", TW_TYPE_DOUBLE, &DECLIVE,
		" label='Declive' min=0.0 max=1 step=0.05 keyIncr=s keyDecr=S help='Declive' ");

	TwAddVarRW(bar, "Gravity", TW_TYPE_DOUBLE, &GRAVITYVALUE,
		" label='Gravity' min=-200 max=200 step=1 keyIncr=s keyDecr=S help='Gravity' ");
	

	TwAddVarRW(bar, "Viscosity", TW_TYPE_DOUBLE, &VISCOSITY,
		" label='Viscosity' min=0.1 max=20 step=0.1 keyIncr=s keyDecr=S help='Fluid Viscosity' ");

	TwAddVarRW(bar, "Rest Density", TW_TYPE_DOUBLE, &RESTDENSITY,
		" label='Rest Density' min=100 max=10000 step=100 keyIncr=s keyDecr=S help='Rest density' ");

	TwAddVarRW(bar, "Restitution", TW_TYPE_DOUBLE, &RESTITUTION,
		" label='Restitution' min=0.1 max=2 step=0.05 keyIncr=s keyDecr=S help='Determines the velocity 'lost' when a particle hits an object' ");

	TwAddVarRW(bar, "SurfaceTension", TW_TYPE_DOUBLE, &SURFACETENSION,
		" label='Surface Tension' min=0.0001 max=2 step=0.0001 keyIncr=s keyDecr=S help='Determines the surface tension of the fluid' ");

	TwAddVarRW(bar, "SurfaceTensionThreshold", TW_TYPE_DOUBLE, &THRESHOLD,
		" label='Surface Tension Threshold' min=1 max=100 step=0.005 keyIncr=s keyDecr=S help='Threshold where surface tension is calculated' ");
	
	TwAddVarRW(bar, "KernelSize", TW_TYPE_DOUBLE, &H,
		" label='Kernel Size' min=0.0001 max=100 step=0.0001 keyIncr=s keyDecr=S help='Kernel size' ");

	TwAddVarRW(bar, "Stiff", TW_TYPE_DOUBLE, &STIFF,
		" label='Stiff' min=1 max=100 step=1 keyIncr=s keyDecr=S help='Determines how close the particles of the simulation can be (1- close 5- they start pushing other particles)' ");

	TwAddVarRW(bar, "Mass", TW_TYPE_DOUBLE, &MASS ,
		" label='Mass' min=0.01 max=100 step=0.01 keyIncr=s keyDecr=S help='Mass of the particle' ");
	
	
	
	
	
	double H = 0.0457; // Com o H a 10 ja apanha algumas particulas vizinhas
	double STIFF = 3.0;
	double DECLIVE = 0.0;
	double MASS = 0.02;
	double RESTDENSITY = 998.29;
	double VISCOSITY = 3.5;

	// OpenGL settings 
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	setLighting();
	// enter GLUT's main loop
	glutMainLoop();

	return 1;
}