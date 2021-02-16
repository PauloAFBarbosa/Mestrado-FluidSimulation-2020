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

#include <iostream>
#include <fstream>
#include <iomanip> 
using namespace std;

#endif

#define _USE_MATH_DEFINES
#include <math.h>
#include "HashMap.h"


#define NOQUERY 0
#define QUERYV1 1

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
int pontos_lado = 10;
int countPoints = 0;
//usado para inserir multiplos batches de pontos
int pontos_inseridos = 0;


HashMap* h;
Octree* o;
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

//Da maneira que fiz agora ele verifica 1734 pontos. Em vez de 10 000 pontos.



/**
 * @brief calculates the cam values (used in GluLookAt function) from the alteration of alfa, beta and radius
 *
 */

float length(float vec[3]) {
	return (sqrt(pow(vec[0],2)+ pow(vec[1], 2)+ pow(vec[2], 2)));
}

void normalize(float in[3], float out[3]) {
	
	float len = length(in);
	out[0] = in[0] / len;
	out[1] = in[1] / len;
	out[2] = in[2] / len;
	
}

float dot(float v1[3], float v2[3]) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void spherical2Cartesian() {
	camX = radius * cos(beta) * sin(alfa) + (XMIN+XMAX)/2;
	camY = radius * sin(beta) + (YMIN + YMAX) / 2;
	camZ = radius * cos(beta) * cos(alfa) + (ZMIN + ZMAX) / 2;
}

void restart() {
	std::vector<Point*> points; //Guarda os pontos temporariamente
	pontos_inseridos = 0;
	int countPoints = 0;
	std::chrono::steady_clock::time_point begin, end;
	int x_pontos = 0;
	int y_pontos = 0;
	int z_pontos = 0;
	

	int bucketSizes[hashSize];

	for (size_t i = 0; i < hashSize; i++)
	{
		bucketSizes[i] = 0;
	}
	
	
	Point* p;
	
	
	
	for (float x = -particleRadius * pontos_lado; x_pontos < pontos_lado; x += particleDiameter, x_pontos++) {
		for (float y = -particleRadius * pontos_lado, y_pontos = 0; y_pontos < pontos_lado; y += particleDiameter, y_pontos++) {
			for (float z = -particleRadius * pontos_lado, z_pontos = 0; z_pontos < pontos_lado; z += particleDiameter, z_pontos++) {
				p = new Point(x, y, z, countPoints);
				points.push_back(p);
				bucketSizes[h->hashFunction(p->pos, H, hashSize)]++;
				
				countPoints++;
			}
		}
	}

	h = new HashMap(2, hashSize, countPoints);


	begin = std::chrono::steady_clock::now();
	for each (Point * p in points)
	{
		float pos[3];
		pos[0] = p->pos[0];
		pos[1] = p->pos[1];
		pos[2] = p->pos[2];
		unsigned int bucket = h->hashFunction(pos, H,hashSize);
		int offset = 0;
		//vai correr o array com os tamanhos dos buckets até chegar ao bucket atual
		for (unsigned int j = 0; j < bucket; j++)
		{
			//exemplo, se tivermos um array com [p1,p2,p3,p4]
			//p1 e p2 sao do bucket 0 e p3 e p4 do bucket 1
			//o bucket 1 começa no indice 2 que é o bucketsize do bucket 0
			offset += bucketSizes[j];
		}

		h->addParticle(p, H, offset);

	}
	h->computeOffsets();

	end = std::chrono::steady_clock::now();
	
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[micro seconds]" << std::endl;
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[nano seconds]" << std::endl;
}

void processKeys(unsigned char c, int xx, int yy) {
	ofstream hashSizeFile;
	ofstream positionFile;
	float step = 1.0;
	
	int x_pontos = 0;
	int y_pontos = 0;
	int z_pontos = 0;
	
	std::chrono::steady_clock::time_point begin, end;
	Point* p;

	//hashmap
	int bucketSizes[hashSize];

	for (size_t i = 0; i < hashSize; i++)
	{
		bucketSizes[i] = 0;
	}

	std::vector<Point*> points; //Guarda os pontos temporariamente
	int total = 0;

	float sizex = XMAX - XMIN;
	if (sizex < 0)
		sizex = sizex * (-1);

	float sizez = YMAX - YMIN;
	if (sizez < 0)
		sizez = sizez * (-1);

	switch (c)
	{
	case 'o':
		start = !start;
		break;
	case 'w':
		center[1] += step;
		break;
	case 'd':
		center[0] += step;
		break;
	case 's':center[1] -= step;
		break;
	case 'a':
		center[0] -= step;
		break;
	case 'q':
		center[2] += step;
		break;
	case 'e':
		center[2] -= step;
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
	case 'f':
		simulateFrame = true;
		break;
	case 'p':
		countPoints = 0;
		
		for (float x = -particleRadius * pontos_lado; x_pontos < pontos_lado ; x += particleDiameter,x_pontos++) {
			for (float y = -particleRadius * pontos_lado, y_pontos = 0; y_pontos < pontos_lado; y += particleDiameter, y_pontos++) {
				for (float z = -particleRadius * pontos_lado, z_pontos = 0; z_pontos < pontos_lado; z += particleDiameter, z_pontos++) {
					p = new Point(x, y, z, countPoints);
					//printf("positions %f %f %f \n", x, y, z);
					points.push_back(p);
					bucketSizes[ h->hashFunction(p->pos,H,hashSize)]++;
					total++;
					countPoints++;
				}
			}
		}
		
		printf("Criei %d pontos\n", countPoints);

		h = new HashMap(2, hashSize, countPoints);
		

		begin = std::chrono::steady_clock::now();
		for each (Point * p in points)
		{
			float pos[3];
			pos[0]=p->pos[0];
			pos[1] = p->pos[1];
			pos[2] = p->pos[2];
			unsigned int bucket = h->hashFunction(pos, H,hashSize);
			int offset = 0;
			//vai correr o array com os tamanhos dos buckets até chegar ao bucket atual
			for (int j = 0; j < bucket; j++)
			{
				//exemplo, se tivermos um array com [p1,p2,p3,p4]
				//p1 e p2 sao do bucket 0 e p3 e p4 do bucket 1
				//o bucket 1 começa no indice 2 que é o bucketsize do bucket 0
				offset += bucketSizes[j];
			}

			h->addParticle(p, H, offset);
			
		}
		
		h->computeOffsets();

		
		hashSizeFile.open("hashSizeFile.txt");
		for (int i = 0; i < hashSize; i++)
		{
			
			hashSizeFile << h->bucketSizes[i]<<"\n";
		}
		hashSizeFile.close();

		printf("Criou o ficheiro hashsizefile\n");
		
		positionFile.open("positionFile.txt");
		positionFile << setprecision(10)<<fixed;
		for (int i = 0; i < countPoints; i++)
		{

			positionFile << h->particles[i].pos[0] <<" "<< h->particles[i].pos[1] << " "<<h->particles[i].pos[2] <<" 0.0"<<"\n";
		}
		positionFile.close();

		printf("Criou o ficheiro positionfile\n");

		end = std::chrono::steady_clock::now();
		
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
	float ratio = w * 1.0 / h;

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
	glTranslatef(center[0], center[1], center[2]);
	glutWireCube(H);

	glPopMatrix();
}


void updatePoints() {
	std::vector<Point*> points; //Guarda os pontos temporariamente
	float x, y, z;
	x = o->pos[0];
	y = o->pos[1];
	z = o->pos[2];
	delete o;
	o = new Octree(x, y, z, 0, max_points, max_depth);
	for (size_t i = 0; i < points.size(); i++)
	{
		//points.at(i)->force[0] = 0;
		//points.at(i)->force[1] = 0;
		//points.at(i)->force[2] = 0;
		o->insertPoint(points.at(i));

	}
}

//-------------------V 2
//Kernels
float useDefaultKernel(float distVector[3], float supportRadius) {
	float dist = length(distVector);
	
	if (dist > supportRadius) {
		
		return 0.0;
	}
	else {
		return (315 / (64 * M_PI * powf(supportRadius, 9.0f))) * powf(supportRadius * supportRadius - dist * dist, 3.0f);
	}
}

void useDefaultKernel_gradient(float distVector[3], float supportRadius,float ret[3]) {
	float dist = length(distVector);
	if (dist > supportRadius) {
		ret[0] = 0;
		ret[1] = 0;
		ret[2] = 0;
	}
	else {
		ret[0] = -(distVector[0] * (945 / (32 * M_PI * powf(supportRadius, 9.0f))) * powf(supportRadius * supportRadius - dist * dist, 2.0f));
		ret[1] = -(distVector[1] * (945 / (32 * M_PI * powf(supportRadius, 9.0f))) * powf(supportRadius * supportRadius - dist * dist, 2.0f));
		ret[2] = -(distVector[2] * (945 / (32 * M_PI * powf(supportRadius, 9.0f))) * powf(supportRadius * supportRadius - dist * dist, 2.0f));
	}
}

float useDefaultKernel_laplacian(float distVector[3], float supportRadius) {
	float dist = length(distVector);
	if (dist > supportRadius)
		return 0.0f;
	else
		return -(945 / (32 * M_PI * powf(supportRadius, 9.0f))) * (supportRadius * supportRadius - dist * dist) * (3 * supportRadius * supportRadius - 7 * dist * dist);
}

void usePressureKernel_gradient(float distVector[3], float supportRadius, float ret[3]) {
	float dist = length(distVector);
	if (dist > supportRadius) {
		
		ret[0] = 0;
		ret[1] = 0;
		ret[2] = 0;
		
	}
	else
	{
		
		float normalized[3];
		normalize(distVector,normalized);
		ret[0]= -(normalized[0] * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
		ret[1] = -(normalized[1] * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
		ret[2] = -(normalized[2] * (45 / (M_PI * powf(supportRadius, 6.0f))) * powf(supportRadius - dist, 2.0f));
		
	}
}

float useViscosityKernel_laplacian(float distVector[3], float supportRadius) {
	float dist = length(distVector);
	if (dist > supportRadius)
		return 0.0f;
	else
		return (45 / (M_PI * powf(supportRadius, 6.0f))) * (supportRadius - dist);
}

float computeDensity(float position[3],int query) {
	float sum = 0;
	int realcount = 0;
	
	if (query == NOQUERY) {
		
		for (int i = 0; i < h->numbParticles; i++)
		{
			float arg[3];
			arg[0] = position[0] - h->particles[i].pos[0];
			arg[1] = position[1] - h->particles[i].pos[1];
			arg[2] = position[2] - h->particles[i].pos[2];
			//if (length(arg) <= H){
			//		realcount++;
			//		printf("no getadj temos a particula com %f %f %f pos\n", h->particles[ i].pos[0], h->particles[ i].pos[1], h->particles[i].pos[2]);
			//}

			sum += MASS * useDefaultKernel(arg, H);
		}
		//printf("Sem usar query, vamos ter %d particulas a contar\n", realcount);
	}
	else if (query == QUERYV1) {
		//Versão com query
		//Point* retPoint = (Point*)malloc(sizeof(Point) * h->numbParticles);
		//int particlesFound = 0;
		//particlesFound=h->getAdjBruteForce(position, H, retPoint);
		//printf("Encontrou %d particulas com o brute force\n", particlesFound);
		//for (int i = 0; i < particlesFound; i++)
		//{
		//	printf("Particula na pos %f %f %f\n", retPoint[i].pos[0], retPoint[i].pos[1], retPoint[i].pos[2]);
		//}
		//printf("Vai fazer o meu getadj\n---------\n");
		unsigned int ret[27];
		int retSize;
		retSize=h->getAdj(position, H, ret);

		int count = 0;
		int dentroDoBucket = 0;
		int dentroDoBucketNaoAfeta = 0;
		for (int j = 0; j < retSize; j++)
		{
			unsigned int bucket = ret[j];
	
			int offset = h->offsets[bucket];
			int bucketSize = h->bucketSizes[bucket];
			for (unsigned int i = 0; i < bucketSize; i++)
			{

				float arg[3];
				arg[0] = position[0] - h->particles[offset + i].pos[0];
				arg[1] = position[1] - h->particles[offset + i].pos[1];
				arg[2] = position[2] - h->particles[offset + i].pos[2];
				
				float ker_res = useDefaultKernel(arg, H);

				if (ker_res != 0) {
					dentroDoBucket++;
				}
				else
				{
					dentroDoBucketNaoAfeta ++ ;
				}

				sum += MASS * ker_res;

				//if (length(arg) <= H){
				//	realcount++;
				//	printf("no getadj temos a particula com %f %f %f pos\n", h->particles[offset + i].pos[0], h->particles[offset + i].pos[1], h->particles[offset + i].pos[2]);
				//}
				
				count++;
			}
			printf("Dentro do bucket %d encontrou %d particulas que afetam e %d que nao afetam\n", bucket, dentroDoBucket, dentroDoBucketNaoAfeta);
			dentroDoBucketNaoAfeta = 0;
			dentroDoBucket = 0;
		}
		printf("Fim de uma particula\n");
		//printf("Com o meu getadj encontrou %d particulas mas apenas %d realmente importaram\n-------\n", count,realcount);
	}
	
	//printf("Sum dentro do computedensity %f\n", sum);
	return sum;
}

float computePressure(float density) {
	return STIFF * (density - RESTDENSITY);
}

void computeForce( float density, float pressure, float position[3], float ret[3],int query) {
	float sum[3] = { 0,0,0 };
	int count = 0;
	int realcount = 0;
	if (query == NOQUERY) {
		for (int i = 0; i < h->numbParticles; i++)
		{
			if (position[0] == h->particles[i].pos[0] && position[1] == h->particles[i].pos[1] && position[2] == h->particles[i].pos[2])
				continue;
			float ret_kernel[3] = { 0,0,0 };

			float arg[3];
			arg[0] = position[0] - h->particles[i].pos[0];
			arg[1] = position[1] - h->particles[i].pos[1];
			arg[2] = position[2] - h->particles[i].pos[2];
			usePressureKernel_gradient(arg, H, ret_kernel);

			//if (ret_kernel[0] == 0 && ret_kernel[1] == 0 && ret_kernel[2] == 0) {

			//}
			//else {
				//printf("Ponto que vai contar %f %f %f\n", h->particles[i].pos[0], h->particles[i].pos[1], h->particles[i].pos[2]);
				//realcount++;
			//}

			//printf("RetKernel %f %f %f \n", ret_kernel[0], ret_kernel[1], ret_kernel[2]);
			
			sum[0] += ret_kernel[0] * (pressure / (density * density) + h->particles[i].pressure / (h->particles[i].density * h->particles[i].density)) * MASS;
			sum[1] += ret_kernel[1] * (pressure / (density * density) + h->particles[i].pressure / (h->particles[i].density * h->particles[i].density)) * MASS;
			sum[2] += ret_kernel[2] * (pressure / (density * density) + h->particles[i].pressure / (h->particles[i].density * h->particles[i].density)) * MASS;
		}
		//printf("Sem query ele tem realcount %d\n", realcount);
		
	}
	else if (query == QUERYV1) {

		//Point* retPoint = (Point*)malloc(sizeof(Point) * h->numbParticles);
		//int particlesFound = 0;
		//particlesFound=h->getAdjBruteForce(position, H, retPoint);
		//printf("Encontrou %d particulas com o brute force\n", particlesFound);
		//for (int i = 0; i < particlesFound; i++)
		//{
		//	printf("Particula na pos %f %f %f\n", retPoint[i].pos[0], retPoint[i].pos[1], retPoint[i].pos[2]);
		//}
		unsigned int ret[27];
		int retSize;
		//printf("vou mandar a posicao %f %f %f\n", position[0], position[1], position[2]);
		retSize=h->getAdj(position, H, ret);
		//printf("retSize deu %d\n", retSize);
		for (int j = 0; j < retSize; j++)
		{
			unsigned int bucket = ret[j];
			int offset = h->offsets[bucket];
			int bucketSize = h->bucketSizes[bucket];
			for (unsigned int i = 0; i < bucketSize; i++)
			{
				if (position[0] == h->particles[offset + i].pos[0] && position[1] == h->particles[offset + i].pos[1] && position[2] == h->particles[offset + i].pos[2])
					continue;
				
				float ret_kernel[3] = { 0,0,0 };

				float arg[3];
				arg[0] = position[0] - h->particles[offset+i].pos[0];
				arg[1] = position[1] - h->particles[offset+i].pos[1];
				arg[2] = position[2] - h->particles[offset+i].pos[2];
				usePressureKernel_gradient(arg, H, ret_kernel);
				//if (ret_kernel[0] == 0 && ret_kernel[1] == 0 && ret_kernel[2] == 0) {

				//}
				//else {
				//	printf("Ponto que vai contar %f %f %f\n", h->particles[offset + i].pos[0], h->particles[offset + i].pos[1], h->particles[offset + i].pos[2]);
				//	realcount++;
				//}
				//printf("RetKernel %f %f %f \n", ret_kernel[0], ret_kernel[1], ret_kernel[2]);
				sum[0] += ret_kernel[0] * (pressure / (density * density) + h->particles[offset+i].pressure / (h->particles[offset + i].density * h->particles[offset + i].density)) * MASS;
				sum[1] += ret_kernel[1] * (pressure / (density * density) + h->particles[offset+i].pressure / (h->particles[offset + i].density * h->particles[offset + i].density)) * MASS;
				sum[2] += ret_kernel[2] * (pressure / (density * density) + h->particles[offset + i].pressure / (h->particles[offset + i].density * h->particles[offset + i].density)) * MASS;
				count++;
			}
		}
		//printf("---O real count do force deu %d\n", realcount);
	}
	//printf("Count %d \n", count);
	
	
	ret[0] = -(sum[0] * density);
	ret[1] = -(sum[1] * density);
	ret[2] = -(sum[2] * density);
	
}

void computeViscosity( float velocity[3], float position[3],float ret[3], int query) {
	float sum[3] = { 0,0,0 };

	if (query == NOQUERY) {
		for (int i = 0; i < h->numbParticles; i++)
		{
			if (position[0] == h->particles[i].pos[0] && position[1] == h->particles[i].pos[1] && position[2] == h->particles[i].pos[2])
				continue;
			float arg[3];
			arg[0] = position[0] - h->particles[i].pos[0];
			arg[1] = position[1] - h->particles[i].pos[1];
			arg[2] = position[2] - h->particles[i].pos[2];

			sum[0] += (h->particles[i].velocity[0] - velocity[0]) * (MASS / h->particles[i].density) * useViscosityKernel_laplacian(arg, H);
			sum[1] += (h->particles[i].velocity[1] - velocity[1]) * (MASS / h->particles[i].density) * useViscosityKernel_laplacian(arg, H);
			sum[2] += (h->particles[i].velocity[2] - velocity[2]) * (MASS / h->particles[i].density) * useViscosityKernel_laplacian(arg, H);
		}
	}
	else if (query == QUERYV1) {
		unsigned int ret[27];
		int retSize;
		//printf("vou mandar a posicao %f %f %f\n", position[0], position[1], position[2]);
		retSize = h->getAdj(position, H, ret);
		//printf("retSize deu %d\n", retSize);
		for (int j = 0; j < retSize; j++)
		{
			unsigned int bucket = ret[j];
			int offset = h->offsets[bucket];
			int bucketSize = h->bucketSizes[bucket];
			for (unsigned int i = 0; i < bucketSize; i++)
			{
				if (position[0] == h->particles[offset+i].pos[0] && position[1] == h->particles[offset + i].pos[1] && position[2] == h->particles[offset + i].pos[2])
					continue;
				float arg[3];
				arg[0] = position[0] - h->particles[offset + i].pos[0];
				arg[1] = position[1] - h->particles[offset + i].pos[1];
				arg[2] = position[2] - h->particles[offset + i].pos[2];

				sum[0] += (h->particles[offset + i].velocity[0] - velocity[0]) * (MASS / h->particles[offset + i].density) * useViscosityKernel_laplacian(arg, H);
				sum[1] += (h->particles[offset + i].velocity[1] - velocity[1]) * (MASS / h->particles[offset + i].density) * useViscosityKernel_laplacian(arg, H);
				sum[2] += (h->particles[offset + i].velocity[2] - velocity[2]) * (MASS / h->particles[offset + i].density) * useViscosityKernel_laplacian(arg, H);
			}
		}
	}
	ret[0] = sum[0]*VISCOSITY;
	ret[1] = sum[1] * VISCOSITY;
	ret[2] = sum[2] * VISCOSITY;
}


//Calcula tudo que usa vizinhos no mesmo ciclo
//Evita fazer a procura de vizinhos muitas vezes
void computeForceViscositySurfaceNormal(float velocity[3], float density, float pressure, float position[3], float retForce[3], float retViscosity[3], float retSurfaceNormal[3], int query) {
	float sumViscosity[3] = { 0,0,0 };
	float sumForce[3] = { 0,0,0 };
	float sumSurfaceNormal[3] = { 0,0,0 };
	int count = 0;
	int realcount = 0;

	if (query == NOQUERY) {
		count = 0;
		for (int i = 0; i < h->numbParticles; i++)
		{
			if (position[0] == h->particles[i].pos[0] && position[1] == h->particles[i].pos[1] && position[2] == h->particles[i].pos[2])
				continue;
			//Viscosity
			float arg[3];

			arg[0] = position[0] - h->particles[i].pos[0];
			arg[1] = position[1] - h->particles[i].pos[1];
			arg[2] = position[2] - h->particles[i].pos[2];

			sumViscosity[0] += (h->particles[i].velocity[0] - velocity[0]) * (MASS / h->particles[i].density) * useViscosityKernel_laplacian(arg, H);
			sumViscosity[1] += (h->particles[i].velocity[1] - velocity[1]) * (MASS / h->particles[i].density) * useViscosityKernel_laplacian(arg, H);
			sumViscosity[2] += (h->particles[i].velocity[2] - velocity[2]) * (MASS / h->particles[i].density) * useViscosityKernel_laplacian(arg, H);

			//Force
			float ret_kernel[3] = { 0,0,0 };
			arg[0] = position[0] - h->particles[i].pos[0];
			arg[1] = position[1] - h->particles[i].pos[1];
			arg[2] = position[2] - h->particles[i].pos[2];
			usePressureKernel_gradient(arg, H, ret_kernel);

			//if (ret_kernel[0] != 0 || ret_kernel[1] != 0 || ret_kernel[2] != 0)
			//	count++;

			sumForce[0] += ret_kernel[0] * (pressure / (density * density) + h->particles[i].pressure / (h->particles[i].density * h->particles[i].density)) * MASS;
			sumForce[1] += ret_kernel[1] * (pressure / (density * density) + h->particles[i].pressure / (h->particles[i].density * h->particles[i].density)) * MASS;
			sumForce[2] += ret_kernel[2] * (pressure / (density * density) + h->particles[i].pressure / (h->particles[i].density * h->particles[i].density)) * MASS;

			//SurfaceNormal
			float ret2[3] = { 0,0,0 };
			
			arg[0] = position[0] - h->particles[i].pos[0];
			arg[1] = position[1] - h->particles[i].pos[1];
			arg[2] = position[2] - h->particles[i].pos[2];
			useDefaultKernel_gradient(arg, H, ret2);
			sumSurfaceNormal[0] += ret2[0] * (MASS / h->particles[i].density);
			sumSurfaceNormal[1] += ret2[1] * (MASS / h->particles[i].density);
			sumSurfaceNormal[2] += ret2[2] * (MASS / h->particles[i].density);
		}
		//printf("Esta particula tem %d vizinhos\n", count);
	}
	else if (query == QUERYV1) {
		unsigned int ret[27];
		int retSize;
		//printf("vou mandar a posicao %f %f %f\n", position[0], position[1], position[2]);
		retSize = h->getAdj(position, H, ret);
		//printf("retSize deu %d\n", retSize);
		count = 0;
		for (int j = 0; j < retSize; j++)
		{
			
			unsigned int bucket = ret[j];
			int offset = h->offsets[bucket];
			int bucketSize = h->bucketSizes[bucket];
			for (unsigned int i = 0; i < bucketSize; i++)
			{
				if (position[0] == h->particles[offset + i].pos[0] && position[1] == h->particles[offset + i].pos[1] && position[2] == h->particles[offset + i].pos[2])
					continue;
				//Viscosity
				float arg[3];
				arg[0] = position[0] - h->particles[offset + i].pos[0];
				arg[1] = position[1] - h->particles[offset + i].pos[1];
				arg[2] = position[2] - h->particles[offset + i].pos[2];

				sumViscosity[0] += (h->particles[offset + i].velocity[0] - velocity[0]) * (MASS / h->particles[offset + i].density) * useViscosityKernel_laplacian(arg, H);
				sumViscosity[1] += (h->particles[offset + i].velocity[1] - velocity[1]) * (MASS / h->particles[offset + i].density) * useViscosityKernel_laplacian(arg, H);
				sumViscosity[2] += (h->particles[offset + i].velocity[2] - velocity[2]) * (MASS / h->particles[offset + i].density) * useViscosityKernel_laplacian(arg, H);
				

				
				//Force
				float ret_kernel[3] = { 0,0,0 };
				arg[0] = position[0] - h->particles[offset + i].pos[0];
				arg[1] = position[1] - h->particles[offset + i].pos[1];
				arg[2] = position[2] - h->particles[offset + i].pos[2];
				usePressureKernel_gradient(arg, H, ret_kernel);

				sumForce[0] += ret_kernel[0] * (pressure / (density * density) + h->particles[offset + i].pressure / (h->particles[offset + i].density * h->particles[offset + i].density)) * MASS;
				sumForce[1] += ret_kernel[1] * (pressure / (density * density) + h->particles[offset + i].pressure / (h->particles[offset + i].density * h->particles[offset + i].density)) * MASS;
				sumForce[2] += ret_kernel[2] * (pressure / (density * density) + h->particles[offset + i].pressure / (h->particles[offset + i].density * h->particles[offset + i].density)) * MASS;
				
				

				//SurfaceNormal
				float ret2[3] = { 0,0,0 };
				
				arg[0] = position[0] - h->particles[offset + i].pos[0];
				arg[1] = position[1] - h->particles[offset + i].pos[1];
				arg[2] = position[2] - h->particles[offset + i].pos[2];
				useDefaultKernel_gradient(arg, H, ret2);
				sumSurfaceNormal[0] += ret2[0] * (MASS / h->particles[offset + i].density);
				sumSurfaceNormal[1] += ret2[1] * (MASS / h->particles[offset + i].density);
				sumSurfaceNormal[2] += ret2[2] * (MASS / h->particles[offset + i].density);
			}
		}
		//printf("Esta particula tem %d vizinhos retSize %d\n", count,retSize);
	}

	
	retViscosity[0] = sumViscosity[0] * VISCOSITY;
	retViscosity[1] = sumViscosity[1] * VISCOSITY;
	retViscosity[2] = sumViscosity[2] * VISCOSITY;

	retForce[0] = -(sumForce[0] * density);
	retForce[1] = -(sumForce[1] * density);
	retForce[2] = -(sumForce[2] * density);

	retSurfaceNormal[0] = sumSurfaceNormal[0];
	retSurfaceNormal[1] = sumSurfaceNormal[1];
	retSurfaceNormal[2] = sumSurfaceNormal[2];
}

void computeGravity(float density,float ret[3]) {
	ret[0]= GRAVITY[0] * density;
	ret[1] = GRAVITY[1] * density;
	ret[2] = GRAVITY[2] * density;
}

void computeSurfaceNormal(float position[3],float ret[3],int query) {
	float sum[3] = { 0,0,0 };
	if (query == NOQUERY) {
		for (int i = 0; i < h->numbParticles; i++)
		{
			float ret2[3] = { 0,0,0 };
			float arg[3];
			arg[0] = position[0] - h->particles[i].pos[0];
			arg[1] = position[1] - h->particles[i].pos[1];
			arg[2] = position[2] - h->particles[i].pos[2];
			useDefaultKernel_gradient(arg, H, ret2);

			sum[0] += ret2[0] * (MASS / h->particles[i].density);
			sum[1] += ret2[1] * (MASS / h->particles[i].density);
			sum[2] += ret2[2] * (MASS / h->particles[i].density);
		}
	}
	else if (query == QUERYV1) {
		unsigned int ret[27];
		int retSize;
		//printf("vou mandar a posicao %f %f %f\n", position[0], position[1], position[2]);
		retSize = h->getAdj(position, H, ret);
		//printf("retSize deu %d\n", retSize);
		for (int j = 0; j < retSize; j++)
		{
			unsigned int bucket = ret[j];
			int offset = h->offsets[bucket];
			int bucketSize = h->bucketSizes[bucket];
			for (unsigned int i = 0; i < bucketSize; i++)
			{
				float ret2[3] = { 0,0,0 };

				float arg[3];
				arg[0] = position[0] - h->particles[offset+i].pos[0];
				arg[1] = position[1] - h->particles[offset + i].pos[1];
				arg[2] = position[2] - h->particles[offset + i].pos[2];
				useDefaultKernel_gradient(arg, H, ret2);

				sum[0] += ret2[0] * (MASS / h->particles[offset + i].density);
				sum[1] += ret2[1] * (MASS / h->particles[offset + i].density);
				sum[2] += ret2[2] * (MASS / h->particles[offset + i].density);
			}
		}
	}
	
	ret[0] = sum[0];
	ret[1] = sum[1];
	ret[2] = sum[2];
}

void computeSurfaceTension(float surfaceNormal[3], float position[3], float ret[3],int query) {
	float sum = 0.0f;
	if (query == NOQUERY) {
		for (int i = 0; i < h->numbParticles; i++)
		{
			float arg[3];
			arg[0] = position[0] - h->particles[i].pos[0];
			arg[1] = position[1] - h->particles[i].pos[1];
			arg[2] = position[2] - h->particles[i].pos[2];
			sum += (MASS / h->particles[i].density) * useDefaultKernel_laplacian(arg, H);
		}
	}
	else if (query == QUERYV1) {
		unsigned int ret[27];
		int retSize;
		//printf("vou mandar a posicao %f %f %f\n", position[0], position[1], position[2]);
		retSize = h->getAdj(position, H, ret);
		//printf("retSize deu %d\n", retSize);
		for (int j = 0; j < retSize; j++)
		{
			unsigned int bucket = ret[j];
			int offset = h->offsets[bucket];
			int bucketSize = h->bucketSizes[bucket];
			for (unsigned int i = 0; i < bucketSize; i++)
			{
				float arg[3];
				arg[0] = position[0] - h->particles[offset + i].pos[0];
				arg[1] = position[1] - h->particles[offset + i].pos[1];
				arg[2] = position[2] - h->particles[offset + i].pos[2];
				sum += (MASS / h->particles[offset + i].density) * useDefaultKernel_laplacian(arg, H);
			}
		}
	}
	

	float surfaceNormalNormalized[3];
	normalize(surfaceNormal, surfaceNormalNormalized);
	ret[0] = -(surfaceNormalNormalized[0] * SURFACETENSION * sum);
	ret[1] = -(surfaceNormalNormalized[1] * SURFACETENSION * sum);
	ret[2] = -(surfaceNormalNormalized[2] * SURFACETENSION * sum);
	
}

bool detectCollision(Point * particle, float contactPoint[3], float unitSurfaceNormal[3]) {
	//o cubo esta centrado em 0,0,0
	float newx = particle->pos[0] + XMAX;
	float temp = (XMAX + XMAX) - newx;
	temp = temp / (XMAX + XMAX); //devolve 1 quando newx é 0, ou seja, a particula esta encostada a parede esquerda
								// devolve 0 quando esta encostada a parede direita

	float newy = YMIN + (temp * DECLIVE);

	if (particle->pos[0] <= XMAX && particle->pos[0] >= XMIN && particle->pos[1] <= YMAX && particle->pos[1] >= newy && particle->pos[2] <= ZMAX && particle->pos[2] >= ZMIN)
		return false;

	char maxComponent = 'x';
	float maxValue = abs(particle->pos[0]);
	//Por causa do declive temos de ter isso em conta ao encontrar o maxvalue. (se nao fizer + temp*declive as vezes da como max component o Z quando na realidade deveria ter sido o Y, so nao foi por causa do declive)
	if (maxValue < abs(particle->pos[1])+(temp*DECLIVE)) {
		maxComponent = 'y';
		maxValue = abs(particle->pos[1]) + (temp * DECLIVE);
	}
	if (maxValue < abs(particle->pos[2])) {
		maxComponent = 'z';
		maxValue = abs(particle->pos[2]);
	}
	// 'unitSurfaceNormal' is based on the current position component with the largest absolute value
	
	switch (maxComponent) {
	case 'x':
		if (particle->pos[0] < XMIN) {
			contactPoint[0] = XMIN;
			contactPoint[1] = particle->pos[1]; 
			contactPoint[2] = particle->pos[2];

			if (particle->pos[1] < newy)     contactPoint[1] = newy;
			else if (particle->pos[1] > YMAX) contactPoint[1] = YMAX;
			if (particle->pos[2] < ZMIN)     contactPoint[2] = ZMIN;
			else if (particle->pos[2] > ZMAX) contactPoint[2] = ZMAX;
			
			unitSurfaceNormal[0] = 1;
			unitSurfaceNormal[1] = 0;
			unitSurfaceNormal[2] = 0;
			
		}
		else if (particle->pos[0] > XMAX) {
			contactPoint[0] = XMAX;
			contactPoint[1] = particle->pos[1];
			contactPoint[2] = particle->pos[2];

			if (particle->pos[1] < newy)     contactPoint[1] = newy;
			else if (particle->pos[1] > YMAX) contactPoint[1] = YMAX;
			if (particle->pos[2] < ZMIN)     contactPoint[2] = ZMIN;
			else if (particle->pos[2] > ZMAX) contactPoint[2] = ZMAX;
			
			unitSurfaceNormal[0] = -1;
			unitSurfaceNormal[1] = 0;
			unitSurfaceNormal[2] = 0;
			

		}
		break;
	case 'y':
		if (particle->pos[1] < newy) {
			contactPoint[0] = particle->pos[0];
			contactPoint[1] = newy;
			contactPoint[2] = particle->pos[2];

			if (particle->pos[0] < XMIN)     contactPoint[0] = XMIN;
			else if (particle->pos[0] > XMAX) contactPoint[0] = XMAX;
			if (particle->pos[2] < ZMIN)     contactPoint[2] = ZMIN;
			else if (particle->pos[2] > ZMAX) contactPoint[2] = ZMAX;

			unitSurfaceNormal[0] = DECLIVE;
			unitSurfaceNormal[1] = 1-DECLIVE;
			unitSurfaceNormal[2] = 0;
			
		}
		else if (particle->pos[1] > YMAX) {
			contactPoint[0] = particle->pos[0];
			contactPoint[1] = YMAX;
			contactPoint[2] = particle->pos[2];

			if (particle->pos[0] < XMIN)     contactPoint[0] = XMIN;
			else if (particle->pos[0] > XMAX) contactPoint[0] = XMAX;
			if (particle->pos[2] < ZMIN)     contactPoint[2] = ZMIN;
			else if (particle->pos[2] > ZMAX) contactPoint[2] = ZMAX;

			unitSurfaceNormal[0] = 0;
			unitSurfaceNormal[1] = -1;
			unitSurfaceNormal[2] = 0;
			
		}
		break;
	case 'z':
		if (particle->pos[2] < ZMIN) {
			contactPoint[0] = particle->pos[0];
			contactPoint[1] = particle->pos[1];
			contactPoint[2] = ZMIN;

			if (particle->pos[0] < XMIN)     contactPoint[0] = XMIN;
			else if (particle->pos[0] > XMAX) contactPoint[0] = XMAX;
			if (particle->pos[1] < newy)     contactPoint[1] = newy;
			else if (particle->pos[1] > YMAX) contactPoint[1] = YMAX;
			unitSurfaceNormal[0] = 0;
			unitSurfaceNormal[1] = 0;
			unitSurfaceNormal[2] = 1;
			
		}
		else if (particle->pos[2] > ZMAX) {
			contactPoint[0] = particle->pos[0];
			contactPoint[1] = particle->pos[1];
			contactPoint[2] = ZMAX;

			if (particle->pos[0] < XMIN)     contactPoint[0] = XMIN;
			else if (particle->pos[0] > XMAX) contactPoint[0] = XMAX;
			if (particle->pos[1] < newy)     contactPoint[1] = newy;
			else if (particle->pos[1] > YMAX) contactPoint[1] = YMAX;
			unitSurfaceNormal[0] = 0;
			unitSurfaceNormal[1] = 0;
			unitSurfaceNormal[2] = -1;
			
		}
		break;
	}

		
	return true;
}

void updateVelocity(float velocity[3], float unitSurfaceNormal[3], float penetrationDepth, float ret[3]) {
	//ret = velocity - unitSurfaceNormal * (1 + RESTITUTION * penetrationDepth / (TIMESTEP * glm::length(velocity))) * glm::dot(velocity, unitSurfaceNormal);
	
	ret[0] = velocity[0] - unitSurfaceNormal[0] * (1 + RESTITUTION * penetrationDepth / (TIMESTEP * length(velocity))) * dot(unitSurfaceNormal, velocity);
	ret[1] = velocity[1] - unitSurfaceNormal[1] * (1 + RESTITUTION * penetrationDepth / (TIMESTEP * length(velocity))) * dot(unitSurfaceNormal, velocity);
	ret[2] = velocity[2] - unitSurfaceNormal[2] * (1 + RESTITUTION * penetrationDepth / (TIMESTEP * length(velocity))) * dot(unitSurfaceNormal, velocity);
}

void simulate() {
	//density and pressure
	for (int i = 0; i < h->numbParticles; i++){
		
		h->particles[i].density = computeDensity(h->particles[i].pos, QUERYV1);
		h->particles[i].pressure = computePressure(h->particles[i].density);
		
	}

	//Internal forces && External forces
	//Calcular as forças internas/externas no mesmo ciclo para ser mais eficiente
	for (int i = 0; i < h->numbParticles; i++) {
		//Em vez de fazer
		//ComputeForce()
		//ComputeDensity()
		//Faz apenas uma função e assim faz apenas uma vez a busca dos vizinhos, e com isso calcula tudo que precisa calcular
		
		float retForce[3] = { 0,0,0 };
		float retViscosity[3] = { 0,0,0 };
		float retSurfaceNormal[3] = { 0,0,0 };
		
		//Esta função calcula a force e viscosity no mesmo ciclo, assim é mais eficiente
		computeForceViscositySurfaceNormal(h->particles[i].velocity, h->particles[i].density, h->particles[i].pressure, h->particles[i].pos, retForce, retViscosity,retSurfaceNormal, QUERYV1);
		
		h->particles[i].force[0] = retForce[0];
		h->particles[i].force[1] = retForce[1];
		h->particles[i].force[2] = retForce[2];

		h->particles[i].viscosity[0] = retViscosity[0];
		h->particles[i].viscosity[1] = retViscosity[1];
		h->particles[i].viscosity[2] = retViscosity[2];

		float retGravity[3] = { 0,0,0 };
		computeGravity(h->particles[i].density, retGravity);
		h->particles[i].gravity[0] = retGravity[0];
		h->particles[i].gravity[1] = retGravity[1];
		h->particles[i].gravity[2] = retGravity[2];

		h->particles[i].surfaceNormal[0] = retSurfaceNormal[0];
		h->particles[i].surfaceNormal[1] = retSurfaceNormal[1];
		h->particles[i].surfaceNormal[2] = retSurfaceNormal[2];

		if (length(h->particles[i].surfaceNormal) >= THRESHOLD) {
			float ret[3] = { 0,0,0 };
			//Tenho de fazer esta função aqui, porque preciso de ter a surface normal calculada
			computeSurfaceTension(h->particles[i].surfaceNormal, h->particles[i].pos, ret, QUERYV1);
			h->particles[i].surfaceTension[0] = ret[0];
			h->particles[i].surfaceTension[1] = ret[1];
			h->particles[i].surfaceTension[2] = ret[2];
		}
		else {
			h->particles[i].surfaceTension[0] = 0;
			h->particles[i].surfaceTension[1] = 0;
			h->particles[i].surfaceTension[2] = 0;
		}
	}

	// Time integration and collision handling
	static float time = 0.0f;
	time += TIMESTEP;
	float totalForce[3] = { 0,0,0 };


	for (int i = 0; i < h->numbParticles; i++) 
	{
		
		totalForce[0] = h->particles[i].force[0] + h->particles[i].viscosity[0] + h->particles[i].gravity[0] + h->particles[i].surfaceTension[0];
		totalForce[1] = h->particles[i].force[1] + h->particles[i].viscosity[1] + h->particles[i].gravity[1] + h->particles[i].surfaceTension[1];
		totalForce[2] = h->particles[i].force[2] + h->particles[i].viscosity[2] + h->particles[i].gravity[2] + h->particles[i].surfaceTension[2];

		
		//employEulerIntegrator
		h->particles[i].acceleration[0] = totalForce[0] / h->particles[i].density;
		h->particles[i].acceleration[1] = totalForce[1] / h->particles[i].density;
		h->particles[i].acceleration[2] = totalForce[2] / h->particles[i].density;

		h->particles[i].velocity[0] = h->particles[i].velocity[0] + h->particles[i].acceleration[0] * TIMESTEP;
		h->particles[i].velocity[1] = h->particles[i].velocity[1] + h->particles[i].acceleration[1] * TIMESTEP;
		h->particles[i].velocity[2] = h->particles[i].velocity[2] + h->particles[i].acceleration[2] * TIMESTEP;

		
		
		h->particles[i].pos[0] = h->particles[i].pos[0] + h->particles[i].velocity[0] * TIMESTEP;
		h->particles[i].pos[1] = h->particles[i].pos[1] + h->particles[i].velocity[1] * TIMESTEP;
		h->particles[i].pos[2] = h->particles[i].pos[2] + h->particles[i].velocity[2] * TIMESTEP;

		float contactPoint[3] = { 0,0,0 };
		float unitSurfaceNormal[3] = { 0,0,0 };
		if (detectCollision(&h->particles[i], contactPoint, unitSurfaceNormal)) {
			
			float ret[3] = { 0,0,0 };
			float arg[3];
			arg[0] = h->particles[i].pos[0] - contactPoint[0];
			arg[1] = h->particles[i].pos[1] - contactPoint[1];
			arg[2] = h->particles[i].pos[2] - contactPoint[2];
			updateVelocity(h->particles[i].velocity, unitSurfaceNormal, length(arg), ret);
			
			h->particles[i].velocity[0] = ret[0];
			h->particles[i].velocity[1] = ret[1];
			h->particles[i].velocity[2] = ret[2];
			

			h->particles[i].pos[0] = contactPoint[0];
			h->particles[i].pos[1] = contactPoint[1];
			h->particles[i].pos[2] = contactPoint[2];
		}
	}
	//printf("---------\n");

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

	static float t = 0;

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
	else if (simulateFrame) {
		printf("vai simular um frame %d\n",framesSimulated);
		auto started = std::chrono::high_resolution_clock::now();
		simulate();
		auto done = std::chrono::high_resolution_clock::now();
		
		std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(done - started).count();

		simulateFrame = false;
		framesSimulated++;
	}
	if (h)
		h->updateHashMap(H);

	//Draw stuff
	//drawAxis();

	//Calcula o raio das esferas
	float sphereRadius = powf((3 * MASS) / (4 * M_PI * RESTDENSITY), 1.0f / 3.0f);
	
	for (int i = 0; h!=NULL && i < h->numbParticles; i++) 
	{
			h->particles[i].draw(sphereRadius); 
			//printf("Render scene Pos da particula %f %f %f\n", h->particles[i].pos[0], h->particles[i].pos[1], h->particles[i].pos[2]);
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
	float side = 2.5;
	XMIN = -side;
	XMAX = side;
	YMIN = -side;
	YMAX = side;
	ZMIN = -side;
	ZMAX = side;

	GRAVITY[0] = 0.0;
	GRAVITY[1] = -9.8;
	GRAVITY[2] = 0.0;

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
	HashMap * h2;
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

	TwAddVarRW(bar, "TimeStep", TW_TYPE_FLOAT, &TIMESTEP,
		" label='TimeStep' min=0.001 max=2 step=0.001 keyIncr=s keyDecr=S help='Time step used in the simulation' ");

	TwAddVarRW(bar, "Particles", TW_TYPE_INT32, &pontos_lado,
		" label='Particles' min=1 max=10000000 step=1000 key=w help='Number of particles ^3 .' ");

	TwAddVarRW(bar, "Declive", TW_TYPE_FLOAT, &DECLIVE,
		" label='Declive' min=-1.0 max=1 step=0.05 keyIncr=s keyDecr=S help='Declive' ");

	TwAddVarRW(bar, "Gravity", TW_TYPE_FLOAT, &GRAVITYVALUE,
		" label='Gravity' min=-200 max=200 step=1 keyIncr=s keyDecr=S help='Gravity' ");
	

	TwAddVarRW(bar, "Viscosity", TW_TYPE_FLOAT, &VISCOSITY,
		" label='Viscosity' min=3 max=25 step=0.1 keyIncr=s keyDecr=S help='Fluid Viscosity' ");

	TwAddVarRW(bar, "Rest Density", TW_TYPE_FLOAT, &RESTDENSITY,
		" label='Rest Density' min=100 max=10000 step=100 keyIncr=s keyDecr=S help='Rest density' ");

	TwAddVarRW(bar, "Restitution", TW_TYPE_FLOAT, &RESTITUTION,
		" label='Restitution' min=0 max=10 step=0.05 keyIncr=s keyDecr=S help='Determines the velocity 'lost' when a particle hits an object' ");

	TwAddVarRW(bar, "SurfaceTension", TW_TYPE_FLOAT, &SURFACETENSION,
		" label='Surface Tension' min=0.0001 max=2 step=0.0001 keyIncr=s keyDecr=S help='Determines the surface tension of the fluid' ");

	TwAddVarRW(bar, "SurfaceTensionThreshold", TW_TYPE_FLOAT, &THRESHOLD,
		" label='Surface Tension Threshold' min=1 max=100 step=0.005 keyIncr=s keyDecr=S help='Threshold where surface tension is calculated' ");
	
	TwAddVarRW(bar, "KernelSize", TW_TYPE_FLOAT, &H,
		" label='Kernel Size' min=0.0001 max=100 step=0.0001 keyIncr=s keyDecr=S help='Kernel size' ");

	TwAddVarRW(bar, "Stiff", TW_TYPE_FLOAT, &STIFF,
		" label='Stiff' min=1 max=100 step=1 keyIncr=s keyDecr=S help='Determines how close the particles of the simulation can be (1- close 5- they start pushing other particles)' ");

	TwAddVarRW(bar, "Mass", TW_TYPE_FLOAT, &MASS ,
		" label='Mass' min=0.01 max=100 step=0.01 keyIncr=s keyDecr=S help='Mass of the particle' ");
	

	// OpenGL settings 
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	setLighting();
	// enter GLUT's main loop
	glutMainLoop();

	return 1;
}