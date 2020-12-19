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

double PFORCE_FREQ = 8.0;
double PFORCE_MIN = 0;
//testar com 10 aqui 
double PFORCE_MAX = 0;
double H = 0.01; // Com o H a 10 ja apanha algumas particulas vizinhas
double m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow(H, 9));	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4
double m_SpikyKern = -45.0f / (3.141592 * pow(H, 6));			// Laplacian of viscocity (denominator): PI h^6
double m_LapKern = 45.0f / (3.141592 * pow(H, 6));

double GRAVITYVALUE = -9.8;
glm::vec3 GRAVITY = glm::vec3(0, GRAVITYVALUE, 0);
//mins e max para o volume da simula��o
double XMIN = 0, XMAX = 200;
double YMIN = 0, YMAX = 200;
double ZMIN = 0, ZMAX = 200;

double m_Time = 0;							// Start at T=0
//m_dt -> delta time
double m_DT = 0.003;

double PINTSTIFF = 1.5;
double PEXTSTIFF = 50000.0;
double PEXTDAMP = 100.0;
double PACCEL_LIMIT = 150.0;	// m / s^2
double PVEL_LIMIT = 3.0;		// m / s
double PMAX_FRAC = 1.0;

double EPSILON = 0.0001;
double RADIUS = 0.02;
double DECLIVE = 0.0;
double SIMULATIONSCALE = 0.005;
double MASS = 0.00020543;

double RESTDENSITY = 600;
double VISCOSITY = 0.35;
double K = 1.5;
//se o spacing for muito grande as particulas ignoram as colisoes com a caixa. nao sei porque ----- o max que posso ter � 3
//Spacing vai ser determinado pela rest density 
double SPACING = 3; //espaco usado para separar particulas

int points_p = 5000;
int pontos_inseridos = 0;
Octree* o;
double size = 100;
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

double alfa = 0.0f, beta = 0.5f, radius = 200.0f;
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


void SetupSpacing()
{
	
		// Determine spacing from density
		double PDIST =pow(MASS / RESTDENSITY, 1 / 3.0);
		SPACING = PDIST * 0.87 / SIMULATIONSCALE;

		H = SPACING; 

		//vai fazer update aos kernels

		m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow(H, 9));	// Wpoly6 kernel (denominator part) - 2003 Muller, p.4
		m_SpikyKern = -45.0f / (3.141592 * pow(H, 6));			// Laplacian of viscocity (denominator): PI h^6
		m_LapKern = 45.0f / (3.141592 * pow(H, 6));
		printf("Spacing -> %f \n", SPACING);
	
	
}

void restart() {
	points.clear();
	pontos_inseridos = 0;
	int countPoints = 0;
	std::chrono::steady_clock::time_point begin, end;
	bool done = false;
	double step = 1.0;
	double r1, r2, r3;
	Point* p;
	double sizex = XMAX - XMIN;
	if (sizex < 0)
		sizex = sizex * (-1);

	double sizez = YMAX - YMIN;
	if (sizez < 0)
		sizez = sizez * (-1);
	//controla o Y
	for (size_t i = 0; i < (o->size * 2) / SPACING && done == false; i++)
	{
		//controla o Z
		for (size_t j = 0; j < (sizez) / SPACING && done == false; j++)
		{
			//controla o X
			for (size_t k = 0; k < (sizex) / SPACING && done == false; k++)
			{
				if (countPoints >= points_p) {
					done = true;
					break;
				}
				//vai adicionar um bocado de random aos pontos criados
				double rx = ((double)rand() / (RAND_MAX));
				double rz = ((double)rand() / (RAND_MAX));
				r1 = (XMIN ) + k * SPACING;
				//o -10 � para nao estar a cirar as particulas mesmo em cima do cubo
				r2 = 15 + i * SPACING;
				r3 = (ZMIN)+j * SPACING;

				//isto é usado para nao gerar particulas perto das paredes
				r1 = r1 * 0.8 +5;
				r3 = r3 * 0.8+5;

				p = new Point(r1, r2, r3 );
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
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[�s]" << std::endl;
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
}

void processKeys(unsigned char c, int xx, int yy) {
	double step = 1.0;
	double r1, r2, r3;

	int countPoints = 0;
	bool done = false;
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

		//controla o Y
		for (size_t i = 0; i < (o->size * 2) / SPACING && done == false; i++)
		{
			//controla o Z
			for (size_t j = 0; j < (sizez) / SPACING && done == false; j++)
			{
				//controla o X
				for (size_t k = 0; k < (sizex) / SPACING && done == false; k++)
				{
					if (countPoints >= points_p) {
						done = true;
						break;
					}
					//vai adicionar um bocado de random aos pontos criados
					double rx = ((double)rand() / (RAND_MAX));
					double rz = ((double)rand() / (RAND_MAX));

					r1 = XMIN + k * SPACING;
					//o -10 � para nao estar a cirar as particulas mesmo em cima do cubo
					r2 = (15) + i * SPACING;
					r3 = ZMIN + j * SPACING;
					//isto é usado para nao gerar particulas perto das paredes
					r1 = r1 * 0.8 + 5;
					r3 = r3 * 0.8 + 5;
					p = new Point(r1+rx, r2, r3+rz);
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
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[�s]" << std::endl;
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;

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

	case GLUT_KEY_PAGE_DOWN: radius -= 2.0f;
		if (radius < 1.0f)
			radius = 1.0f;
		break;

	case GLUT_KEY_PAGE_UP: radius += 2.0f; break;
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



void kernelPoly6(double h, double* ret) {

	//*ret= (315 / (64 * 3.14159265359 * pow(h, 9))) * pow((pow(h, 2) - pow(glm::l1Norm(a - b), 2)), 3);
	*ret = 315.0f / (64.0f * 3.141592 * pow(h, 9));
	return;

}

void kernelSpiky(double h, double* ret) {

	*ret = ((double)-45.0 / (double)(3.14159265359 * pow(h, 6)));


	return;
}

void kernelSpikyPositive(double h, double* ret) {

	*ret = ((double)45.0 / (double)(3.14159265359 * pow(h, 6)));


	return;
}

void pressureGrid() {
	int i, j, cnt = 0;
	int nbr;
	float dx, dy, dz, sum, dsq, c;
	float d = SIMULATIONSCALE;
	float d2 = d * d;

	float radius = H / SIMULATIONSCALE;

	glm::vec3 dst;
	
	int nbrcnt = 0;
	int srch = 0;

	

	for each (Point * p in points ) {
		
		sum = 0.0;

		std::vector<Point*> query;
		int count = 0;
		o->queryOctree(p->pos, H, count, &query);

		for each (Point * p2 in query) {
			//se o ponto for o mesmo que o ponto atual, passa a frente
			if (p->pos.x == p2->pos.x && p->pos.y == p2->pos.y && p->pos.z == p2->pos.z)
				continue;
			dst = p2->pos;
			dst -= p->pos;
			dsq = d2 * (dst.x * dst.x + dst.y * dst.y + dst.z * dst.z);
			if (dsq <= (H*H)) {
				c = (H * H) - dsq;
				sum += c * c * c;
				
				
			}
			
			
		}
		//printf("Sum %f\n", sum);
		p->density = sum * MASS * m_Poly6Kern;
		p->pressure = (p->density - RESTDENSITY) * PINTSTIFF;
		
		p->density = 1.0f / p->density;

		//printf("Density %f\nPressure %f\n", p->density, p->pressure);

	}

}

void forceGrid() {
	
	register float pterm, vterm, dterm;
	int i, j, nbr;
	float c, d;
	float dx, dy, dz;
	float mR, mR2, visc;

	d = SIMULATIONSCALE;
	mR = H;
	visc = VISCOSITY;
	glm::vec3	jpos;
	float		jdist;
	float		jpress;
	float		jdensity;
	glm::vec3	jveleval;

	float		dsq;
	float		d2 = d * d;
	
	for each (Point * p in points) {
		p->force = glm::vec3(0);

		std::vector<Point*> query;
		int count = 0;
		o->queryOctree(p->pos, H, count, &query);
		//printf("2- A query apanhou %d elementos\n", count);
		for each (Point * p2 in query)
		{
			if (p->pos.x == p2->pos.x && p->pos.y == p2->pos.y && p->pos.z == p2->pos.z)
				continue;
			
			jpos = p2->pos;
			dx = (p->pos.x - jpos.x);		// dist in cm
			dy = (p->pos.y - jpos.y);
			dz = (p->pos.z - jpos.z);
			dsq = d2 * (dx * dx + dy * dy + dz * dz);
			if (dsq <= (H*H)) {

				jdist = sqrt(dsq);

				jpress = p2->pressure;
				jdensity = p2->density;
				jveleval = p2->velocityEval;
				dx = (p->pos.x - jpos.x);		// dist in cm
				dy = (p->pos.y - jpos.y);
				dz = (p->pos.z - jpos.z);
				c = (mR - jdist);
				pterm = d * -0.5f * c * m_SpikyKern * (p->pressure + jpress) / jdist;
				dterm = c * (p->density) * jdensity;
				vterm = m_LapKern * visc;
				p->force.x += (pterm * dx + vterm * (jveleval.x - p->velocityEval.x)) * dterm;
				p->force.y += (pterm * dy + vterm * (jveleval.y - p->velocityEval.y)) * dterm;
				p->force.z += (pterm * dz + vterm * (jveleval.z - p->velocityEval.z)) * dterm;

				//printf("Force %f %f %f\n", p->force.x, p->force.y, p->force.z);
			}
			
		}
		//printf("------\n");
	}
}

void move() {
	glm::vec3 norm, z;
	glm::vec3 dir, accel;
	glm::vec3 vnext;
	glm::vec3 bmin, bmax;
	glm::vec3 clr;
	double adj;
	float AL, AL2, SL, SL2, ss, radius;
	float stiff, damp, speed, diff;

	AL = PACCEL_LIMIT;	AL2 = AL * AL;
	SL = PVEL_LIMIT;	SL2 = SL * SL;

	stiff = PEXTSTIFF;
	damp = PEXTDAMP;
	radius = RADIUS;
	bmin = glm::vec3(XMIN, YMIN, ZMIN);
	bmax = glm::vec3(XMAX, YMAX, ZMAX);;
	ss = SIMULATIONSCALE;

	for each (Point * p in points) {
		// Compute Acceleration		
		accel = p->force;
		accel *= MASS;

		// Boundary Conditions
		// Y-axis walls
		diff = radius - (p->pos.y - (bmin.y + (p->pos.x - bmin.x) * DECLIVE)) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(-DECLIVE, 1 - DECLIVE, 0);

			adj = stiff * diff - damp * glm::dot(norm, p->velocityEval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}
		diff = radius - (bmax.y - p->pos.y) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(0, -1, 0);
			adj = stiff * diff - damp * glm::dot(norm, p->velocityEval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}

		// X-axis walls
		if (true) {
			diff = radius - (p->pos.x - (bmin.x + (sin(m_Time * PFORCE_FREQ) + 1) * 0.5 * PFORCE_MIN)) * ss;
			//diff = 2 * radius - ( p->pos.x - min.x + (sin(m_Time*10.0)-1) * m_Param[FORCE_XMIN_SIN] )*ss;	
			if (diff > EPSILON) {
				norm = glm::vec3(1.0, 0, 0);
				adj = (PFORCE_MIN + 1) * stiff * diff - damp * glm::dot(norm, p->velocityEval);
				accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
			}

			diff = radius - ((bmax.x - (sin(m_Time * PFORCE_FREQ) + 1) * 0.5 * PFORCE_MAX) - p->pos.x) * ss;
			if (diff > EPSILON) {
				norm = glm::vec3(-1, 0, 0);
				adj = (PFORCE_MAX + 1) * stiff * diff - damp * glm::dot(norm, p->velocityEval);
				accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
			}
		}

		// Z-axis walls
		diff = radius - (p->pos.z - bmin.z) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(0, 0, 1);
			adj = stiff * diff - damp * glm::dot(norm, p->velocityEval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}
		diff = radius - (bmax.z - p->pos.z) * ss;
		if (diff > EPSILON) {
			norm = glm::vec3(0, 0, -1);
			adj = stiff * diff - damp * glm::dot(norm, p->velocityEval);
			accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
		}




		// Plane gravity
		accel += GRAVITY;



		// Acceleration limiting 
		speed = accel.x * accel.x + accel.y * accel.y + accel.z * accel.z;
		if (speed > AL2) {
			accel *= AL / sqrt(speed);
		}

		// Velocity limiting 
		speed = p->velocity.x * p->velocity.x + p->velocity.y * p->velocity.y + p->velocity.z * p->velocity.z;
		if (speed > SL2) {
			speed = SL2;
			p->velocity *= SL / sqrt(speed);
		}

		// Leapfrog Integration ----------------------------
		vnext = accel;
		vnext *= m_DT;
		vnext += p->velocity;						// v(t+1/2) = v(t-1/2) + a(t) dt

		p->velocityEval = p->velocity;
		p->velocityEval += vnext;
		p->velocityEval *= 0.5;					// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
		p->velocity = vnext;
		vnext *= m_DT / ss;
		p->pos += vnext;						// p(t+1) = p(t) + v(t+1/2) dt

		/*if ( m_Param[PCLR_MODE]==1.0 ) {
			adj = fabs(vnext.x)+fabs(vnext.y)+fabs(vnext.z) / 7000.0;
			adj = (adj > 1.0) ? 1.0 : adj;
			*pclr = COLORA( 0, adj, adj, 1 );
		}
		if ( m_Param[PCLR_MODE]==2.0 ) {
			float v = 0.5 + ( *ppress / 1500.0);
			if ( v < 0.1 ) v = 0.1;
			if ( v > 1.0 ) v = 1.0;
			*pclr = COLORA ( v, 1-v, 0, 1 );
		}*/


		// Euler integration -------------------------------
		/* accel += m_Gravity;
		accel *= m_DT;
		p->vel += accel;				// v(t+1) = v(t) + a(t) dt
		p->vel_eval += accel;
		p->vel_eval *= m_DT/d;
		p->pos += p->vel_eval;
		p->vel_eval = p->vel;  */



	}

}



void renderScene(void) {
	char fpss[200];

	static double t = 0;



	static double size = 1.5;
	static bool grow = true;
	static double speed = 0.01;
	static double maxsize = 2;
	static double minsize = 0.1;
	static double g_RotateS = 0.1;
	static double g_RotateX = 0.0;

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


	//Draw stuff
	//drawAxis();

	o->drawOut();

	//Vai afetar o vector force de cada particula


	//Vai calcular a density de cada particula
	pressureGrid();

	forceGrid();

	//usado para mover os pontos
	move();

	//Vai reconstruir a octree
	updatePoints();

	int count = 0;
	//count � usado para contar quantos pontos sao "checados" a ver se isto � eficiente
	//o->queryOctree(center, size_query,count);
	drawQueryVolume();


	//glColor3f(1.0f, 1.0f, 1.0f);
	glPointSize(3);
	glBegin(GL_POINTS);
	for (size_t i = 0; i < points.size(); i++)
	{
		points.at(i)->draw();
	}
	glVertex3f(XMIN, YMIN, ZMIN);
	glVertex3f(XMAX, YMIN, ZMIN);
	glVertex3f(XMIN, YMAX, ZMIN);
	glVertex3f(XMIN, YMIN, ZMAX);

	glVertex3f(XMAX, YMAX, ZMAX);
	glVertex3f(XMIN, YMAX, ZMAX);
	glVertex3f(XMAX, YMIN, ZMAX);
	glVertex3f(XMAX, YMAX, ZMIN);

	glEnd();

	glColor3f(0.0f, 0.0f, 1.0f);

	TwDraw();

	// End of frame
	glutSwapBuffers();
}

void SetupDualWave() {
	XMIN = -50;
	XMAX = 50;
	YMIN = 0;
	YMAX = 100;
	ZMIN = -50;
	ZMAX = 50;

	PFORCE_MIN= 20.0;
	PFORCE_MAX= 20.0;
	GRAVITY = glm::vec3(0.0f, -9.8f, 0.0f);

}

void SetupNormal() {
	XMIN = -50;
	XMAX = 50;
	YMIN = 0;
	YMAX = 100;
	ZMIN = -50;
	ZMAX = 50;

	GRAVITY = glm::vec3(0.0f, -9.8f, 0.0f);

}

int main(int argc, char** argv) {
	o = new Octree(0, 100, 0, size, max_points, max_depth);

	SetupDualWave();
	//SetupNormal();

	SetupSpacing();

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

	// Add 'speed' to 'bar': it is a modifable (RW) variable of type TW_TYPE_DOUBLE. Its key shortcuts are [s] and [S].
	TwAddVarRW(bar, "Delta Time", TW_TYPE_DOUBLE, &m_DT,
		" label='Delta Time' min=0.001 max=2 step=0.001 keyIncr=s keyDecr=S help='Delta Time' ");

	// Add 'wire' to 'bar': it is a modifable variable of type TW_TYPE_BOOL32 (32 bits boolean). Its key shortcut is [w].
	TwAddVarRW(bar, "Particles", TW_TYPE_INT32, &points_p,
		" label='Particles' min=1 max=10000000 step=1000 key=w help='Number of particles.' ");
	//Se aumentar o H para 40 da para ver umas coisas engra�adas, talvez de para perceber o que esta a acontecer de errado
	TwAddVarRW(bar, "Kernel H", TW_TYPE_DOUBLE, &H,
		" label='Kernel H' min=0.1 max=50 step=0.5 keyIncr=s keyDecr=S help='Kernel H size' ");

	TwAddVarRW(bar, "Declive", TW_TYPE_DOUBLE, &DECLIVE,
		" label='Declive' min=0.0 max=1 step=0.05 keyIncr=s keyDecr=S help='Declive' ");

	//Apenas controla as intera��es. Desenha sempre com radius a 1
	TwAddVarRW(bar, "Particle Radius", TW_TYPE_DOUBLE, &RADIUS,
		" label='Particle Radius' min=0.01 max=50 step=0.01 keyIncr=s keyDecr=S help='Particle Radius - Only controlls the interactions, the draw will have a fixed radius' ");

	TwAddVarRW(bar, "Particle Spacing", TW_TYPE_DOUBLE, &SPACING,
		" label='Particle Spacing' min=1 max=10 step=0.5 keyIncr=s keyDecr=S help='Spacing between particles' ");

	TwAddVarRW(bar, "Simulation Scale", TW_TYPE_DOUBLE, &SIMULATIONSCALE,
		" label='Simulation Scale' min=0.001 max=3 step=0.001 keyIncr=s keyDecr=S help='Simulation Scale' ");

	TwAddVarRW(bar, "Gravity", TW_TYPE_DOUBLE, &GRAVITYVALUE,
		" label='Gravity' min=-200 max=200 step=1 keyIncr=s keyDecr=S help='Gravity' ");
	TwAddVarRW(bar, "Damp", TW_TYPE_DOUBLE, &PEXTDAMP,
		" label='Damp' min=1 max=1000 step=10 keyIncr=s keyDecr=S help='Atenution when particle hits wall' ");

	TwAddVarRW(bar, "Viscosity", TW_TYPE_DOUBLE, &VISCOSITY,
		" label='Viscosity' min=0.1 max=1 step=0.1 keyIncr=s keyDecr=S help='Fluid Viscosity' ");
	

	// OpenGL settings 
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	// enter GLUT's main loop
	glutMainLoop();

	return 1;
}