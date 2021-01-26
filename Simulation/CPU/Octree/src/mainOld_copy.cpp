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

double GRAVITYVALUE = -9.8;
glm::vec3 GRAVITY= glm::vec3(0, GRAVITYVALUE, 0);
//mins e max para o volume da simulação
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
double RADIUS = 1;
double DECLIVE = 0.0;
double SIMULATIONSCALE = 0.002;
double MASS = 0.00020543;

double H =3; // Com o H a 10 ja apanha algumas particulas vizinhas
double RESTDENSITY= 1500;
double VISCOSITY=0.35;
double K = 1.5;
//se o spacing for muito grande as particulas ignoram as colisoes com a caixa. nao sei porque ----- o max que posso ter é 3
double SPACING = 3; //espaco usado para separar particulas

int points_p = 5000;
int pontos_inseridos = 0;
Octree* o;
double size = 100;
int max_points = 100;
int max_depth = 5;

int frame=0;
double fps=0;
glm::vec3 center= glm::vec3( 0, 0, 0 );
double size_query=H*2;
std::vector<Point *> points;

int mytime;
int timebase;
int maxTime = 10000;
double stepT = 100.0f / maxTime;

double alfa = 0.0f, beta = 0.5f, radius = 500.0f;
double camX = 50, camY = 50, camZ = -50;
double rot = 0;
int winid = 0;

//Da maneira que fiz agora ele verifica 1734 pontos. Em vez de 10 000 pontos.



/**
 * @brief calculates the cam values (used in GluLookAt function) from the alteration of alfa, beta and radius
 * 
 */
void spherical2Cartesian() {
	camX = radius * cos(beta) * sin(alfa) + 100;
	camY = radius * sin(beta) +100;
	camZ = radius * cos(beta) * cos(alfa)+100;
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

	//controla o Y
	for (size_t i = 0; i < (o->size * 2) / SPACING && done == false; i++)
	{
		//controla o Z
		for (size_t j = 0; j < (o->size * 2) / SPACING && done == false; j++)
		{
			//controla o X
			for (size_t k = 0; k < (o->size * 2) / SPACING && done == false; k++)
			{
				if (countPoints >= points_p) {
					done = true;
					break;
				}
				r1 = k * SPACING;
				//o -10 é para nao estar a cirar as particulas mesmo em cima do cubo
				r2 = (size * 2 - 10) - i * SPACING;
				r3 = j * SPACING;

				p = new Point(r1, r2, r3);
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
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
}

void processKeys(unsigned char c, int xx, int yy) {
	double step = 1.0;
	double r1,r2,r3;
	
	int countPoints = 0;
	bool done = false;
	std::chrono::steady_clock::time_point begin, end;
	Point* p;
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
			for (size_t i = 0; i < (o->size*2)/SPACING && done ==false ; i++)
			{
				//controla o Z
				for (size_t j = 0; j < (o->size * 2) / SPACING && done == false; j++)
				{
					//controla o X
					for (size_t k = 0; k < (o->size * 2) / SPACING && done == false; k++)
					{
						if (countPoints >= points_p) {
							done = true;
							break;
						}
						r1 = k * SPACING;
						//o -10 é para nao estar a cirar as particulas mesmo em cima do cubo
						r2 = (size*2 -10) - i * SPACING;
						r3 = j * SPACING;

						p = new Point(r1, r2, r3);
						points.push_back(p);
						countPoints++;
					}

				}
			}

		//separei a parte de inserir os pontos porque quero medir apenas a parte de inserir.
		//no programa ele vai gerar os pontos 1 vez e vai inserir numa nova arvore cada frame
		begin = std::chrono::steady_clock::now();
		for (int i = 0; i < points_p; i++) {
			o->insertPoint(points.at(pontos_inseridos+i));
		}

		end = std::chrono::steady_clock::now();
		pontos_inseridos += points_p;
		std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
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

void move() {
	
	double dif;
	for each (Point * p in points)
	{
		glm::vec3 acceleration = p->force;
		acceleration *= MASS;
		
		//vai ver se o ponto saiu para baixo
		//estao a cair do cubo.
		//vai calcular a colisão com o centro da particula --> se dif for positivo esta a baixo da base
		double slope_x = p->pos.x - XMIN;
		dif = RADIUS - (p->pos.y - (YMIN + (slope_x)*DECLIVE)) * SIMULATIONSCALE;

		if (dif > EPSILON) {
			//vai fazer a normal com slope, caso exista
			glm::vec3 normal = glm::vec3(-DECLIVE,1- DECLIVE,0);

			double adj = PINTSTIFF * dif - PEXTDAMP * glm::dot(normal, p->velocityEval);
			
			acceleration.x += adj * normal.x ;
			acceleration.y += adj * normal.y ;
			acceleration.z += adj * normal.z ;
		}

		
		// Nao vai fazer bound na parede de cima
		dif = RADIUS - (YMAX - p->pos.y) * SIMULATIONSCALE;
		if (dif > EPSILON) {
			glm::vec3 normal = glm::vec3(0, -1, 0);
			double adj = PINTSTIFF * dif - PEXTDAMP * glm::dot(normal, p->velocityEval);
			acceleration.x += adj * normal.x ;
			acceleration.y += adj * normal.y ;
			acceleration.z += adj * normal.z ;
			
		}

		//X-azis walls
		dif = RADIUS - (p->pos.x - XMIN) * SIMULATIONSCALE;
		if (dif > EPSILON) {
			glm::vec3 normal = glm::vec3(1, 0, 0);
			double adj = PINTSTIFF * dif - PEXTDAMP * glm::dot(normal, p->velocityEval);
			acceleration.x += adj * normal.x;
			acceleration.y += adj * normal.y;
			acceleration.z += adj * normal.z;
		}
		dif = RADIUS - (XMAX - p->pos.x) * SIMULATIONSCALE;
		if (dif > EPSILON) {
			glm::vec3 normal = glm::vec3(-1, 0, 0);
			double adj = PINTSTIFF * dif - PEXTDAMP * glm::dot(normal, p->velocityEval);
			acceleration.x += adj * normal.x;
			acceleration.y += adj * normal.y;
			acceleration.z += adj * normal.z;
		}

		// Z-axis walls
		dif = RADIUS - (p->pos.z - ZMIN) * SIMULATIONSCALE;
		if (dif > EPSILON) {
			glm::vec3 normal = glm::vec3(0, 0, 1);
			double adj = PINTSTIFF * dif - PEXTDAMP * glm::dot(normal, p->velocityEval);
			acceleration.x += adj * normal.x;
			acceleration.y += adj * normal.y;
			acceleration.z += adj * normal.z;
		}
		dif = RADIUS - (ZMAX - p->pos.z) * SIMULATIONSCALE;
		if (dif > EPSILON) {
			glm::vec3 normal = glm::vec3(0, 0, -1);
			double adj = PINTSTIFF * dif - PEXTDAMP * glm::dot(normal, p->velocityEval);
			acceleration.x += adj * normal.x;
			acceleration.y += adj * normal.y;
			acceleration.z += adj * normal.z;
		}
		
		// Plane gravity
		acceleration += GRAVITY;

		// Acceleration limiting 
		double speed = acceleration.x * acceleration.x + acceleration.y * acceleration.y + acceleration.z * acceleration.z;
		if (speed > PACCEL_LIMIT* PACCEL_LIMIT) {
			acceleration *= PACCEL_LIMIT / sqrt(speed);
		}

		// Velocity limiting 
		speed = p->velocity.x * p->velocity.x + p->velocity.y * p->velocity.y + p->velocity.z * p->velocity.z;
		if (speed > PVEL_LIMIT* PVEL_LIMIT) {
			speed = PVEL_LIMIT;
			p->velocity *= PVEL_LIMIT / sqrt(speed);
		}
		
		// Leapfrog Integration ----------------------------
		glm::vec3 vnext = acceleration;
		vnext *= m_DT;
		vnext += p->velocity;						// v(t+1/2) = v(t-1/2) + a(t) dt

		p->velocityEval = p->velocity;
		p->velocityEval += vnext;
		p->velocityEval *= 0.5;					// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
		p->velocity = vnext;
		vnext *= m_DT / SIMULATIONSCALE;
		p->pos += vnext;						// p(t+1) = p(t) + v(t+1/2) dt

		
	}
}

void kernelPoly6(double h,double * ret) {
	
	//*ret= (315 / (64 * 3.14159265359 * pow(h, 9))) * pow((pow(h, 2) - pow(glm::l1Norm(a - b), 2)), 3);
	*ret = 315.0f / (64.0f * 3.141592 * pow(h, 9));
	return;
	
}

 void kernelSpiky(double h, double * ret) {

	*ret = ((double)-45.0 / (double) (3.14159265359 * pow(h, 6))) ;
	
	
	return ;
}

 void kernelSpikyPositive(double h, double * ret) {

	 *ret = ((double)45.0 / (double)(3.14159265359 * pow(h, 6)));


	 return;
 }

void pressureGrid() {
	double sum;
	//calcula todas as densidades
	for each (Point * p in points)
	{
		sum = 0;
		//Vai buscar vizinhos
		std::vector<Point*> query;
		int count = 0;
		o->queryOctree(p->pos, H, count,&query);
		
		//printf("Apanhou %d na query\n", query.size());
		//Density
		p->density = 0;
		p->pressure = 0;
		for each (Point * p2 in query)
		{
			//printf("---------- Entrou no for da query\n");
			//Se for a mesma particula passa a frante
			//printf("---------- Paticula 1 pos %f %f %f\n",p->pos.x, p->pos.y, p->pos.z);
			//printf("---------- Paticula 2 pos %f %f %f\n", p2->pos.x, p2->pos.y, p2->pos.z);
			if (p->pos.x == p2->pos.x && p->pos.y == p2->pos.y && p->pos.z == p2->pos.z)
				continue;
			//printf("---------- Passou o if que verifica se é a mesma particula\n");
			//vai calcular a distancia entre os dois pontos
			double d = SIMULATIONSCALE * SIMULATIONSCALE *(p2->pos.x * p2->pos.x + p2->pos.y * p2->pos.y + p2->pos.z * p2->pos.z);

			//printf("---------- D-> %f\n", d);

			if (d <= (H * H)) {
				//printf("\n\n-------------entrou aqui-------------- \n\n");
				double dif = (H * H) - d;
				sum += dif * dif * dif; //^3
			}
		}
		double kernel;
		
		kernelPoly6(H ,&kernel );
		//density
		
		p->density = sum * MASS * kernel;
		p->pressure = (p->density - RESTDENSITY) * K;
		p->density = 1 / p->density;
		//printf("Kernel %lf \nSUM %f\nParticula Density %f \n",kernel,sum,p->density);

		//printf("Particula pressure %f \n", p->pressure);
		
		
		//printf("Density %f \n", p->density);
	}
	
}

void forceGrid() {

	
	//calcula todas as densidades
	for each (Point * p in points)
	{
		p->force = glm::vec3(0);
		//Vai buscar vizinhos
		std::vector<Point*> query;
		int count = 0;
		o->queryOctree(p->pos, H, count,&query);
		//printf("2- A query apanhou %d elementos\n", count);
		for each (Point * p2 in query)
		{
			//Se for a mesma particula passa a frante
			if (p->pos.x == p2->pos.x && p->pos.y == p2->pos.y && p->pos.z == p2->pos.z)
				continue;
			//calcular distancia
			double dx, dy, dz,dist;
			dx = p->pos.x - p2->pos.x;
			dy = p->pos.y - p2->pos.y;
			dz = p->pos.z - p2->pos.z;

			dist = SIMULATIONSCALE * SIMULATIONSCALE * (dx * dx + dy * dy + dz * dz);
			//printf("Dist %f\n", dist);
			if (dist <= (H * H)) {
				double sqrdist = sqrt(dist);

				float p2Pressure = p2->pressure;
				float p2Density = p2->density;
				glm::vec3 p2VelEval = p2->velocityEval;

				//printf("sqrDist %f\n", sqrdist);
				dx = p->pos.x - p2->pos.x;
				dy = p->pos.y - p2->pos.y;
				dz = p->pos.z - p2->pos.z;
				
				double dif = H - sqrdist;
				//printf("dif %f\n", dif);
				double kernel;
				
				kernelSpiky(H, &kernel);
				//posso usar p2Pressure ou p2->pressure
				double pressureTerm = SIMULATIONSCALE * (-0.5f)* dif * kernel * (p->pressure + p2Pressure) / sqrdist;
				double densityTerm = dif * p->density * p2->density;
				
				kernelSpikyPositive(H, &kernel);
				double viscosityTerm = kernel * VISCOSITY;

				

				p->force.x += (pressureTerm * dx + viscosityTerm * (p2->velocityEval.x - p->velocityEval.x)) * densityTerm;
				p->force.y += (pressureTerm * dy + viscosityTerm * (p2->velocityEval.y - p->velocityEval.y)) * densityTerm;
				p->force.z += (pressureTerm * dz + viscosityTerm * (p2->velocityEval.z - p->velocityEval.z)) * densityTerm;
				//printf("Pressure term %f\n", pressureTerm);
				//printf("Visc term %f\n", viscosityTerm);
				//printf("densityTerm %f\n", densityTerm);
				//printf("Forca da particula %f %f %f\n", p->force.x, p->force.y, p->force.z);
			}
		}

		//printf("Density %f \n", p->density);
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
		100, 100, 100,
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

	o->draw();

	//Vai afetar o vector force de cada particula
	

	//Vai calcular a density de cada particula
	pressureGrid();

	forceGrid();

	//usado para mover os pontos
	move();

	//Vai reconstruir a octree
	updatePoints();
	
	int count=0;
	//count é usado para contar quantos pontos sao "checados" a ver se isto é eficiente
	//o->queryOctree(center, size_query,count);
	drawQueryVolume();


	//glColor3f(1.0f, 1.0f, 1.0f);
	glPointSize(3);
	glBegin(GL_POINTS);
	for (size_t i = 0; i < points.size(); i++)
	{
		points.at(i)->draw();
	}
	glVertex3f(0, 0, 0);
	glVertex3f(200, 0, 0);
	glVertex3f(200, 0, 200);
	glVertex3f(0, 0, 200);

	glVertex3f(0, 200, 0);
	glVertex3f(200, 200, 0);
	glVertex3f(200, 200, 200);
	glVertex3f(0, 200, 200);
	glEnd();

	glColor3f(0.0f, 0.0f, 1.0f);

	TwDraw();

	// End of frame
	glutSwapBuffers();
}

int main(int argc, char** argv) {
	o= new Octree(100,100,100,size,max_points,max_depth);

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
	//Se aumentar o H para 40 da para ver umas coisas engraçadas, talvez de para perceber o que esta a acontecer de errado
	TwAddVarRW(bar, "Kernel H", TW_TYPE_DOUBLE, &H,
		" label='Kernel H' min=0.1 max=50 step=0.5 keyIncr=s keyDecr=S help='Kernel H size' ");

	TwAddVarRW(bar, "Declive", TW_TYPE_DOUBLE, &DECLIVE,
		" label='Declive' min=0.0 max=1 step=0.05 keyIncr=s keyDecr=S help='Declive' ");

	//Apenas controla as interações. Desenha sempre com radius a 1
	TwAddVarRW(bar, "Particle Radius", TW_TYPE_DOUBLE, &RADIUS,
		" label='Particle Radius' min=0.01 max=50 step=0.5 keyIncr=s keyDecr=S help='Particle Radius - Only controlls the interactions, the draw will have a fixed radius' ");

	TwAddVarRW(bar, "Particle Spacing", TW_TYPE_DOUBLE, &SPACING,
		" label='Particle Spacing' min=1 max=10 step=0.5 keyIncr=s keyDecr=S help='Spacing between particles' ");

	TwAddVarRW(bar, "Simulation Scale", TW_TYPE_DOUBLE, &SIMULATIONSCALE,
		" label='Simulation Scale' min=0.001 max=3 step=0.001 keyIncr=s keyDecr=S help='Simulation Scale' ");

	TwAddVarRW(bar, "Gravity", TW_TYPE_DOUBLE, &GRAVITYVALUE,
		" label='Gravity' min=-200 max=200 step=1 keyIncr=s keyDecr=S help='Gravity' ");
	
	// OpenGL settings 
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	// enter GLUT's main loop
	glutMainLoop();

	return 1;
}