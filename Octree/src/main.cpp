#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include "Octree.h"
#include<GL/glut.h>

#endif

#define _USE_MATH_DEFINES
#include <math.h>

std::vector<Point *> points;
Octree * o;
float size = 10;
int max_points = 1;
int max_depth = 2;

int time;
int timebase;
int maxTime = 10000;
float stepT = 100.0f / maxTime;

float alfa = 0.0f, beta = 0.5f, radius = 20.0f;
float camX = 0, camY = 5, camZ = 20;
float rot = 0;
int winid = 0;

/**
 * @brief calculates the cam values (used in GluLookAt function) from the alteration of alfa, beta and radius
 * 
 */
void spherical2Cartesian() {
	camX = radius * cos(beta) * sin(alfa);
	camY = radius * sin(beta);
	camZ = radius * cos(beta) * cos(alfa);
}

void processKeys(unsigned char c, int xx, int yy) {
	float step = 0.1;
	float r1,r2,r3;
	Point* p;
	switch (c)
	{
	case 'p':

		r1 = (float) (rand()) / ((float) (RAND_MAX ));
		r2 = (float)(rand()) / ((float)(RAND_MAX ));
		r3 = (float)(rand()) / ((float)(RAND_MAX ));
		
		r1=r1* (size * 2) - size;
		r2 = r2 * (size * 2) - size;
		r3 = r3 * (size * 2) - size;
		p = new Point(r1, r2, r3);

		

		o->insertPoint(p);
		points.push_back(p);
		
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


void renderScene(void) {
	static float t = 0;
	float pos[3], der[3];
	float m[4][4];

	static float size = 1.5;
	static bool grow = true;
	static float speed = 0.01;
	static float maxsize = 2;
	static float minsize = 0.1;
	static float g_RotateS = 0.1;
	static float g_RotateX = 0.0;

	// clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// set the camera
	glLoadIdentity();
	gluLookAt(camX, camY, camZ,
		0.0, 0.0, 0.0,
		0.0f, 1.0f, 0.0f);
	//Draw stuff
	//drawAxis();

	o->draw();
	glColor3f(1.0f, 1.0f, 1.0f);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (size_t i = 0; i < points.size(); i++)
	{
		points.at(i)->draw();
	}
	glEnd();

	glColor3f(0.0f, 0.0f, 1.0f);

	// End of frame
	glutSwapBuffers();
}

int main(int argc, char** argv) {
	o= new Octree(0,0,0,size,max_points,max_depth);

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

	// enter GLUT's main loop
	glutMainLoop();

	return 1;
}