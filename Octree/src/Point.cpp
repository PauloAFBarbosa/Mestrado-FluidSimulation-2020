#include "Point.h"

#include "../toolkits/glut/GL/glut.h"

Point::Point(float x, float y, float z,int id)
{
	this->id = id;
	this->originalPos[0] = x;
	this->originalPos[1] = y;
	this->originalPos[2] = z;
	this->pos[0] = x;
	this->pos[1] = y;
	this->pos[2] = z;
	
	this->color[0] = 1;
	this->color[1] = 1;
	this->color[2] = 1;

	this->force[0] = 0;
	this->force[1] = 0;
	this->force[2] = 0;

	
	this->velocity[0] = 0;
	this->velocity[1] = 0;
	this->velocity[2] = 0;
	this->velocityEval[0] = 0;
	this->velocityEval[1] = 0;
	this->velocityEval[2] = 0;
	this->force[0] = 0;
	this->force[1] = 0;
	this->force[2] = 0;
	this->viscosity[0] = 0;
	this->viscosity[1] = 0;
	this->viscosity[2] = 0;
	this->acceleration[0] = 0;
	this->acceleration[1] = 0;
	this->acceleration[2] = 0;
	this->gravity[0] = 0;
	this->gravity[1] = 0;
	this->gravity[2] = 0;
	this->surfaceNormal[0] = 0;
	this->surfaceNormal[1] = 0;
	this->surfaceNormal[2] = 0;
	this->surfaceTension[0] = 0;
	this->surfaceTension[1] = 0;
	this->surfaceTension[2] = 0;

	this->density = 0;
	this->pressure = 0;
	
}

Point::~Point()
{
	
}

void Point::draw(double r)
{
	
	glColor3f(this->color[0], this->color[1], this->color[2]);
	glPushMatrix();
	glTranslatef(this->pos[0], this->pos[1], this->pos[2]);
	glutSolidSphere(r, 5, 5);
	
	glPopMatrix();
	//glVertex3f(this->pos.x, this->pos.y, this->pos.z);

	this->color[0] = 1;
	this->color[1] = 1;
	this->color[2] = 1;

}
