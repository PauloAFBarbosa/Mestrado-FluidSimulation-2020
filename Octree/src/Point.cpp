#include "Point.h"

#include "../toolkits/glut/GL/glut.h"

Point::Point(float x, float y, float z)
{
	this->originalPos = glm::vec3(x,y,z);
	this->pos = glm::vec3(x, y, z);
	
	this->color = glm::vec3(1, 1, 1);

	this->force = glm::vec3(0, 0, 0);

	
	this->velocity = glm::vec3 (0);
	this->velocityEval = glm::vec3(0);
	this->density = 0;
	this->pressure= 0;
	this->force = glm::vec3(0);
	this->viscosity = glm::vec3(0);
	this->acceleration = glm::vec3(0);
	this->gravity = glm::vec3(0);
	this->surfaceNormal = glm::vec3(0);
	this->surfaceTension = glm::vec3(0);
	
}

Point::~Point()
{
	
}

void Point::draw(double r)
{
	
	glColor3f(this->color.x, this->color.y, this->color.z);
	glPushMatrix();
	glTranslatef(this->pos.x, this->pos.y, this->pos.z);
	glutSolidSphere(r, 4, 4);
	
	glPopMatrix();
	//glVertex3f(this->pos.x, this->pos.y, this->pos.z);

	this->color.x = 1;
	this->color.y = 1;
	this->color.z = 1;

}
