#include "Point.h"
#include <vector>
#include "../toolkits/glut/GL/glut.h"

Point::Point(float x, float y, float z)
{
	std::vector<float> temp (3);
	temp.at(0) = x;
	temp.at(1) = y;
	temp.at(2) = z;
	this->pos = temp;

	std::vector<float> temp2(3);
	temp2.at(0) = 1;
	temp2.at(1) = 1;
	temp2.at(2) = 1;

	this->color = temp2;
	
}

Point::~Point()
{
	
}

void Point::draw()
{
	glColor3f(this->color.at(0), this->color.at(1), this->color.at(2));
	glVertex3f(this->pos.at(0), this->pos.at(1), this->pos.at(2));

	this->color.at(0) = 1;
	this->color.at(1) = 1;
	this->color.at(2) = 1;

}
