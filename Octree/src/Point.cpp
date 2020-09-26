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
	
}

void Point::draw()
{
	glVertex3f(this->pos.at(0), this->pos.at(1), this->pos.at(2));
}
