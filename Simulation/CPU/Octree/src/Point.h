#pragma once

#include "../glm/vec3.hpp"

class Point {
	

public:
	int id;
	float originalPos[3];
	float pos[3];
	float color[3];
	float force[3];
	float velocity[3];
	float velocityEval[3];
	float acceleration[3];
	float gravity[3];
	float surfaceNormal[3];
	float surfaceTension[3];
	
	float density;
	float pressure;
	
	float viscosity[3];

	unsigned int mortonCell;

	

	//Cria um ponto na posição xyz, sem qualquer força aplicada ao mesmo
	Point(float x, float y, float z,int id,unsigned int morton_cell);
	~Point();
	void draw(float r);
private:
};