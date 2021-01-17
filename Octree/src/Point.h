#pragma once

#include "../glm/vec3.hpp"

class Point {
	

public:
	int id;
	double originalPos[3];
	double pos[3];
	double color[3];
	double force[3];
	double velocity[3];
	double velocityEval[3];
	double acceleration[3];
	double gravity[3];
	double surfaceNormal[3];
	double surfaceTension[3];
	
	double density;
	double pressure;
	
	double viscosity[3];

	

	//Cria um ponto na posição xyz, sem qualquer força aplicada ao mesmo
	Point(float x, float y, float z,int id);
	~Point();
	void draw(double r);
private:
};