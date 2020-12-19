#pragma once

#include "../glm/vec3.hpp"

class Point {
	

public:
	glm::vec3 originalPos;
	glm::vec3 pos;
	glm::vec3 color;
	glm::vec3 force;
	glm::vec3 velocity;
	glm::vec3 velocityEval;
	glm::vec3 acceleration;
	double density;
	double pressure;
	
	glm::vec3 viscosity;
	//Cria um ponto na posição xyz, sem qualquer força aplicada ao mesmo
	Point(float x, float y, float z);
	~Point();
	void draw();
private:
};