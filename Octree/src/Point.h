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
	double density;
	double pressure;
	
	glm::vec3 viscosity;
	//Cria um ponto na posi��o xyz, sem qualquer for�a aplicada ao mesmo
	Point(float x, float y, float z);
	~Point();
	void draw();
private:
};