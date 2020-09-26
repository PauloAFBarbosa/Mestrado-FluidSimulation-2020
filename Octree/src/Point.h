#pragma once
#include <vector>

class Point {
	

public:
	std::vector<float> pos;
	
	//Cria um ponto na posição xyz, sem qualquer força aplicada ao mesmo
	Point(float x, float y, float z);
	void draw();
private:
};