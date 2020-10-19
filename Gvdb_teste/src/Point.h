#pragma once
#include <vector>

class Point {
	

public:
	std::vector<float> pos;
	std::vector<float> color;
	//Cria um ponto na posição xyz, sem qualquer força aplicada ao mesmo
	Point(float x, float y, float z);
	~Point();
	void draw();
private:
};