#pragma once
#include <vector>

class Point {
	

public:
	std::vector<float> pos;
	std::vector<float> color;
	//Cria um ponto na posi��o xyz, sem qualquer for�a aplicada ao mesmo
	Point(float x, float y, float z);
	~Point();
	void draw();
private:
};