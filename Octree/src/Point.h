#pragma once
#include <vector>

class Point {
	

public:
	std::vector<float> pos;
	
	//Cria um ponto na posi��o xyz, sem qualquer for�a aplicada ao mesmo
	Point(float x, float y, float z);
	void draw();
private:
};