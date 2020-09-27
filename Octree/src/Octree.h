#pragma once
#include <vector>
#include "../src/Point.h"

class Octree {


	public:
		std::vector<float> pos;
		//metade do comprimento do lado do cubo
		float size;
		bool devided;
		int max_points;
		int max_depth;
		std::vector<Point *> points;

		Octree * up_front_left;
		Octree * up_front_right;
		Octree * up_back_left;
		Octree * up_back_right;

		Octree * down_front_left;
		Octree * down_front_right;
		Octree * down_back_left;
		Octree * down_back_right;

		Octree(float x, float y, float z, float size,int max_points,int max_depth);
		~Octree();
		void Octree::destructRecursive();

		bool insertPoint(Point * p);
		bool isOutside(std::vector<float> pos,float size,Point * p);
		void subdivide();
		
		//o argumento count é usado para saber quantos pontos sao testados
		std::vector<Point*> queryOctree(std::vector<float> center, float size,int & count);
		bool intersects(std::vector<float> center, float size);

		void draw();

	private:
		float right, left, top, bot, front, back;
};