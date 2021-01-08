#pragma once

#include "../src/Point.h"
#include <vector>

class Octree {


	public:
		double pos[3];
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

		Octree(double x, double y, double z, double size,int max_points,int max_depth);
		~Octree();
		void Octree::destructRecursive();

		bool insertPoint(Point * p);
		bool isOutside(double pos[3],float size,Point * p);
		void subdivide();
		
		//o argumento count é usado para saber quantos pontos sao testados
		void Octree::queryOctree(double center[3], float size, int& count, std::vector<Point*> * ret);
		bool intersects(double center[3], float size);

		void draw();
		void drawOut();

	private:
		float right, left, top, bot, front, back;
};