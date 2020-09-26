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

		bool insertPoint(Point * p);
		bool isOutside(std::vector<float> pos,float size,Point * p);
		void subdivide();

		void draw();

	private:
};