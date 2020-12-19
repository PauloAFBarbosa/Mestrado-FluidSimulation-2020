#pragma once
#include "../glm/vec3.hpp"
#include "../src/Point.h"
#include <vector>

class Octree {


	public:
		glm::vec3 pos;
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
		bool isOutside(glm::vec3 pos,float size,Point * p);
		void subdivide();
		
		//o argumento count é usado para saber quantos pontos sao testados
		void Octree::queryOctree(glm::vec3 center, float size, int& count, std::vector<Point*> * ret);
		bool intersects(glm::vec3 center, float size);

		void draw();
		void drawOut();

	private:
		float right, left, top, bot, front, back;
};