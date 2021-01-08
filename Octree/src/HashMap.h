#pragma once

#include "../src/Point.h"
#include <vector>

class HashMap {


	public:
		int numbParticles;
		int size;
		int * bucketSizes;	// vay ser um array com tamanho size que vai dizer quantas particulas tem em cada bucket
		Point * particles;	//array que vai guardar todas as particulas


		HashMap::HashMap(int version,int size, int particleNumber);
		unsigned int HashMap::hashFunction(double pos[3],double H);
		void HashMap::addParticle(Point* p,double H,int offset);

	private:
		float right, left, top, bot, front, back;
};