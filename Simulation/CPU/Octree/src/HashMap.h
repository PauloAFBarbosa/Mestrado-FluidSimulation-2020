#pragma once

#include "../src/Point.h"
#include <vector>

class HashMap {


	public:
		int numbParticles;
		int size;
		int * bucketSizes;	// vay ser um array com tamanho size que vai dizer quantas particulas tem em cada bucket
		int * offsets; // vai guardar o offset para aceder a cada bucket. Basicamente é a soma do bucketsizes, ou seja, offsets[0]=0 offsets[1]=bucketsizes[0] offsets[2]=offsets[1]+bucketsizes[1] ... 
		Point * particles;	//array que vai guardar todas as particulas


		HashMap::HashMap(int version,int size, int particleNumber);
		unsigned int HashMap::hashFunction(float pos[3],float H,int size);
		unsigned int HashMap::hashFunctionV2(float pos[3], float H, int size);
		void HashMap::addParticle(Point* p,float H,int offset);
		int HashMap::getAdj(float pos[3], float H, unsigned int ret[27]);
		int HashMap::getAdjV2(float pos[3], float H, unsigned int ret[27]);
		int HashMap::getAdjBruteForce(float pos[3], float H, Point* ret);
		void HashMap::computeOffsets();
		void HashMap::updateHashMap(float H);

	private:
		
};