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
		unsigned int HashMap::hashFunction(double pos[3],double H,int size);
		unsigned int HashMap::hashFunctionV2(double pos[3], double H, int size);
		void HashMap::addParticle(Point* p,double H,int offset);
		int HashMap::getAdj(double pos[3], double H, unsigned int ret[27]);
		int HashMap::getAdjV2(double pos[3], double H, unsigned int ret[27]);
		int HashMap::getAdjBruteForce(double pos[3], double H, Point* ret);
		void HashMap::computeOffsets();
		void HashMap::updateHashMap(double H);

	private:
		
};