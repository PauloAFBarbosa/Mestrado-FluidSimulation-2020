#include "HashMap.h"
#include "../toolkits/glut/GL/glut.h"

HashMap::HashMap(int version,int size, int particleNumber)
{
	//versão 1 -> cada bucket pode guardar todas as particulas
	if (version == 1) {

	}
	//versão 2 -> Cada bucket guarda apenas as suas particulas
	//Ao criar as particulas, vai ter de correr 2 vezes. Uma para contar quantas particulas vai ter cada bucker, outra para as meter dentro do hashmap
	else if (version == 2) {
		this->size = size;
		this->bucketSizes = new int[size];
		for (int i = 0; i < size; i++)
		{
			this->bucketSizes[i] = 0;
		}
		printf("Size de 1 point e %d\n",sizeof(Point));
		printf("Size double %d\n", sizeof(double));
		printf("Size vec3 %d\n", sizeof(glm::vec3));
		this->particles = (Point*)malloc(sizeof(Point) * particleNumber);
		this->numbParticles = particleNumber;
	}
	
}

unsigned int HashMap::hashFunction(double pos[3],double H) {
	unsigned int m = 100;
	unsigned int p1 = 2693;
	unsigned int p2 = 3163;
	unsigned int p3 = 4091;

	int part1 = ((pos[0] / H) * p1);
	int part2 = ((pos[1] / H) * p2);
	int part3 = ((pos[2] / H) * p3);

	unsigned int ret = (part1 ^ part2 ^ part3) % m;
	
	return ret;
}

void HashMap::addParticle(Point* p,double H,int offset) {
	double density;
	double pressure;

	glm::vec3 viscosity;

	unsigned int bucket = this->hashFunction(p->pos, H);

	//[offset + bucket_size]
	//bucket_size dentro do h conta quantos elementos tem atualmente dentro do bucket
	//printf("offset %d bucketSizes %d bucket %d\n", offset, h->bucketSizes[bucket],bucket);
	//pos
	this->particles[offset + this->bucketSizes[bucket]].pos[0] = p->pos[0];
	this->particles[offset + this->bucketSizes[bucket]].pos[1] = p->pos[1];
	this->particles[offset + this->bucketSizes[bucket]].pos[2] = p->pos[2];

	//originalPos
	this->particles[offset + this->bucketSizes[bucket]].originalPos[0] = p->originalPos[0];
	this->particles[offset + this->bucketSizes[bucket]].originalPos[1] = p->originalPos[1];
	this->particles[offset + this->bucketSizes[bucket]].originalPos[2] = p->originalPos[2];

	//color
	this->particles[offset + this->bucketSizes[bucket]].color[0] = p->color[0];
	this->particles[offset + this->bucketSizes[bucket]].color[1] = p->color[1];
	this->particles[offset + this->bucketSizes[bucket]].color[2] = p->color[2];

	this->particles[offset + this->bucketSizes[bucket]].force[0] = p->force[0];
	this->particles[offset + this->bucketSizes[bucket]].force[1] = p->force[1];
	this->particles[offset + this->bucketSizes[bucket]].force[2] = p->force[2];

	this->particles[offset + this->bucketSizes[bucket]].velocity[0] = p->velocity[0];
	this->particles[offset + this->bucketSizes[bucket]].velocity[1] = p->velocity[1];
	this->particles[offset + this->bucketSizes[bucket]].velocity[2] = p->velocity[2];

	this->particles[offset + this->bucketSizes[bucket]].velocityEval[0] = p->velocityEval[0];
	this->particles[offset + this->bucketSizes[bucket]].velocityEval[1] = p->velocityEval[1];
	this->particles[offset + this->bucketSizes[bucket]].velocityEval[2] = p->velocityEval[2];

	this->particles[offset + this->bucketSizes[bucket]].acceleration[0] = p->acceleration[0];
	this->particles[offset + this->bucketSizes[bucket]].acceleration[1] = p->acceleration[1];
	this->particles[offset + this->bucketSizes[bucket]].acceleration[2] = p->acceleration[2];

	this->particles[offset + this->bucketSizes[bucket]].gravity[0] = p->gravity[0];
	this->particles[offset + this->bucketSizes[bucket]].gravity[1] = p->gravity[1];
	this->particles[offset + this->bucketSizes[bucket]].gravity[2] = p->gravity[2];

	this->particles[offset + this->bucketSizes[bucket]].surfaceNormal[0] = p->surfaceNormal[0];
	this->particles[offset + this->bucketSizes[bucket]].surfaceNormal[1] = p->surfaceNormal[1];
	this->particles[offset + this->bucketSizes[bucket]].surfaceNormal[2] = p->surfaceNormal[2];

	this->particles[offset + this->bucketSizes[bucket]].surfaceTension[0] = p->surfaceTension[0];
	this->particles[offset + this->bucketSizes[bucket]].surfaceTension[1] = p->surfaceTension[1];
	this->particles[offset + this->bucketSizes[bucket]].surfaceTension[2] = p->surfaceTension[2];

	this->particles[offset + this->bucketSizes[bucket]].viscosity[0] = p->viscosity[0];
	this->particles[offset + this->bucketSizes[bucket]].viscosity[1] = p->viscosity[1];
	this->particles[offset + this->bucketSizes[bucket]].viscosity[2] = p->viscosity[2];

	this->particles[offset + this->bucketSizes[bucket]].density = p->density;
	this->particles[offset + this->bucketSizes[bucket]].pressure = p->pressure;

	this->bucketSizes[bucket]++;

}