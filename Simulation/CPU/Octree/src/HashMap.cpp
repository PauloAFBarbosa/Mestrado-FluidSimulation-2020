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
		this->offsets = new int[size];
		for (int i = 0; i < size; i++)
		{
			this->offsets[i] = 0;
		}
		this->particles = (Point*)malloc(sizeof(Point) * particleNumber);
		this->numbParticles = particleNumber;
	}
	
}

void HashMap::computeOffsets(){
	this->offsets[0] = 0;

	for (int i = 1; i < this->size; i++)
	{
		this->offsets[i] = this->offsets[i - 1] + this->bucketSizes[i - 1];
	}
}

unsigned int HashMap::hashFunction(float pos[3],float H,int size) {
	unsigned int m = 100;
	unsigned int p1 = 2693;
	unsigned int p2 = 3163;
	unsigned int p3 = 4091;

	int part1 = ((int)(pos[0] / H) * p1);
	int part2 = ((int)(pos[1] / H) * p2);
	int part3 = ((int)(pos[2] / H) * p3);
	
	unsigned int ret = (unsigned int)(part1 ^ part2 ^ part3) % size;
	
	return ret;
}

unsigned int HashMap::hashFunctionV2(float pos[3], float H, int size) {

	int newx = ((int)(pos[0] / H));
	int newy = ((int)(pos[1] / H));
	int newz = ((int)(pos[2] / H));

	unsigned int ret = (unsigned int)(newx+newy*size+newz*size)%size;

	return ret;
}

bool contains(unsigned int arr[27], int size, unsigned int member) {
	bool ret = false;

	for (int i = 0; i < size; i++)
	{
		if (arr[i] == member)
			ret = true;
	}
	return ret;
}

//Vai devolver um array com 27 indices, dos buckets adjacentes ao bucket que é mandado, inclusive manda o bucket que recebeu.
//Com as colisões varias posições podem ir para o mesmo bucket, quando isso acontece eu estou a adicionar varias vezes o mesmo bucket.
//Assim algumas particulas vao ter um peso dobradou (ou mais)
//Tenho de remover buckets repetidos
//Vai devolver o tamanho to array que devolve---> devido as colisões as vezes o array vai ter menos que 27, para nao adicionar buckets repetidos
int HashMap::getAdj(float pos[3], float H, unsigned int ret[27]) {
	
	unsigned int p1 = 2693;
	unsigned int p2 = 3163;
	unsigned int p3 = 4091;
	int part1, part2, part3;
	int cell[3];
	cell[0] = (int)(pos[0] / H);
	cell[1] = (int)(pos[1] / H);
	cell[2] = (int)(pos[2] / H);

	unsigned int tempBucket;
	int retSize = 0;

	//printf("Recebi a posicao %f %f %f\n", pos[0], pos[1], pos[2]);
	//printf("calculei a cell %d %d %d\n", cell[0], cell[1], cell[2]);
	
	//bucket recebido
	part1 = ((int)cell[0] * p1);
	part2 = ((int)cell[1] * p2);
	part3 = ((int)cell[2] * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	ret[0] = tempBucket;
	
	retSize++;
	//bucket right
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1]) * p2);
	part3 = ((int)(cell[2]) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret,retSize,tempBucket)==false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1]) * p2);
	part3 = ((int)(cell[2]) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket right front
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1]) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket right back
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1]) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket front
	part1 = ((int)(cell[0]) * p1);
	part2 = ((int)(cell[1]) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket back
	part1 = ((int)(cell[0]) * p1);
	part2 = ((int)(cell[1]) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left front
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1]) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left back 
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1]) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//-----O mesmo mas para os buckets de cima
	//bucket recebido cima
	part1 = ((int)(cell[0]) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2]) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket right cima
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2]) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left cima
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2]) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket right front cima
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket right back cima
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket front cima
	part1 = ((int)(cell[0]) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket back cima
	part1 = ((int)(cell[0]) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left front cima
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left back cima 
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1] + 1) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}


	//-----O mesmo mas para os buckets de baixo
	//bucket recebido baixo
	part1 = ((int)(cell[0]) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2]) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket right baixo
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2]) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left baixo
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2]) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket right front baixo
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket right back baixo
	part1 = ((int)(cell[0] + 1) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket front baixo
	part1 = ((int)(cell[0]) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket back baixo
	part1 = ((int)(cell[0]) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left front baixo
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2] - 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	//bucket left back baixo 
	part1 = ((int)(cell[0] - 1) * p1);
	part2 = ((int)(cell[1] - 1) * p2);
	part3 = ((int)(cell[2] + 1) * p3);
	tempBucket = (unsigned int)(part1 ^ part2 ^ part3) % this->size;
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = tempBucket;
		retSize++;
	}

	

	return retSize;
	//for (int i = 0; i < 27; i++)
	//{
	//	printf("Os buckets que vai devolver sao %d\n", ret[i]);
	//}
}


int HashMap::getAdjV2(float pos[3], float H, unsigned int ret[27]) {
	
	float tempPos[3];
	unsigned int tempBucket;
	int retSize = 0;
	//bucket recebido
	ret[0] = this->hashFunction(pos, H, this->size);
	retSize++;
	//bucket right
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1];
	tempPos[2] = pos[2];
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1];
	tempPos[2] = pos[2];
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket right front
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1];
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket right back
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1];
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket front
	tempPos[0] = pos[0];
	tempPos[1] = pos[1];
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket back
	tempPos[0] = pos[0];
	tempPos[1] = pos[1];
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left front
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1];
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left back 
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1];
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//-----O mesmo mas para os buckets de cima
	//bucket recebido cima
	tempPos[0] = pos[0];
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2];
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket right cima
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2];
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left cima
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2];
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket right front cima
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket right back cima
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket front cima
	tempPos[0] = pos[0];
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket back cima
	tempPos[0] = pos[0];
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left front cima
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left back cima 
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1] + H;
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}


	//-----O mesmo mas para os buckets de baixo
	//bucket recebido baixo
	tempPos[0] = pos[0];
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2];
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket right baixo
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2];
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left baixo
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2];
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket right front baixo
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket right back baixo
	tempPos[0] = pos[0] + H;
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket front baixo
	tempPos[0] = pos[0];
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket back baixo
	tempPos[0] = pos[0];
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left front baixo
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2] - H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	//bucket left back baixo 
	tempPos[0] = pos[0] - H;
	tempPos[1] = pos[1] - H;
	tempPos[2] = pos[2] + H;
	tempBucket = this->hashFunction(tempPos, H, this->size);
	if (contains(ret, retSize, tempBucket) == false) {
		ret[retSize] = this->hashFunction(tempPos, H, this->size);
		retSize++;
	}

	return retSize;
	//for (int i = 0; i < 27; i++)
	//{
	//	printf("Os buckets que vai devolver sao %d\n", ret[i]);
	//}
}
void HashMap::addParticle(Point* p,float H,int offset) {
	float density;
	float pressure;

	glm::vec3 viscosity;

	unsigned int bucket = this->hashFunction(p->pos, H,this->size);

	//[offset + bucket_size]
	//bucket_size dentro do h conta quantos elementos tem atualmente dentro do bucket
	//printf("offset %d bucketSizes %d bucket %d\n", offset, h->bucketSizes[bucket],bucket);
	this->particles[offset + this->bucketSizes[bucket]].id = p->id;
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

void HashMap::updateHashMap(float H) {
	
	int * bucketSizes  = new int[this->size];
	Point* particlesTemp = (Point*)malloc(sizeof(Point) * this->numbParticles);
	//vai inicializar o bucketsize temporario a zero (o que vai contar qual o tamanho de cada bucket, para saber o offset)
	//Da reset aos bucketsizes do hashmap
	for (size_t i = 0; i < this->size; i++)
	{
		bucketSizes[i] = 0;
		this->bucketSizes[i] = 0;
	}
	//copia as particulas do hashmap para um array temporario, e vai contando os bucketsizes
	//tenho de mudar a particula para um array temporario porque ao mover a particula A para o lugar da particula B, a particula B vai se perder.
	for (int i = 0; i < this->numbParticles; i++)
	{
		particlesTemp[i] = this->particles[i];
		bucketSizes[this->hashFunction(particlesTemp[i].pos, H,this->size)]++;
	}

	//Agora vai correr todas as particulas de novo e vai adicionar ao hashmap
	// e no final vai hazer o compute dos offsets
	for (int i =0 ; i <this->numbParticles ; i ++)
	{
		float pos[3];
		pos[0] = particlesTemp[i].pos[0];
		pos[1] = particlesTemp[i].pos[1];
		pos[2] = particlesTemp[i].pos[2];
		unsigned int bucket = this->hashFunction(pos, H, this->size);
		
		int offset = 0;
		//vai correr o array com os tamanhos dos buckets até chegar ao bucket atual
		for (int j = 0; j < bucket; j++)
		{
			//exemplo, se tivermos um array com [p1,p2,p3,p4]
			//p1 e p2 sao do bucket 0 e p3 e p4 do bucket 1
			//o bucket 1 começa no indice 2 que é o bucketsize do bucket 0
			offset += bucketSizes[j];
		}

		this->addParticle(&particlesTemp[i], H, offset);

	}
	this->computeOffsets();

	delete bucketSizes;
	delete particlesTemp;
}

int HashMap::getAdjBruteForce(float pos[3], float H, Point * ret) {
	int pointsAdded = 0;
	for (int i = 0; i < this->numbParticles; i++)
	{
		Point p1 = this->particles[i];
		float dist[3];
		dist[0] = pos[0] - p1.pos[0];
		dist[1] = pos[1] - p1.pos[1];
		dist[2] = pos[2] - p1.pos[2];
		float length = (sqrt(pow(dist[0], 2) + pow(dist[1], 2) + pow(dist[2], 2)));
		if (length < H) {
			ret[pointsAdded] = p1;
			pointsAdded++;
		}
	}
	return pointsAdded;
}