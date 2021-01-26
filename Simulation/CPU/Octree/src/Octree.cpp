#include "Octree.h"
#include "../toolkits/glut/GL/glut.h"

Octree::Octree(double x, double y, double z, double size, int max_points,int max_depth)
{

	
	this->pos[0] = x;
	this->pos[0] = y;
	this->pos[0] = z;
	this->size = size;
	this->devided = false;
	this->max_points=max_points;
	
	std::vector<Point *> temp2;
	this->points=temp2;

	this->up_front_left=NULL;
	this->up_front_right=NULL;
	this->up_back_left=NULL;
	this->up_back_right=NULL;

	this->down_front_left=NULL;
	this->down_front_right=NULL;
	this->down_back_left=NULL;
	this->down_back_right=NULL;

	this->max_depth = max_depth;
}

void Octree::destructRecursive() {
	if (this->devided) {
		this->up_front_left->destructRecursive();
		this->up_front_right->destructRecursive();
		this->up_back_left->destructRecursive();
		this->up_back_right->destructRecursive();

		this->down_front_left->destructRecursive();
		this->down_front_right->destructRecursive();
		this->down_back_left->destructRecursive();
		this->down_back_right->destructRecursive();
	}
	//delete this;
}

Octree::~Octree()
{
	if (this->devided) {
		delete this->up_front_left;
		delete this->up_front_right;
		delete this->up_back_left;
		delete this->up_back_right;

		delete this->down_front_left;
		delete this->down_front_right;
		delete this->down_back_left;
		delete this->down_back_right;
	}
	
	
}

bool Octree::insertPoint(Point * p){
	//Posso calcular estas boundaries logo no construtor 
	if (this->isOutside(this->pos,this->size,p))
	{
		return false;
	}
	//se tem menos pontos que o maximo, vai adicionar
	if (this->max_depth==0 || this->points.size()<this->max_points){
		
		this->points.push_back(p);
		return true;
	}
	else{
		
		if (this->devided==false){
			this->subdivide();
		}

		if (this->up_front_left->insertPoint(p)){
			return true;
		}
		else if(this->up_front_right->insertPoint(p)){
			return true;
		}
		else if(this->up_back_left->insertPoint(p)){
			return true;
		}
		else if(this->up_back_right->insertPoint(p)){
			return true;
		}
		else if(this->down_front_left->insertPoint(p)){
			return true;
		}
		else if(this->down_front_right->insertPoint(p)){
			return true;
		}
		else if(this->down_back_left->insertPoint(p)){
			return true;
		}
		else if(this->down_back_right->insertPoint(p)){
			return true;
		}	
	}
}
//return true if point is outside
bool Octree::isOutside(double pos[3],float size,Point * p){
	return (p->pos[0] > (pos[0] + size) || p->pos[0] < (pos[0] - size) ||
		p->pos[1] > (pos[1] + size) || p->pos[1] < (pos[1] - size) ||
		p->pos[2] > (pos[2] + size) || p->pos[2] < (pos[2] - size));
}

void Octree::subdivide() {

	this->up_front_left = new Octree(this->pos[0] - (this->size / 2), this->pos[1] + (this->size / 2), this->pos[2] - (this->size / 2), this->size / 2, this->max_points, this->max_depth - 1);
	this->up_front_right = new Octree(this->pos[0] + (this->size / 2), this->pos[1] + (this->size / 2), this->pos[2] - (this->size / 2), this->size / 2, this->max_points, this->max_depth - 1);
	this->up_back_left = new Octree(this->pos[0] - (this->size / 2), this->pos[1] + (this->size / 2), this->pos[2] + (this->size / 2), this->size / 2, this->max_points, this->max_depth - 1);
	this->up_back_right = new Octree(this->pos[0] + (this->size / 2), this->pos[1] + (this->size / 2), this->pos[2] + (this->size / 2), this->size / 2, this->max_points, this->max_depth - 1);

	this->down_front_left = new Octree(this->pos[0] - (this->size / 2), this->pos[1] - (this->size / 2), this->pos[2] - (this->size / 2), this->size / 2, this->max_points, this->max_depth - 1);
	this->down_front_right = new Octree(this->pos[0] + (this->size / 2), this->pos[1] - (this->size / 2), this->pos[2] - (this->size / 2), this->size / 2, this->max_points, this->max_depth - 1);
	this->down_back_left = new Octree(this->pos[0] - (this->size / 2), this->pos[1] - (this->size / 2), this->pos[2] + (this->size / 2), this->size / 2, this->max_points, this->max_depth - 1);
	this->down_back_right = new Octree(this->pos[0] + (this->size / 2), this->pos[1] - (this->size / 2), this->pos[2] + (this->size / 2), this->size / 2, this->max_points, this->max_depth - 1);

	this->devided = true;
}



void Octree::queryOctree(double center[3], float size,int & count, std::vector<Point*> * ret)
{

	if (!this->intersects(center, size)) {

		return;
	}
	else {
		for (size_t i = 0; i < this->points.size() ; i++)
		{
			
			
			
			if (!this->isOutside(center, size, this->points.at(i))) {
				
				//puts the point red
				//this->points.at(i)->color.x = 1;
				//this->points.at(i)->color.y = 0;
				//this->points.at(i)->color.z = 0;
				ret->push_back(this->points.at(i));
				count++;
			}
		}
		
		if (this->devided) {
			
			this->up_front_left->queryOctree(center, size,count,ret);
			//ret.insert(ret.end(),ret2.begin(),ret2.end());
			//ret2.clear();

			this->up_front_right->queryOctree(center, size, count,ret);
			//ret.insert(ret.end(), ret2.begin(), ret2.end());
			//ret2.clear();

			this->up_back_left->queryOctree(center, size, count, ret);
			//ret.insert(ret.end(), ret2.begin(), ret2.end());
			//ret2.clear();

			this->up_back_right->queryOctree(center, size, count, ret);
			//ret.insert(ret.end(), ret2.begin(), ret2.end());
			//ret2.clear();

			this->down_front_left->queryOctree(center, size, count, ret);
			//ret.insert(ret.end(), ret2.begin(), ret2.end());
			//ret2.clear();

			this->down_front_right->queryOctree(center, size, count, ret);
			//ret.insert(ret.end(), ret2.begin(), ret2.end());
			//ret2.clear();

			this->down_back_left->queryOctree(center, size, count, ret);
			//ret.insert(ret.end(), ret2.begin(), ret2.end());
			//ret2.clear();

			this->down_back_right->queryOctree(center, size, count, ret);
			//ret.insert(ret.end(), ret2.begin(), ret2.end());
			//ret2.clear();

			return ;
		}
		
		return ;
	}
}

//devolve true se houver interseçao
bool Octree::intersects(double center[3], float size)
{	//aplicar a lei de De Morgan faz isto ficar mais rapido
	//em vez de A ou b ou c fica
	// !A e !B e !C - assim mal um falhe ele nao precisa procurar mais nada
	return (!(center[0] - size > this->pos[0] + this->size) && //esq dir
		!(center[0] + size < this->pos[0] - this->size) && //dir esq
		!(center[1] - size > this->pos[1] + this->size) && //baixo cima
		!(center[1] + size < this->pos[1] - this->size) && //cima baixo
		!(center[2] - size > this->pos[2] + this->size) && //frente tras
		!(center[2] + size < this->pos[2] - this->size)  //tras frente
		);
}

void Octree::draw(){
	
	glPushMatrix();
	glTranslatef(this->pos[0],this->pos[1],this->pos[2]);
	glutWireCube(this->size*2);
	glPopMatrix();

	if (this->devided == true) {
		this->up_front_left->draw();
		this->up_front_right->draw();
		this->up_back_left->draw();
		this->up_back_right->draw();

		this->down_front_left->draw();
		this->down_front_right->draw();
		this->down_back_left->draw();
		this->down_back_right->draw();
	}
}


void Octree::drawOut() {

	glPushMatrix();
	glTranslatef(this->pos[0], this->pos[1], this->pos[2]);
	glutWireCube(this->size * 2);
	glPopMatrix();

}

