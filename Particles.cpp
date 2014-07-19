#include "Particles.h"

Particles::Particles(double x, double y) {
   
    this->x = x;
    this->y = y;
}

Particles::~Particles() {
}

bool Particles:: isEqual(double x,double y)
{
	if(this->x == x && this->y == y ){
		return true;
	}
	return false;	
}
