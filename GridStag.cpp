#include "GridStag.h"
extern double zoomFactor;

void GridStag :: initGridStag()
{
	nX = nY = GRID_SIZE;
	//nY = 8;
	dx = (double)zoomFactor/nX;
	dy = 1;//(double)gridSizeHL/nY;

	u.resize(nY,nX+1);
	v.resize(nY+1,nX);
	p.resize(nY,nX);
	d.resize(nY,nX);

	u0.resize(nY,nX+1);
	v0.resize(nY+1,nX);
	p0.resize(nY,nX);
	d0.resize(nY,nX);
	cellType.resize(nY,nX);
	isFluidBoundary.resize(nY, nX);
	distanceLevelSet.resize(nY, nX);

	u.clear();
	v.clear();
	p.clear();
	d.clear();
	u0.clear();
	v0.clear();
	p0.clear();
	d0.clear();
	cellType.clear();
	isFluidBoundary.clear();
	distanceLevelSet.clear();
	//fluidParticles.clear(); 
}

bool GridStag :: isParticlePresent(double mc_x, double mc_y) 
{
	for (unsigned int i=0;i< this->fluidParticles.size();i++)
		if (fluidParticles.at(i)->isEqual(mc_x,mc_y) == true )
			return true;
	return false;		
}
