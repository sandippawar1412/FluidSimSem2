#include "GridStag.h"
//#include "commonData.h"
extern double zoomFactor;
GridStag :: GridStag()
{
	nX = 30;
	nY = 30;
	
	dx = (double)(zoomFactor-0)/nX;
	dy = (double)(zoomFactor-0)/nY;

		
	u.resize(nY,nX+1);
	v.resize(10,10);
	p.resize(nX,nY);
	d.resize(nX,nY);
	
	u0.resize(nY,nX+1);
	v0.resize(nY+1,nX);
	p0.resize(nX,nY);
	distanceLevelSet.resize(nX, nY);
	d0.resize(nX,nY);
	u.clear();
	v.clear();
	p.clear();
	d.clear();
	
	u0.clear();
	v0.clear();
	p0.clear();
	d0.clear();
	distanceLevelSet.clear();
	cellType.resize(nX,nY);
	cellType.clear();
	isFluidBoundary.resize(nY, nX);
	isFluidBoundary.clear();
}

void GridStag :: initGridStag()
{
	nX = nY = 64;
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
	u.clear();
	v.clear();
	p.clear();
	d.clear();
	
	u0.clear();
	v0.clear();
	p0.clear();
	d0.clear();
	
	cellType.resize(nY,nX);
	cellType.clear();
	isFluidBoundary.resize(nY, nX);
	isFluidBoundary.clear();
	distanceLevelSet.resize(nY, nX);
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
