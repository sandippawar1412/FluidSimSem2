//GridStag.h

#ifndef _GRIDSTAG_H_
#define _GRIDSTAG_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <vector>
#include "Particles.h"
using namespace boost::numeric::ublas;

#include "ParameterFLAGS.hpp"


class GridStag
{
public:
	int gridSizeVL;//grid Size Vertical line - physical
	int gridSizeHL; //grid Size Horizontal Lines - Physical
	double stepX ; //associated with gridSizeVL
	double stepY; //associated with gridSizeHL
	
	int nX; //Number of vertical lines
	int nY; //Number if Horizontal Lines..
	
	double dx;
	double dy;
	
	matrix<double> u; //VL - 
	matrix<double> v; //HL
	
	matrix<double> p; //pressure
	matrix<double> d; //density
	matrix<double> cellType; //type of cell..liquid/air/solid
	matrix<double> isFluidBoundary; //type of cell..liquid/air/solid

	matrix<double> u0; //VL - 
	matrix<double> v0; //HL
	
	matrix<double> p0; //pressure
	matrix<double> d0; //density
	matrix<double> distanceLevelSet;
	std::vector< Particles* > fluidParticles ;
	void initGridStag();
	bool isParticlePresent(double mc_x, double mc_y) ;
};



#endif
