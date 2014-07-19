
//Renderer.h
#ifndef _RENDERER_H_
#define _RENDERER_H_

#include "GridStag.h"
//#include "Particles.h"

#define MAX(a,b) a>b?a:b
#define SIGN(a) a<0?1:0

#define LIQUID 1
#define SOLID 2
#define AIR 0


#define U(i,j) ((u(i,j)+u(i,j+1))/2)
#define V(i,j) ((v(i,j)+v(i+1,j))/2)

#define ERRF 0.01
#define FLUID 1
#define SOLID 2
#define AIR 0

class Renderer{

public:
	GridStag* sGrid;
	
	double gridLBX;
	double gridLBY;
	
	double gridRTX;
	double gridRTY;
	
	double stepX;
	double stepY;
		
	
	void updateGrid(GridStag* );
	void initRenderer();
	void init(GridStag* );
	
	Renderer();

	void renderGrid();	//rener nX*nY grid wire frame...
	void renderBoundary();
	void renderSurfaceBoundary();
	void drawSquare(double ,double, double, double, int); //used to render a cell
	void drawSquareFilled(double ,double, double, double, double); //used to render a cell
	void drawSquareFilled(int, int, int, int, double); //used to render a cell
	
	void renderVector2D(matrix<double> u, matrix<double> v);// render vector with arrow..
	//void renderCellContent2D(); //can use renderMat(sGrid->cellType)
	
	
	double getMax(matrix<double> mat); //helper to renderVector2D
	void drawVectorAngnMag(double lbx,double lby,double rtx,double rty,bool usign,bool vsign,double mag,double max,double ang);//helper to renderVector2D

	void renderMat(matrix<double> , int ); //display matrix on the sGrid..mode = same color,intensity
	
	/*---for staggered with stam---*/
	
	void renderVel2D_Stam();
	void renderDen2D_Stam();
	void renderDensity(matrix<double>); 
	/*--particles---*/
	void renderParticles();
	
	/*debug the advect*/
	void renderParticle(double x, double y,double R,double G,double B,int psize) ;
		
};


#endif
