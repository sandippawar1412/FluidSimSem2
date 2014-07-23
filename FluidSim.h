//FluidSim.h

#ifndef _FLUIDSIM_H
#define _FLUIDSIM_H

#include "GridStag.h"
#include <vector>
#include "Particles.h"
#include "Vec.h"

#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"

using namespace std;

//double force = 5.0f ;
//double density = 100.0f ;
#define PARTICLE_PER_CELL 4
#define FLUID 1
#define SOLID 2
#define AIR 0

//#define
#define DENSITY 100
#define FORCE -9.8

#define gravity -9.8

#define VELOCITY_BC 1
#define VELOCITY_BC2 4
#define DENSITY_BC 2
#define PRESSURE_BC 3

//enum ProjectFlag{GSRflag = 0, PCGflag, PCGBridsonflag} ;

//ProjectFlag ProjectRoutine1 = PCGBridsonflag;

#define GSRflag 0
#define PCGflag 1
#define PCGBridsonflag 2
#define NoPressure 3

#define U(i,j) ((u(i,j)+u(i,j+1))/2)
#define V(i,j) ((v(i,j)+v(i+1,j))/2)
#define ON 1
#define OFF 0
#define forEach(i,start,end) for(int i=start;i<=end;i++)
class FluidSim{

public:
	GridStag *sGrid; //can be used as call by name...declare GridStag& sGrid
	matrix<int> uValid , vValid ; //stores the flag for valid velocities bordering the liquid cells
	double dt;
	double getMinDistance(double,double,double);
	void calculateLevelSetDistance();

	void init(GridStag* );
	void simulate(double);
	void applyBoundaryConditions(int );

	double cfl();
	matrix<double> advect2DSelf(matrix<double> q, double,matrix<double> u,matrix<double> v,int ); //passing u and v compenents...
	matrix<double> addForce(matrix<double> dest, double dt, matrix<double> src) ;
	//void projectPCG();	//Preconditioned Conjugate Gradient..
	//matrix<double> projectGSr(matrix<double> u,matrix<double> v, double dt);	//Gauss-Seidel relaxation//applies pressure and returns as well

	void initFluidBody_Helper(int bx,int tx, int by,int ty,matrix<double>& mat,double val);
	void initFluidBody(int ); //set markers...add density/particles in region

	void advectParticles(std::vector< Particles* > & plst, matrix<double> ,matrix<double>, double dt );
	
	matrix<double> addGravity(matrix<double>, double ); //ua=un+dt*g	
	
	void initSolidBoundary(int ); //set markers..applyboundarycondition..can change the boundary...by giving choice
	
	void markFluidCells();	
	//void extrapolate();
	//void extrapolate4();
	void extrapolate2D(matrix<double> &grid, matrix<int> &valid);
	//helper to extrapolate
	int isInBound(int i, int j, int l, int nY, int nX, int val);
	
	Vec2d trace_rk2(const Vec2d& position, double dt) ;
	Vec2d get_Velocity(const Vec2d& position) ;
	void advect(double);

	double getVelInterpolated(double x,double y, matrix<double> mat);
	
	void solvePressureBridson(float dt);
	//Solver data
   PCGSolver<double> solver;
   SparseMatrixd matrix1;
   std::vector<double> rhs;
   std::vector<double> pressure;
   
   //for tracking advection
   void advect2DCell(int i,int j, double dt, matrix<double> u, matrix<double> v, int component);
   void advect2DParticle(double i,double j, double dt, matrix<double> u, matrix<double> v);
   void RK2(double &posx, double &posy,matrix<double> u, matrix<double> v, double dt);
   
   
   void setValidVelocity(int val);
   void addViscosity(double,double) ;
   matrix<double> addVisc_Helper(matrix<double>, double, double , int );
   int getNeighbours(matrix<double>mat, int comp, int i, int j, double &nval);

   void setVelocityError(matrix<double> &mat);
};

#endif
