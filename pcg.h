//pcg.h
#ifndef _PCG_H
#define _PCG_H
#include "opengl.h"
#include "GridStag.h"
#define LIQUID 1
#define SOLID 2
#define AIR 0

//#define sqr(a) (a*a)

class PCG{
public:
	GridStag* sGrid;
	double dt;
	
	matrix<double> rhs ;
	matrix<double> Adiag ;
	matrix<double> Aplusi ;
	matrix<double> Aplusj ;
	matrix<double> preCon ;
	matrix<double> s,z,q ;
	matrix<double> b;
	
	//void init(GridStag*);// sets sGrid.. allocates memory to all matrices
	
	void reserveMem();
	
	void solvePressure2D(GridStag*, double);
	
	//PCG algo
 
    void calculateDivergence();  //fig 4.2..result in rhs 
    void formA();  //fig 4.4 
    bool solvePressure(); //fig 4.5 //main PCG algo
	void formPreconditioner(); //fig 4.6
    void applyPreconditioner();//fig 4.7 //in 4.5
    void applyPressure();//fig 4.1 //update the u's and v's 
    void applyA();    
    //helper functions
    double dotProduct(matrix<double> ,matrix<double> );
    
    double getMaxR();
  
	//helper
	double sqr(double )	;
	
	void computeResidual();
    
};

#endif
