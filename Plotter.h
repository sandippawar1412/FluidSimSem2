/*
 * Plotter.h
 *
 *  Created on: 30-Sep-2013
 *      Author: srp
 */

#ifndef PLOTTER_H_
#define PLOTTER_H_

#include "GridStag.h"
#include "Printer.h"
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

#define forEach(i,start,end) for(int i=start;i<=end;i++)
#define SET2DINT(arr, val , noofFiles, ni) forEach(i,0,noofFiles-1) forEach(j,0,ni-1) arr[i][j]=val ;

#define MAX_VEL 10000
#define MIN_VEL -10000
class Plotter{

public:
	double maxU,maxV,minU,minV;
	int ni ;
	int noofFiles ;
	int maxIter ;
	GridStag* sGrid ;
	Printer* prnt;
	void init(GridStag*sGrid,int ni,int noofFiles,int maxIterations) ;
	void prepareData(matrix<double>u, matrix<double>v, int it) ;
	void printMaxMin();
	void createDatFile(); //number if iterations considered
	void createDatFile_Helper(char fname[], double min, double max, int *arr);
	int getIndex(double key,double min,double max, int ni);
};



#endif /* PLOTTER_H_ */
