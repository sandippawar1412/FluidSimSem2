/*
 * Plotter.cpp
 *
 *  Created on: 30-Sep-2013
 *      Author: srp
 */

#include "Plotter.h"
#include <cstring>
void Plotter :: init(GridStag*sGrid,int ni, int noofFiles,int maxIter) {
	this->sGrid = sGrid ;
	this->ni = ni ;
	this->noofFiles = noofFiles ;
	this->maxIter = maxIter ;

	this->prnt = new Printer;
	prnt->init(sGrid);
}


void Plotter :: prepareData(matrix<double>u, matrix<double>v, int it) {
	 static double maxU=MIN_VEL,maxV=MIN_VEL ;
	 static double minU=MAX_VEL,minV=MAX_VEL ;
	 int modFac = maxIter/(noofFiles-1);
	//double maxUV=MAX_VEL,minUV=MIN_VEL ;
	char str[50]="\0";
	forEach(i, 0, (int)u.size1()-1 ){
		forEach(j, 0, (int)u.size2()-1 ){
			double temp = (u(i,j));//+u(i,j+1))/2 ;

			maxU = maxU<temp?temp:maxU ;
			minU = minU>temp?temp:minU ;
		}
	}
	if(!(it%modFac)){
		sprintf(str,"datFiles//.dataU%d",it/modFac) ;
			//prnt->toFile(u,str,'w');
		prnt->toFile(u,str,'w','n') ;//nonzero
	}

	forEach(i, 0, (int)v.size1()-1 ){
		forEach(j, 0, (int)v.size2()-1 ){
			double temp = (v(i,j));//+v(i+1,j))/2.0 ;
//			cout<<"V="<<v(i,j)<<" "<<v(i+1,j)<<" "<<temp<<endl;
			maxV = maxV < temp ? temp : maxV ;
			minV = minV > temp ? temp : minV ;
		}
	}
	if(!(it%modFac)){
		sprintf(str,"datFiles//.dataV%d",it/modFac) ;
		prnt->toFile(v,str,'w','n');
	}
	this->maxU = maxU ;
	this->maxV = maxV ;
	this->minU = minU ;
	this->minV = minV ;

}


void Plotter :: printMaxMin(){
	cout<<"U : "<<minU<<" "<<maxU<<endl;
	cout<<"V : "<<minV<<" "<<maxV<<endl;
}


void Plotter :: createDatFile(){
	char str[50]="\0";
	int **cntArr ;

	cntArr = new int*[noofFiles] ;
	forEach(i,0,noofFiles-1)
		cntArr[i] = new int[ni] ;
	//int cntArr[6][100] ;

	SET2DINT(cntArr, 0 , noofFiles, ni);

	double step = (maxU-minU)/ni ;
	forEach(i,0,noofFiles-1){ //timestep cnt
		sprintf(str,"datFiles//.dataU%d",i) ;
		createDatFile_Helper(str,minU,maxU,cntArr[i]);//process file using min and max and populate cntArr
	}
	strcpy(str,"datFiles//_dataU.dat") ;
	char msg[1000];
	forEach(i,0,ni-1){
		sprintf(msg,"%+lf - %+lf",minU+step*i,minU+step*(i+1)) ;
		forEach(j,0,noofFiles-1){
			sprintf(msg,"%s %3d",msg,cntArr[j][i]) ;
		}
		if(i==0)
			prnt->toFile(msg,str,'w');
		else
			prnt->toFile(msg,str,'a');
	}

	step = (maxV-minV)/ni ;
	//memset(cntArr, 0 , ni*noofFiles*sizeof(int) );
	SET2DINT(cntArr, 0 , noofFiles, ni);
	forEach(i,0,noofFiles-1){ //timestep cnt
		sprintf(str,"datFiles//.dataV%d",i) ;
		createDatFile_Helper(str,minV,maxV,cntArr[i]);//process file using min and max and populate cntArr
	}
	strcpy(str,"datFiles//_dataV.dat") ;
	forEach(i,0,ni-1){
		sprintf(msg,"%+lf - %+lf",minV+step*i,minV+step*(i+1)) ;
		forEach(j,0,noofFiles-1){
			sprintf(msg,"%s %3d",msg,cntArr[j][i]) ;
		}
		if(i==0)
			prnt->toFile(msg,str,'w');
		else
			prnt->toFile(msg,str,'a');
	}

}

void Plotter :: createDatFile_Helper(char str[],double min,double max, int *arr){

	//sprintf(str,"dataV%d.dat",i) ;

	FILE*fp = fopen(str,"r");
	double temp=0;
	while(!feof(fp)){
		fscanf(fp,"%lf",&temp);
		//cout<<temp<<" " ;
		int index = getIndex(temp,min,max,ni);
		arr[index]++;
	}
	fclose(fp);
}

int Plotter :: getIndex(double key, double min,double max,int ni)//number of intervals
{
	double step = (max-min)/ni ;

	forEach(i,1,ni){
		if( key>=(min + step*(i-1)) && key< (min+step*i) )
			return i-1;
		if(i==ni)
			return i-1 ;//i in 1 to 10->0 - 9
	}
	return -1;
}
