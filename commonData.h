
//#include <vector>
#include <iostream>
using namespace std;
double zoomFactor = 1.0;

int flagC = 0;
int sb=0			; //state Boundary layer..
int flagApproximateFlow=0;
int componentFlag = 1;
int traceCellFlag = 0;
#include <vector>
std::vector<int> traceAdvectCells(4) ;
std::vector<double> traceAdvectCellsd(4) ;

//vector<int> traceAdvectCells(4) ;

#define FRAME_INTERVAL_ms 40  //increase to speed down,decrease to speed up
//time interval between 2 frames = 40 ms
//total frames covered in 1s = 1000ms/40ms = 25
