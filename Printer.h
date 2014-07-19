//Printer.h

#ifndef _PRINTER_H
#define _PRINTER_H

#include "opengl.h"
#include "GridStag.h"
#include <fstream>
//using namespace std ;
class Printer{
public:
	GridStag* sGrid;
	void init(GridStag* );
	void matrices(matrix<double> );
	void toFile( matrix<double>,  char* fname,char,char);
	void toFile( char* msg,  char* fname,char);
	void msg(char* );
};


#endif

