//Printer.cpp
#include "Printer.h"

void Printer :: init(GridStag* sGrid)
{
	this->sGrid = sGrid;
}

void Printer :: matrices(matrix<double> mat)
{
	cout<<endl;	
	for( int i = (int)(mat.size1()-1);i >= 0; i--,cout<<endl) //value returned by size1 & size2 is of type unsigned int..
		for( int j=0;j < (int)mat.size2();j++,cout<<" ")
			/*if(mat(i,j)>0 && mat(i,j)<1)
				cout<<"9";
			*/
			cout<<mat(i,j);
	cout<<endl;
}

void Printer :: msg(char* str)
{
	cout<<endl<<str<<endl;
}
void Printer :: toFile( char *msg=(char*)"Dummy1",  char fname[]="defaultPrint.log",char ch='a')
{
	ofstream out;
	if(ch=='a')
		out.open(fname,ios::app) ;
	else
		out.open(fname) ;
	out<<msg<<endl ;
	out.close();
}

void Printer :: toFile( matrix<double> mat,  char fname[]="defaultPrint.log",char ch='a',char ch1='z')//zero allowed
{
	ofstream out;
	if(ch=='a')
		out.open(fname,ios::app) ;
	else
		out.open(fname) ;
	out<<endl;
	for( int i = (int)(mat.size1()-1);i >= 0; i--) //value returned by size1 & size2 is of type unsigned int..
	{
		bool flagEndl = false;
		for( int j=0;j < (int)mat.size2();j++)

			if((ch1=='n' && mat(i,j)!=0) || ch1 !='n'){
				out<<mat(i,j)<<" ";
				flagEndl = true;
			}
		if(flagEndl)
			out<<endl;
	}

	out<<endl;

	out.close();
}
