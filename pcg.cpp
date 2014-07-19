
#include "pcg.h"
#include "Printer.h"
/*
 * will have a one function that takes argument as grid...calculates pressure and 
 * stores in the data member of grid...
 * 
 * all the declaration required by pcg will be done in .h file....
 * 
 * grid needs to have following data member(as matrices)
 * 
 * pressure 
 * mark <define each cell...as liquid/solid/air>
 * 
 * 
 * CFG.h should declare following matrices..
 * 
 * rhs(2D)
 * Adiag(2D)
 * Aplusi(2D)
 * Aplusj(2D)
 * preConditioner(2D)
 * s,z,m(2D)
 * 
 * 2D---is of size (ni,nj)
 */ 

void PCG ::  computeResidual()
{
	matrix<double> p = sGrid->p;
	for (int y=1; y< sGrid->nY-1; y++)
		for (int x=1; x< sGrid->nX-1; x++){
			double temp  = Adiag(y,x)* p(y,x) + Aplusj(y,x)*p(y,x+1)
											  +	Aplusi(y,x)*p(y+1,x)
											  +	Aplusj(y,x-1)*p(y,x-1)
											  +	Aplusi(y-1,x)*p(y-1,x) ;
			b(y,x)-=temp;
	//		b(y,x) = b(y,x)<1e-3? 0: b(y,x);								  
		}
}
void PCG :: reserveMem()
{
	int nY = sGrid->nY;
	int nX = sGrid->nX;

	rhs.resize(nY,nX);
	Adiag.resize(nY,nX); 
	Aplusi.resize(nY,nX); 
	Aplusj.resize(nY,nX); 
	preCon.resize(nY,nX); 
	s.resize(nY,nX);
	z.resize(nY,nX);
	q.resize(nY,nX);
	
	rhs.clear();
	Adiag.clear();
	Aplusi.clear();
	Aplusj.clear();
	preCon.clear();
	s.clear();
	z.clear();
	q.clear();
}


void PCG :: solvePressure2D(GridStag* sGrid, double dt)
{
	this->sGrid = sGrid;	
	this->dt = dt;
	reserveMem();
	
    calculateDivergence();
    formA();
    formPreconditioner();
    solvePressure();
    applyPressure();
    computeResidual();
    
    //Printer * prnt = new Printer;
    //prnt->matrices(sGrid->cellType);
    //prnt->matrices(b);
    
    
}

//PCG algo
void PCG :: calculateDivergence()  //fig 4.2..result in rhs 
{
	rhs.clear();
    //double scale = 1/sGrid->dx;
    double scale = 1;
    for(int i = 0; i <= sGrid->nY-1; i++) 
		for(int j = 0; j <= sGrid->nX-1; j++) 
			if(sGrid->cellType(i,j)==LIQUID)
		{
			rhs(i,j) = -scale * ( (sGrid->u(i,j+1)- sGrid->u(i,j)) + (sGrid->v(i+1,j)- sGrid->v(i,j)) ) *(double)sGrid->dx/dt ;
		}		
	b = rhs;	
}	

void PCG :: formA()  //fig 4.4 
{
	
	double scale = 1;// 
	//double scale = dt/(1 * sGrid->dx * sGrid->dx) ;
	Adiag.clear();
    Aplusi.clear();
    Aplusj.clear();
	//cout<<scale<<endl;
    // setting up Adiag, Aplusi, Aplusj
    /*
     * books convention i-->u(i,..)-->Aplusi,,,,,my convention j-->u(...,j)-->Aplusj
     */ 
    for (int y=1; y< sGrid->nY-1; y++)
		for (int x=1; x< sGrid->nX-1; x++)
		{
			if (sGrid->cellType(y,x) == LIQUID && sGrid->cellType(y,x+1) == LIQUID)
			{
				Adiag(y,x) += scale;
				Adiag(y,x+1) += scale;
				Aplusj(y,x) = -scale;
			}	
			else if (sGrid->cellType(y,x) == LIQUID && sGrid->cellType(y,x+1) == AIR)	
				Adiag(y,x) += scale;
			
			if (sGrid->cellType(y,x) == LIQUID && sGrid->cellType(y+1,x) == LIQUID)
			{
				Adiag(y,x) += scale;
				Adiag(y+1,x) += scale;
				Aplusi(y,x) = -scale;
			}	
			else if (sGrid->cellType(y,x) == LIQUID && sGrid->cellType(y+1,x) == AIR)	
				Adiag(y,x) += scale;
			/*if(sGrid->cellType(y+1,x)!=SOLID)
					Adiag(y,x) +=scale;
			if(sGrid->cellType(y-1,x)!=SOLID)
					Adiag(y,x) +=scale;
			if(sGrid->cellType(y,x+1)!=SOLID)
					Adiag(y,x) +=scale;
			if(sGrid->cellType(y,x-1)!=SOLID)
					Adiag(y,x) +=scale;
			if(sGrid->cellType(y+1,x)==LIQUID)
					Aplusi(y,x) = -scale;
			if(sGrid->cellType(y,x+1)==LIQUID)
					Aplusj(y,x) =scale;*/
		}
}	

bool PCG :: solvePressure() //fig 4.5 //main PCG algo
{

	sGrid->p.clear();  //initial guess  p =0

	//following some lines are unclear
	double tol = 1e-5 * getMaxR();
	//cout<<endl<<getMaxR()<<" "<<tol;
	if (getMaxR() == 0.0 )
		return false;
	//follows book
	
	applyPreconditioner()	;
	s = z;
	double row = dotProduct(z,rhs);
	if (row == 0.0)
		return false;
		
	for ( unsigned int iter = 0; iter < 100 ; iter++ )
	{
		applyA();
		double alpha = row/dotProduct( z, s );		
			//cout<<endl<<getMaxR()<< " "<<iter<<" dot,rho = "<<dotProduct( z, s )<<" "<<row<<" "<<tol;
		for (int y=0; y< sGrid->nY; y++)
			for (int x=0; x< sGrid->nX; x++)
			{
				sGrid->p( y , x ) += alpha * s( y , x )	;
				rhs( y , x ) -= alpha * z ( y , x )	;
			}
		
		if (getMaxR() <= tol)	
			return true;
		applyPreconditioner(); //z=applyPreconditioner(r)
		double row_new = dotProduct( z , rhs );
		double beta = row_new / row ;
		
		for (int y=0; y<sGrid->nY; y++)
			for (int x=0; x<sGrid->nX; x++)
				s( y , x ) = z( y , x ) + beta * s( y , x ) ;  //s = z + beta*s
		row = row_new ;
    }
    return false;
}	

void PCG :: applyA()
{
	z.clear();
    for (int y=1; y< sGrid->nY-1; y++)
		for (int x=1; x< sGrid->nX-1; x++)
			if (sGrid->cellType(y,x) == LIQUID)
			{
				z( y , x ) = Adiag( y , x ) * s( y , x ) + 
							 Aplusj( y , x ) * s( y  , x + 1 ) + 
							 Aplusi( y , x ) * s( y + 1 , x  ) + 
							 Aplusi( y-1 , x  ) * s( y - 1 , x  ) +  /*check*/
							 Aplusj( y , x-1 ) * s( y  , x -1  ) ; 
			}
}

void PCG :: formPreconditioner() //fig 4.6
{
		 //double tow = 0.99;
	 double tow = 0.97;
     double row = 0.25;
     double e;
     preCon.clear();
     
    for (int y = 1; y < sGrid->nY - 1; y++)
		for (int x = 1; x < sGrid->nX - 1; x++)
		{
			if(sGrid->cellType( y , x ) == LIQUID )	
			{
				e = Adiag( y , x ) - sqr(Aplusj(y, x-1)*preCon( y, x-1 ) )
								   - sqr(Aplusi(y -1, x )*preCon( y - 1, x ) )
								   - tow * (Aplusj(y,x-1) * Aplusi(y, x-1 ) * sqr(preCon(y,x-1))
								         + Aplusi(y-1 ,x ) * Aplusj(y-1, x) * sqr(preCon(y-1,x)) );
				if (e < row*Adiag(y,x))
					e = Adiag(y,x);              //result changes somewhat..
				preCon(y,x) = 1 / sqrt(e); 
				//preCon(y,x) = 1 / sqrt(e + 1e-6); 
			}
		}
}	

double PCG :: sqr ( double a )
{
	return a*a;
}

void PCG :: applyPreconditioner()//fig 4.7 //in 4.5
{
	//first solve Lq=r (r is rhs)
	q.clear();
	z.clear();
	double t;
	for (int y=1; y< sGrid->nY-1; y++)
		for (int x=1; x< sGrid->nX-1; x++)
			if ( sGrid->cellType( y , x ) == LIQUID )
			{
				t = rhs(y,x) - Aplusj( y  , x-1 ) * preCon( y  , x-1 ) * q( y , x-1 )
							 - Aplusi( y-1 , x ) * preCon( y -1, x ) * q( y-1 , x );
				q( y , x ) = t * preCon( y , x );
			}
	
	//Now solve for L(trans)z=q			
	for ( int y = sGrid->nY - 2 ; y>0 ; y-- )
		for ( int x = sGrid->nX - 2 ; x>0 ; x-- )
				if ( sGrid->cellType( y , x ) == LIQUID )
				{
					t = q( y , x ) - Aplusj( y , x ) * preCon( y , x ) * z ( y  , x+1 )
								   - Aplusi( y , x ) * preCon( y , x ) * z ( y+1 , x  );
					z( y , x ) = t * preCon( y , x );
				}
}	

void PCG :: applyPressure()//fig 4.1 //update the u's and v's 
{
	//double scale = dt / (1 * sGrid->dx); // book define rho value before use					
	//double scale = 1;

	/*for (int y = 2; y <= sGrid->nX-2; y++)  //1 to nX-2 should be ok.. 
		for (int x = 1; x <= sGrid->nY-2; x++)
			if(sGrid->cellType( y-1 , x ) == LIQUID || sGrid->cellType( y , x ) == LIQUID )
			{
				sGrid->v ( y , x ) += scale *( sGrid->p( y , x ) - sGrid->p( y-1 , x ));
			}	
	for (int y = 1; y <= sGrid->nX-2; y++)  //1 to nX-2 should be ok.. 
		for (int x = 2; x <= sGrid->nY-2; x++)
			if(sGrid->cellType( y , x-1 ) == LIQUID || sGrid->cellType( y , x ) == LIQUID )
			{
				sGrid->u ( y , x ) += scale *( sGrid->p( y , x ) - sGrid->p( y, x-1 ));
			}	
			*/
	double scale = dt / (1 * sGrid->dx); // book define rho value before use					
	for(int y=1; y < sGrid->nY-1;y++)		
		for(int x=1; x < sGrid->nX-1;x++){
			if(sGrid->cellType(y,x) == LIQUID){
			//	sGrid->u(y,x) -= scale * sGrid->p(y,x);	//using these commented lines;;;result is not expected..				
			//	sGrid->u(y,x+1) += scale * sGrid->p(y,x);					
				sGrid->u(y,x+1) -= scale * (sGrid->p(y,x+1) - sGrid->p(y,x));					
			//	sGrid->v(y,x) -= scale * sGrid->p(y,x);					
			//	sGrid->v(y+1,x) += scale * sGrid->p(y,x);					
				sGrid->v(y+1,x) -= scale * (sGrid->p(y+1,x) - sGrid->p(y,x));			
			}
		}	
			 
	/*for(int y=0; y <sGrid->nY-1;y++)		
		for(int x=0;y<sGrid->nX-1;x++){
			if(sGrid->cellType(y,x) == SOLID){
				sGrid->u(y,x) = 0;
				sGrid->u(y,x+1) = 0;
				sGrid->v(y,x) = 0;
				sGrid->v(y+1,x) = 0;
			}
		}	 */
}	
        
    //helper functions
double PCG :: dotProduct(matrix<double> mat1 ,matrix<double> mat2)
{
	double product = 0.0;
	for( unsigned int i = 0; i < mat1.size1(); i++ ) 
		for( unsigned int j = 0; j < mat1.size2(); j++ ) 
			product += mat1( i , j ) * mat2( i , j ); //order of i,j not dependent
	return product;
}	
    
double PCG :: getMaxR()
{
	double max = 0.0;
	for( unsigned int i = 0; i < rhs.size1(); i++ ) 
		for( unsigned int j = 0; j < rhs.size2(); j++ ) 
			if(!( fabs(rhs(i, j) ) <= max))  //order of i,j not dependent
				max = fabs(rhs( i , j )) ;
				
	return max;			

}
