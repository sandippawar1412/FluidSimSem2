
#include "Renderer.h"
#include "Printer.h"
#include "opengl.h"


using namespace std;
extern double zoomFactor;
extern int sb;
Renderer :: Renderer()
{
}

void Renderer :: init(GridStag* grid)
{
	this->sGrid = grid;
	int nX = (this->sGrid)->nX;
	int nY = (this->sGrid)->nY;
	
/*--checkIt--*/
	gridRTX = zoomFactor; //zoomFactorX- main.cpp - glOrtho
	gridRTY = zoomFactor; //zoomFactorY
	gridLBX = 0;
	gridLBY = 0;

	stepX = (gridRTX - gridLBX)/nX;
	stepY = (gridRTX - gridLBX)/nY;
}
void Renderer :: updateGrid(GridStag* grid)
{
	this->sGrid = grid;
}

void Renderer :: initRenderer()
{
	int nX = (this->sGrid)->nX;
	int nY = (this->sGrid)->nY;
	
	gridLBX = -9;
	gridLBY = -9;
	stepX = -1*(gridLBX)*2/nX;
	stepY = -1*(gridLBY)*2/nY;
}

void Renderer :: renderGrid()
{
	int nX = (this->sGrid)->nX;
	int nY = (this->sGrid)->nY;

	//cout<<"\nInRenderer:displayGrid"<<nX<< " " <<nY<<" "<<gridLBX<<" "<<gridLBY;
	double gridLBX = this->gridLBX+stepX;
	double gridLBY = this->gridLBY+stepY;
	int fillFlag = 0;
	
	for (int i=1;i< nX-1;i++){
		for (int j=1;j< nY-1;j++){
			glColor3f(1,0,0);//set fill color..used if fillFlag=true
			drawSquare(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,fillFlag);//1: Filled grid,0:wire Frame				
			gridLBX+=stepX;
		}
		gridLBY += stepY;
		gridLBX = this->gridLBX+stepX;
	}
	//renderBoundary();
}

void Renderer :: renderBoundary()
{
	int nX = (this->sGrid)->nX;
	int nY = (this->sGrid)->nY;
	
	//cout<<"\nInRenderer:displayGrid"<<nX<< " " <<nY<<" "<<gridLBX<<" "<<gridLBY;
	double gridLBX = this->gridLBX;
	double gridLBY = this->gridLBY;
	int fillFlag = 0;
	
	for (int i=0;i< nX;i++){
		for (int j=0;j< nY;j++){
			if(i==0+sb||j==0+sb||i==nX-1-sb||j==nY-1-sb)
			{
				glColor3f(0.4,0.4,0.4);
				drawSquare(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,fillFlag=3);//0,1,2has defined colors..for air liquid and solid cells
			}	
			gridLBX+=stepX;
		}
		gridLBY += stepY;
		gridLBX = this->gridLBX;
	}

}

void Renderer :: drawSquare(double lbx,double lby,double rtx,double rty,int fillFlag)
{
		//cout<<"\nRenderer: drawSquare";
		//cout<<"\n.."<<fillFlag;
		switch(fillFlag)
		{
			case AIR:
				glColor3f(0.0,0,0.0);
				break;
			case LIQUID:
				glColor3f(0.1,0.1,0.7); //liquid
				break;
			case SOLID:
				glColor3f(0.2,0.2,0);
				break;
				
		}
		/*following 2 conditions should work..but not working..*/
		//if(fillFlag==LIQUID || fillFlag==SOLID )//uncomment to not render air cells
		//if(fillFlag == LIQUID  )//uncomment to render only liquid..
		
		{
			glBegin(GL_QUADS);		
			glVertex2f(lbx, lby);
			glVertex2f(lbx, rty);
			glVertex2f(rtx, rty);
			glVertex2f(rtx, lby);
			glEnd();
		}
		if(!fillFlag)
		{
			glColor3f(.2,.2,.2); //unComment to show grid lines
			glBegin(GL_LINE_LOOP);		
			glVertex2f(lbx, lby);
			glVertex2f(lbx, rty);
			glVertex2f(rtx, rty);
			glVertex2f(rtx, lby);
			glEnd();
		}
}
void Renderer :: drawSquareFilled(int li, int lj, int ri, int rj, double fillVal)
{
	double lbx = gridLBX + stepX*li;
	double lby = gridLBY + stepY*lj;
	
	double rtx = gridLBX + stepX*ri;
	double rty = gridLBY + stepY*rj;

	glColor3f(fillVal,0,0); 
	glBegin(GL_QUADS);		
	glVertex2f(lbx, lby);
	glVertex2f(lbx, rty);
	glVertex2f(rtx, rty);
	glVertex2f(rtx, lby);
	glEnd();
}
void Renderer :: drawSquareFilled(double lbx,double lby,double rtx,double rty,double fillVal)
{
		
		{
		//	cout<<fillVal<<endl;
			if(fillVal!=-2.0)
				glColor3f(0,0,0-fillVal);
			if(fillVal>=0)
				glColor3f(fillVal/6,fillVal/6,fillVal/6);


			glBegin(GL_QUADS);
			glVertex2f(lbx, lby);
			glVertex2f(lbx, rty);
			glVertex2f(rtx, rty);
			glVertex2f(rtx, lby);
			glEnd();
		}
		if(0)
		{
			glColor3f(.2,.2,.2); //unComment to show grid lines
			glBegin(GL_LINE_LOOP);		
			glVertex2f(lbx, lby);
			glVertex2f(lbx, rty);
			glVertex2f(rtx, rty);
			glVertex2f(rtx, lby);
			glEnd();
		}
}

void Renderer :: renderMat(matrix<double> mat,int mode)
{
	double gridLBX = this->gridLBX;
	double gridLBY = this->gridLBY;
	
	int fillFlag = 0;
	
	switch(mode)
	{
		case 1:
	
		for (unsigned int i=0;i<mat.size1();i++){
			for (unsigned int j=0;j<mat.size2();j++){
				glColor3f(1,0,0);//set fill color..used if fillFlag=true
				if(mat(i,j)!=0)
					fillFlag = 1;
				else
					fillFlag = 0;
				
				//drawSquare(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,fillFlag);//1: Filled grid,0:wire Frame
				drawSquare(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,(int)mat(i,j));//1: Filled grid,0:wire Frame
				gridLBX+=stepX;
			}
			gridLBY += stepY;
			gridLBX = this->gridLBX;
		}
		break;
		case 2 :

	//	int max  = 6;
		for (unsigned int i=0;i<mat.size1();i++){
			for (unsigned int j=0;j<mat.size2();j++){
				double fillVal = -2.0;
				if(i==0||j==0||i==mat.size1()-1||j==mat.size2()-1)
				{
					//glColor3f(0.2,0.2,0.2);//set fill color..used if fillFlag=true
					//drawSquareFilled(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,fillVal);//1: following airs color
				}
				else{
					fillVal  = mat(i,j);//(double)abs(mat(i,j))/max ;
					drawSquareFilled(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,fillVal);//1: Filled grid,0:wire Frame
						//drawSquare(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,(int)mat(i,j));//1: Filled grid,0:wire Frame

				}
				gridLBX+=stepX;
			}
			gridLBY += stepY;
			gridLBX = this->gridLBX;
		}
		break;
				
				
		
	}
	renderBoundary();
}


void Renderer :: renderVector2D(matrix<double> u, matrix<double> v)
{
		//float pi = 22/7;
	#define pi 22/7
	
	double gridLBX = this->gridLBX;
	double gridLBY = this->gridLBY;
	double stepX = this->stepX;
	double stepY = this->stepY;
	int nX = (this->sGrid)->nX;
	int nY = (this->sGrid)->nY;
		
	matrix<double> vel_ang(nY,nX);	
	matrix<double> vel_mag(nY,nX);
		
	for(int i=0;i< nX;i++)
		for(int j=0;j< nY;j++)
		{


//if (V(i,j) != 0.0f && U(i,j) != 0.0f) {
			vel_ang(i,j) = (float)atan((float)U(i,j)/(V(i,j)))*180/(pi);
			vel_mag(i,j) = sqrt(U(i,j)*U(i,j)+V(i,j)*V(i,j));
	//		}
		}
		 
	
	int flag_Quad=0;

	static double maxv=0 ;

	maxv = MAX(maxv,getMax(vel_mag));
	//maxv = getMax(vel_mag) ;
	
	for (int i=0;i<nX;i++)
	{
		for (int j=0;j<nY;j++)
		{
			if(i==0+sb||j==0+sb||i==nX-1-sb||j==nY-1-sb)
			{
				glColor3f(0.3,0.3,0.3);
				//drawSquare(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,flag_Quad=4);//1: following airs color
			}	
			else
			{
				glColor3f(1,0,0);
				flag_Quad=0;
				//drawSquare(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,flag_Quad);//1: on the quad flag
				drawVectorAngnMag(gridLBX,gridLBY,gridLBX+stepX,gridLBY+stepY,SIGN(U(i,j)),SIGN(V(i,j)),vel_mag(i,j),maxv,vel_ang(i,j));
			}
			
			gridLBX+=stepX;
		}
		gridLBY+=stepY;
		gridLBX = this->gridLBX;
	}
	//renderBoundary();
}



void Renderer :: drawVectorAngnMag(double lbx,double lby,double rtx,double rty,bool usign,bool vsign,double mag,double max,double ang)
{//sign=0  +ve, =1  -ve
/*
 * u  v  actual_ang   req_ang    atan(...)
 * +  -    -ve          +ve         u/v
 * -  -    +ve          -ve         u/v
 * +  +    +ve         180-ang      u/v
 * -  +    -ve         180+ang      u/v
 */
	//cout<<"\nang:"<<ang<<endl;
	ang*=(-1);
	if(!usign && !vsign)
		ang+=180;
	if(usign && !vsign)
		ang-=180;
	//cout<<"\nang:"<<ang<<endl;
	{
		float ht = 0.9*(rty-lby) /max *mag;
		
		glPushMatrix();
		glTranslatef ((lbx+rtx)/2, (rty+lby)/2, 0.0);
		glRotatef(ang,0.0,0.0,1.0);
		glTranslatef (-(lbx+rtx)/2, -((rty+lby)/2), 0.0);
		
		//cout<<" "<<max<<" "<<mag<<endl;
		glColor3f(1.0,0,0);
		glBegin(GL_LINES);		
		glVertex2f((lbx+rtx)/2, (rty+lby)/2+ht/2);
		glVertex2f((lbx+rtx)/2, (rty+lby)/2-ht/2);
		
		glVertex2f((lbx+rtx)/2, (rty+lby)/2-ht/2);
		glVertex2f((lbx+rtx)/2+ht/10, (rty+lby)/2-ht/2+ht/10);
		
		glVertex2f((lbx+rtx)/2, (rty+lby)/2-ht/2);
		glVertex2f((lbx+rtx)/2-ht/10, (rty+lby)/2-ht/2+ht/10);
		glEnd();
		glPopMatrix();
	}	
}


double Renderer :: getMax(matrix<double> mat) //get max value from the 2D array
{
	double max = 0.0;
	for(unsigned int i=1;i< mat.size1()-1;i++)
		for(unsigned int j=1;j<mat.size2()-1;j++)
			if(max<mat(i,j))
			  max=mat(i,j);
	return max;
}



void Renderer :: renderDensity(matrix<double> mat)	
{
	int d00,d01,d10,d11;	
	double gridLBX = this->gridLBX;
	double gridLBY = this->gridLBY;
	for (unsigned int i=0;i< mat.size1()-1;i++)
	{
		for (unsigned int j=0;j < mat.size2()-1;j++)
		{
			glBegin ( GL_QUADS );
				d00 = mat(i,j);
				d01 = mat(i,j+1);
				d10 = mat(i+1,j);
				d11 = mat(i+1,j+1);

				glColor3f ( d00, d00, d00 ); glVertex2f ( gridLBX, gridLBY );
				glColor3f ( d10, d10, d10 ); glVertex2f ( gridLBX, gridLBY+stepY );
				glColor3f ( d11, d11, d11 ); glVertex2f ( gridLBX+stepX, gridLBY+stepY );
				glColor3f ( d01, d01, d01 ); glVertex2f (  gridLBX+stepX, gridLBY );
			gridLBX+=stepX;
		glEnd();
		}
		gridLBY+=stepY;
		gridLBX = this->gridLBX;
	}
}


/*---for staggered with stam---*/
	
void Renderer :: renderVel2D_Stam()
{
	double h;
	//h = 1.0f;// /sGrid->nX ;
	h = sGrid->dx;
	double x,y ;
	matrix <double> u = sGrid->u; 
	matrix <double> v = sGrid->v; 
	glColor3f( 1.0f, 1.0f, 1.0f);
	glLineWidth ( 1.0f );
	glBegin (GL_LINES) ;
		for( int i=1; i< sGrid->nY; i++){  //not showing boundary..
			y = (i+0.5f)*h ;
	
			for( int j=1; j< sGrid->nX; j++){
				x = (j+0.5f)*h;

				glVertex2f( x, y );
				glVertex2f( x+U(i,j),  y+V(i,j) ) ;
				
			}
		}
	glEnd();
}

void Renderer :: renderDen2D_Stam()
{
	int i, j;
	double x, y, h, d00, d01, d10, d11;
	
	//h = 1.0f;// /sGrid->nX;
	h = sGrid->dx;

	glBegin ( GL_QUADS );
/* here last grid row and col is not rendered..*/
		for ( i=0 ; i< sGrid->nY - 1 ; i++ ) { //ignored the right/top boundary
			y = (i+0.5f)*h;
			for ( j=0 ; j < sGrid->nX - 1 ; j++ ) {
				x = (j+0.5f)*h; // this 0.5 matter is difficult to understand..if removed smoothness reduces..
				d00 = sGrid->d(i,j);
				d01 = sGrid->d(i,j+1);
				d10 = sGrid->d(i+1,j);
				d11 = sGrid->d(i+1,j+1);

				glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
				glColor3f ( d10, d10, d10 ); glVertex2f ( x, y+h );
				glColor3f ( d11, d11, d11 ); glVertex2f ( x+h, y+h );
				glColor3f ( d01, d01, d01 ); glVertex2f (  x+h, y );
			}
		}
	glEnd ();
}


void Renderer :: renderParticles()
{
//cout<<"renderPar"<<endl;
	glColor3f(0.3,0.3,0.3);
    //renderGrid();
    //glColor3f(0.2,0.2,0.9);
    //renderBoundary();
 	for (unsigned i = 0; i < sGrid->fluidParticles.size() ; i++ ){
	glPointSize(4);
	glBegin(GL_POINTS);
		//glColor3f(1,i*0.001,0);
		glColor3f(0,0,0.6);
		glVertex2d(sGrid->fluidParticles.at(i)->x,sGrid->fluidParticles.at(i)->y);
	glEnd();	
	}
}

void Renderer :: renderParticle(double x, double y,double R,double G,double B,int psize) 
{
	
 	glPointSize(psize);
	glBegin(GL_POINTS);
	glColor3f(R,G,B);
		glVertex2d(x,y);
	glEnd();	
}

void Renderer :: renderSurfaceBoundary() {
    for (int i = sGrid->nX - 1; i > 0; i--)
        for (int j = 0; j < sGrid->nY-1; j++) {
        	int state=0;
        	bool UL = sGrid->cellType(i, j) == FLUID,
        		 UR = sGrid->cellType(i, j+1) == FLUID,
        		 BL = sGrid->cellType(i-1, j) == FLUID,
        		 BR = sGrid->cellType(i-1 ,j+1) == FLUID ;

                if (UL)
                    state = state | 1;
                if (UR)
                    state = state | 2;
                if (BL)
                    state = state | 4;
                if (BR)
                    state = state | 8;

                if(!state)
                	continue;

                double x = gridLBX + (j) * stepX;
                double y = gridLBY + (i) * stepY;
                double incfactor = stepX / 2.0;
                x+=incfactor;
                y+=incfactor;

                glColor3f(1.0, 0.0, 0.0);
                glBegin(GL_LINES);
                switch (state) {
                    case 1:
                        glVertex2f(x, y - incfactor);
                        glVertex2f(x + incfactor, y);
                        break;

                    case 2:
                        glVertex2f(x + incfactor, y);
                        glVertex2f(x + 2 * incfactor, y - incfactor);
                        break;

                    case 3:
                        glVertex2f(x, y - incfactor);
                        glVertex2f(x + 2 * incfactor, y - incfactor);
                        break;

                    case 4:
                        glVertex2f(x, y - incfactor);
                        glVertex2f(x + incfactor, y - 2 * incfactor);
                        break;

                    case 5:
                        glVertex2f(x + incfactor, y);
                        glVertex2f(x + incfactor, y - 2 * incfactor);
                        break;

                    case 6:
                        glVertex2f(x, y - incfactor);
                        glVertex2f(x + incfactor, y);
                        glVertex2f(x + incfactor, y - 2 * incfactor);
                        glVertex2f(x + 2 * incfactor, y - incfactor);
                        break;

                    case 7:
                        glVertex2f(x + incfactor, y - 2 * incfactor);
                        glVertex2f(x + 2 * incfactor, y - incfactor);
                        break;

                    case 8:
                        glVertex2f(x + incfactor, y - 2 * incfactor);
                        glVertex2f(x + 2 * incfactor, y - incfactor);
                        break;

                    case 9:
                        glVertex2f(x + incfactor, y);
                        glVertex2f(x + 2 * incfactor, y - incfactor);
                        glVertex2f(x, y - incfactor);
                        glVertex2f(x + incfactor, y - 2 * incfactor);
                        break;

                    case 10:
                        glVertex2f(x + incfactor, y);
                        glVertex2f(x + incfactor, y - 2 * incfactor);
                        break;

                    case 11:
                        glVertex2f(x, y - incfactor);
                        glVertex2f(x + incfactor, y - 2 * incfactor);
                        break;

                    case 12:
                        glVertex2f(x, y - incfactor);
                        glVertex2f(x + 2 * incfactor, y - incfactor);
                        break;

                    case 13:
                        glVertex2f(x + incfactor, y);
                        glVertex2f(x + 2 * incfactor, y - incfactor);
                        break;

                    case 14:
                        glVertex2f(x, y - incfactor);
                        glVertex2f(x + incfactor, y);
                        break;
                }
              glEnd();

        }
}


