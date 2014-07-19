#include "FluidSim.h"
#include "opengl.h"
#include "pcg.h"
#include "Printer.h"
#include "Vec.h"
#include "pcgsolver/util.h"
#include <cassert>
using namespace std;
extern int sb;
extern int ProjectFLAG;
extern bool AdvectFlag ;

double FluidSim :: getMinDistance(double a, double b, double r){
	double d ;

	if (fabs(a-b) < 1 )
		d = (a + b + sqrt(2-(a-b)))/2.0;
	else
		d = (a>b?b:a) + 1;

	return(d<r?d:r);
}

void FluidSim :: calculateLevelSetDistance(){
	for(int y = 0;y < sGrid->nX; y++)
		for(int x = 0;x < sGrid->nY; x++){
			if(sGrid->cellType(y, x) == FLUID){
				sGrid->distanceLevelSet(y,x) = -1;
			}
			else
				sGrid->distanceLevelSet(y, x) = sGrid->nX*2+10;
		}

	for(int l=0; l<2;l++){

		for(int y = 1;y < sGrid->nX; y++)
			for(int x = 1;x < sGrid->nY; x++){
				if(sGrid->cellType(y, x) != FLUID){
					sGrid->distanceLevelSet(y,x) = getMinDistance(
												sGrid->distanceLevelSet(y,x-1),
												sGrid->distanceLevelSet(y-1,x),
												sGrid->distanceLevelSet(y,x)
												);
				}
			}
		for(int y = sGrid->nX-2;y >=0; y--)
			for(int x = 1;x < sGrid->nY; x++){
				if(sGrid->cellType(y, x) !=FLUID){
					sGrid->distanceLevelSet(y,x) = getMinDistance(
												sGrid->distanceLevelSet(y,x-1),
												sGrid->distanceLevelSet(y+1,x),
												sGrid->distanceLevelSet(y,x)
												);
				}
			}
		for(int y = 1;y < sGrid->nX; y++)
			for(int x = sGrid->nY-2;x>=0; x--){
				if(sGrid->cellType(y, x) !=FLUID){
					sGrid->distanceLevelSet(y,x) = getMinDistance(
												sGrid->distanceLevelSet(y,x+1),
												sGrid->distanceLevelSet(y-1,x),
												sGrid->distanceLevelSet(y,x)
												);
				}
			}
		for(int y = sGrid->nX-2;y >=0; y--)
			for(int x = sGrid->nY-2;x>=0; x--){
				if(sGrid->cellType(y, x) !=FLUID){
					sGrid->distanceLevelSet(y,x) = getMinDistance(
												sGrid->distanceLevelSet(y,x+1),
												sGrid->distanceLevelSet(y+1,x),
												sGrid->distanceLevelSet(y,x)
												);
				}
			}
	}
}

void  FluidSim :: setVelocityError(matrix<double> &mat){
	forEach(i,0,(int)mat.size1()-1)
		forEach(j,0,(int)mat.size2()-1)
			if(fabs( mat(i,j) ) < 0.001)
				mat(i,j)=0;
}

void  FluidSim :: setValidVelocity(int val)
{
	uValid.clear();
	vValid.clear();
	for(int y=1; y < sGrid->nY-1;y++)
		for(int x=1; x < sGrid->nX-1;x++){
			if(sGrid->cellType(y,x) == FLUID){
				uValid(y,x)=val ;
				vValid(y,x)=val ;
				uValid(y,x+1)=val ;
				vValid(y+1,x)=val ;
			}
		}
}

void  FluidSim :: simulate(double timestep)
{
	dt = cfl() > 0 && cfl()< 1 ? cfl():0.1f ;
	{   //mark fluid cells
		advectParticles(sGrid->fluidParticles, sGrid->u,sGrid->v, dt);
		markFluidCells();
			matrix<double > u = sGrid->u;
			matrix<double > v = sGrid->v;
        //advect velocity
		calculateLevelSetDistance();
		sGrid->u = advect2DSelf(sGrid->u,dt,u,v,1);
		sGrid->v = advect2DSelf(sGrid->v,dt,u,v,2);
		applyBoundaryConditions(VELOCITY_BC2);
		//add Gravity
		sGrid->v = addGravity(sGrid->v,dt);
		applyBoundaryConditions(VELOCITY_BC2);
		//add Viscosity
//		addViscosity(0,dt);
//		applyBoundaryConditions(VELOCITY_BC2);
		//apply Pressure
		switch(ProjectFLAG){
		case GSRflag :
			//sGrid->p = projectGSr(sGrid->u, sGrid->v, dt); //used u and v and result stored again in u ,v
			break;
		case PCGflag:
			//projectPCG();
			break;
		case PCGBridsonflag:
			solvePressureBridson((float)dt);
			break;
		case NoPressure:
			break;
		}
		applyBoundaryConditions(VELOCITY_BC2);
        //extrapolation
		extern bool extrapolateFlag;
		if(extrapolateFlag){
			setValidVelocity(1);
			extrapolate2D(sGrid->u,uValid);
			extrapolate2D(sGrid->v,vValid);
			setValidVelocity(0);
		}
		applyBoundaryConditions(VELOCITY_BC2);
        //advect Particles
		//advectParticles(sGrid->fluidParticles, sGrid->u,sGrid->v, dt);
		//markFluidCells();
	}
}


void FluidSim :: applyBoundaryConditions(int bc)//boundary Condition
{
	switch(bc)
	{
	case VELOCITY_BC :
		//bottom boundary
		for(int x = 0+sb; x < sGrid->nX-sb; x++){
			sGrid->v(0+sb,x) =  0;//-sGrid->v(2+sb,x);
			//sGrid->v(1+sb,x) = 0;//sGrid->v(2+sb,x);
			//sGrid->d(0,x) = 0;
		}
		//Top boundary..
		for(int x = 0+sb; x < sGrid->nX-sb; x++){  //allowing in x direction
			sGrid->v(sGrid->nY-1-sb,x) = 0;//sGrid->v(sGrid->nY-2-sb,x);
			//	sGrid->v(sGrid->nY-sb,x) = 0;//-sGrid->v(sGrid->nY-2-sb,x);
			//sGrid->d(sGrid->nY,x) = 0;
		}
		for(int y = 0+sb; y < sGrid->nY-sb; y++){
			//	sGrid->u(y,0+sb) = 0;//-sGrid->u(y,2+sb);
			sGrid->u(y,1+sb) = 0;//sGrid->u(y,2+sb);
			//sGrid->d(y,0) = 0;
		}
		//Right boundary..
		for(int y = 0+sb; y < sGrid->nY-sb; y++){
			sGrid->u(y,sGrid->nX-1-sb) = 0;//-sGrid->u(y,sGrid->nX-2-sb);
			//	sGrid->u(y,sGrid->nX-sb) = 0;//sGrid->u(y,sGrid->nX-2-sb);
			//sGrid->d(y,sGrid->nX) = 0;
		}
		break;
	case VELOCITY_BC2 :
		//bottom boundary
		for(int x = 0+sb; x < sGrid->nX-sb; x++){
			sGrid->v(0+sb,x) =  0;//-sGrid->v(2+sb,x);
			sGrid->v(1+sb,x) = 0;//sGrid->v(2+sb,x);
		}
		//Top boundary..
		for(int x = 0+sb; x < sGrid->nX-sb; x++){  //allowing in x direction
			sGrid->v(sGrid->nY-1-sb,x) = 0;//sGrid->v(sGrid->nY-2-sb,x);
			sGrid->v(sGrid->nY-sb,x) = 0;//-sGrid->v(sGrid->nY-2-sb,x);
		}
		for(int y = 0+sb; y < sGrid->nY-sb; y++){
			sGrid->u(y,0+sb) = 0;//-sGrid->u(y,2+sb);
			sGrid->u(y,1+sb) = 0;//sGrid->u(y,2+sb);
		}
		//Right boundary..
		for(int y = 0+sb; y < sGrid->nY-sb; y++){
			sGrid->u(y,sGrid->nX-1-sb) = 0;//-sGrid->u(y,sGrid->nX-2-sb);
			sGrid->u(y,sGrid->nX-sb) = 0;//sGrid->u(y,sGrid->nX-2-sb);
		}
		break;

	case DENSITY_BC :
		//bottom boundary
		for(int x = 0+sb; x < sGrid->nX-sb; x++){ sGrid->d(0+sb,x) = 0; }
		//Top boundary..
		for(int x = 0+sb; x < sGrid->nX-sb; x++){ sGrid->d(sGrid->nY-1-sb,x) = 0; }
		for(int y = 0+sb; y < sGrid->nY-sb; y++){ sGrid->d(y,0+sb) = 0; }
		//Right boundary..
		for(int y = 0+sb; y < sGrid->nY-sb; y++){ sGrid->d(y,sGrid->nX-1-sb) = 0; }
		break;
	case PRESSURE_BC :
		//bottom boundary
		for(int x = 0+sb; x < sGrid->nX-sb; x++){ sGrid->p(0+sb,x) = 0; }
		//Top boundary..
		for(int x = 0+sb; x < sGrid->nX-sb; x++){ sGrid->p(sGrid->nY-1-sb,x) = 0; }
		for(int y = 0+sb; y < sGrid->nY-sb; y++){ sGrid->p(y,0+sb) = 0; }
		//Right boundary..
		for(int y = 0+sb; y < sGrid->nY-sb; y++){ sGrid->p(y,sGrid->nX-1-sb) = 0; }
		break;
	}
}

void FluidSim :: RK2(double &posx, double &posy,matrix<double> u, matrix<double> v, double dt){
	double dx = sGrid->dx;
	double x = posx*dx;
	double y = posy*dx;
	int nX = sGrid->nX;
	int nY = sGrid->nY;
	double velx = getVelInterpolated(posx,posy-0.5,u);
	double vely = getVelInterpolated(posx-0.5,posy,v);
	posx = x  + 0.5*dt*velx;
	posy = y  + 0.5*dt*vely;

	posx = posx < (1+sb)*dx ? (1+sb)*dx:posx;
	posy = posy < (1+sb)*dx ? (1+sb)*dx:posy;
	posx = posx > (nX-1-sb)*dx ? (nX-1-sb)*dx:posx;
	posy = posy > (nY-1-sb)*dx ? (nY-1-sb)*dx:posy;

	posx/=dx;
	posy/=dx;
	velx = getVelInterpolated(posx,posy-0.5,u);
	vely = getVelInterpolated(posx-0.5,posy,v);
	posx = x  + dt*velx;
	posy = y  + dt*vely;

	posx = posx < (1+sb)*dx ? (1+sb)*dx:posx;
	posy = posy < (1+sb)*dx ? (1+sb)*dx:posy;
	posx = posx > (nX-1-sb)*dx ? (nX-1-sb)*dx:posx;
	posy = posy > (nY-1-sb)*dx ? (nY-1-sb)*dx:posy;

	posx/=dx;
	posy/=dx;
}

matrix<double> FluidSim :: advect2DSelf(matrix<double> q, double dt, matrix<double> u, matrix<double> v,int component)
{

	if(!AdvectFlag){
		matrix<double> temp=q;
		temp.clear();

		int nX = this->sGrid->nX;
		int nY = this->sGrid->nY;
		double dx = sGrid->dx;
		double x, y ;
		for(int i=0+sb;i<=nY-2-sb;i++)
			for(int j=0+sb;j<=nX-2-sb;j++)
			{
				//here advection is done using the cell size...
				//result is ok....
				//not working with considering particle at center i.e. at (i+0.5) and (j+0.5)...
				y = (i)*dx - V(i,j)*dt;
				x = (j)*dx - U(i,j)*dt;
				//cout<<"\n.."<<x<<" "<<y<<" "<<(x<=(0.5+sb)*dx)<<endl;
				if(x<=(0.5+sb)*dx) x=(0.5+sb)*dx;
				if(y<=(0.5+sb)*dx) y=(0.5+sb)*dx;
				if(x>=(nX-2+0.5-sb)*dx) x= (nX-2+0.5-sb)*dx;
				if(y>=(nY-2+0.5-sb)*dx) y= (nY-2+0.5-sb)*dx;
				//cout<<"\n.."<<x<<" "<<y<<" "<<(x<=(0.5+sb)*dx)<<" "<<(0.5+sb)*dx<<endl;

				int xl = (int)(x/dx);
				int yl = (int)(y/dx);
				int xh = xl + 1;
				int yh = yl + 1;
				//cout<<"\n.."<<xl<<" "<<yl<<" "<<xh<<" "<<yh<<endl;

				double alpha = x/dx-xl;
				double beta = y/dx-yl;
				//		temp(i,j) = q((int)y,(int)x);//getting velocity from exact that point won't produce anything...
				//need to interpolate from neighbours...
				//cout<<"\n.."<<xl<<" "<<yl<<" "<<xh<<" "<<yh<<" "<<alpha<<" "<<beta<<endl;
				assert(xl>=0 || yl>=0 || xh<=nX-1 || yh<=nY-1);
				temp(i,j) = (1-alpha)*( (1-beta)*q(yl,xl) + beta*q(yh,xl) ) +\
						alpha *( (1-beta)*q(yl,xh) + beta*q(yh,xh) );
			}
		return temp;
	}
	else{ //proper advection - RK2
		matrix<double> temp=q;
		temp.clear();
		int nX = this->sGrid->nX;
		int nY = this->sGrid->nY;
		double dx = sGrid->dx;
		//dt*=nX;
		double x,y, posx, posy;
		if(component==1){
			for(int i=1+sb;i<=nY-2-sb;i++)
				for(int j=1+sb;j<=nX-1-sb;j++){
					x = (j)*dx;
					y = (i+0.5)*dx;
					posx = x/dx;
					posy = y/dx;
					RK2(posx,posy,sGrid->u,sGrid->v,-dt);
					temp(i,j) = getVelInterpolated(posx,posy-0.5,sGrid->u);
				}
		}
		else{ //component=2 i.e. Vertical component
			for(int i=1+sb;i<=nY-1-sb;i++)
				for(int j=1+sb;j<=nX-2-sb;j++){
					x = (j+0.5)*dx;
					y = (i)*dx;
					posx = x/dx;
					posy = y/dx;
					RK2(posx,posy,sGrid->u,sGrid->v,-dt);
					temp(i,j) = getVelInterpolated(posx-0.5,posy,sGrid->v);
				}
		}
		return temp;
	}
}

#include <vector>
//extern std::vector<int> traceAdvectCells;
extern std::vector<double> traceAdvectCellsd;
void FluidSim :: advect2DCell(int j,int i, double dt, matrix<double> u, matrix<double> v,int component=2)
{
	extern int traceCellFlag;
	traceCellFlag = ON;
	double dx = sGrid->dx;
	//dt*=nX;
	double x,y, posx, posy;
	if(component==1){
		x = (j)*dx;
		y = (i+0.5)*dx;
		posx = x/dx;
		posy = y/dx;
		RK2(posx,posy,sGrid->u,sGrid->v,-dt);
	}
	else{

		x = (j+0.5)*dx;
		y = (i)*dx;
		posx = x/dx;
		posy = y/dx;
		RK2(posx,posy,sGrid->u,sGrid->v,-dt);
	}

	cout <<"Advected From : "<<posx<<" "<<posy<<endl;
	traceAdvectCellsd[0]=x;
	traceAdvectCellsd[1]=y;
	traceAdvectCellsd[2]=posx*dx;
	traceAdvectCellsd[3]=posy*dx;
	glutPostRedisplay();
}

void FluidSim :: advectParticles(std::vector <Particles*> & plist, matrix<double> u, matrix<double>v, double dt)
{

	double dx = sGrid->dx;
	double x, y ;
	extern bool AdvectFlagP;
	if(!AdvectFlagP){
		for ( unsigned i = 0; i < plist.size(); i++){
			int celli = plist.at(i)->y /dx ;
			int cellj = plist.at(i)->x /dx ;

			y = plist.at(i)->y + v(celli,cellj)*dt*0.5; //newPos = oldPos + velocity*dt;
			x = plist.at(i)->x + u(celli,cellj)*dt*0.5;

			celli = y /dx ;
			cellj = x /dx ;

			y = plist.at(i)->y + v(celli,cellj)*dt; //newPos = oldPos + velocity*dt;
			x = plist.at(i)->x + u	(celli,cellj)*dt;

			plist.at(i)->x = x;
			plist.at(i)->y = y;
		}
	}else{
		for ( unsigned i = 0; i < plist.size(); i++){
			double posx = plist.at(i)->x /dx ;
			double posy = plist.at(i)->y /dx ;

			RK2(posx, posy, sGrid->u, sGrid->v, dt);

			plist.at(i)->x = posx*dx;
			plist.at(i)->y = posy*dx;
		}
	}
}

void FluidSim :: advect2DParticle(double x, double y, double dt, matrix<double> u, matrix<double> v)
{

	int nX = this->sGrid->nX;
	int nY = this->sGrid->nY;
	double dx = this->sGrid->dx;

	x*=dx;
	y*=dx;
	double posx = x /dx ;
	double posy = y /dx;

	double velx = getVelInterpolated(posx,posy-0.5,sGrid->u);
	double vely = getVelInterpolated(posx-0.5,posy,sGrid->v);

	posx = x  + 0.5*dt*velx;
	posy = y  + 0.5*dt*vely;

	posx = posx < (1+sb)*dx ? (1+sb)*dx:posx;
	posy = posy < (1+sb)*dx ? (1+sb)*dx:posy;
	posx = posx > (nX-1-sb)*dx ? (nX-1-sb)*dx:posx;
	posy = posy > (nY-1-sb)*dx ? (nY-1-sb)*dx:posy;

	posx/=dx;
	posy/=dx;
	velx = getVelInterpolated(posx,posy-0.5,sGrid->u);
	vely = getVelInterpolated(posx-0.5,posy,sGrid->v);

	posx = x  + dt*velx;
	posy = y  + dt*vely;
	posx = posx < (1+sb)*dx ? (1+sb)*dx:posx;
	posy = posy < (1+sb)*dx ? (1+sb)*dx:posy;
	posx = posx > (nX-1-sb)*dx ? (nX-1-sb)*dx:posx;
	posy = posy > (nY-1-sb)*dx ? (nY-1-sb)*dx:posy;

	posx/=dx;
	posy/=dx;

	cout <<"Advected From : "<<posx<<" "<<posy<<endl;
	traceAdvectCellsd[0]=x;
	traceAdvectCellsd[1]=y;
	traceAdvectCellsd[2]=posx*dx;
	traceAdvectCellsd[3]=posy*dx;
	glutPostRedisplay();

}

void FluidSim :: init(GridStag* sGrid)
{
	this->sGrid = sGrid;
	uValid.resize(sGrid->nY,sGrid->nX+1);
	vValid.resize(sGrid->nY+1,sGrid->nX);
	uValid.clear();
	vValid.clear();
}

double FluidSim ::  cfl()
{
	double maxVel = 0.0;
	for(int y = 0+sb;y < sGrid->nX-sb; y++)
		for(int x = 0+sb;x < sGrid->nY+1-sb; x++){
			if(maxVel < fabs(sGrid->u(y,x)))
				maxVel = fabs(sGrid->u(y,x));
		}
	for(int y = 0+sb;y < sGrid->nX+1-sb; y++)
		for(int x = 0+sb;x < sGrid->nY-sb; x++){
			if(maxVel < fabs(sGrid->v(y,x)))
				maxVel = fabs(sGrid->v(y,x));
		}
	if ( !maxVel )
		return 0;
	return (sGrid->dx/maxVel);
}

Vec2d FluidSim :: get_Velocity(const Vec2d& position) {
	int i, j;
	double fx=0, fy=0;
	double dx = sGrid->dx;
	int nX = this->sGrid->nX;
	int nY = this->sGrid->nY;
	matrix<double> u = sGrid->u;
	matrix<double> v = sGrid->v;

	Vec2d posX = (position / dx - Vec2d(0.0, 0.5));
	Vec2d posY = (position / dx - Vec2d(0.5, 0.0));

	// cout<<position[0]/dx<<" "<<position[1]/dx<<" "<<fx<<" "<<fy<<endl; ;
	//double u_value = getVelInterpolated(posX[0],posX[1],u);
	//double v_value = getVelInterpolated(posY[0],posY[1],v);

	get_barycentric(posX[0], i, fx, 0+sb, nX-sb);
	get_barycentric(posX[1], j, fy, 0+sb, nY-sb);
	//cout<<posX[0]<<" "<<posX[1]<<" "<<fx<<" "<<fy<<endl;

	double u_value = bilerp(
			u(i, j), u(i, j + 1),
			u(i + 1, j), u(i + 1, j + 1),
			fx, fy);

	get_barycentric(posY[0], i, fx, 0+sb, nX-sb);
	get_barycentric(posY[1], j, fy, 0+sb, nY-sb);

	double v_value = bilerp(
			v(i, j), v(i , j + 1 ),
			v(i + 1, j ), v(i + 1, j + 1),
			fx, fy);
	return Vec2d(u_value, v_value);
}

Vec2d FluidSim :: trace_rk2(const Vec2d& position, double dt) {
	Vec2d newPosition = position ;
	//cout<<position[0]<<" "<<newPosition[0]<endl;
	Vec2d velocity = get_Velocity(newPosition);
	velocity = get_Velocity(newPosition + 0.5f * dt * velocity);
	newPosition += dt*velocity;
	return newPosition;
}

matrix<double> FluidSim :: addForce(matrix<double> dest, double dt, matrix<double> src)
{
	for(int i=1+sb;i < sGrid->nY-1-sb;i++) //exclude the boundary cells
		for(int j=1+sb;j<sGrid->nX-1-sb;j++)
			dest( i, j ) = dest( i, j ) + dt* src(i,j);
	return dest;
}


double FluidSim :: getVelInterpolated(double x,double y, matrix<double> mat)
{
	int i = (int)floor(x);
	int j = (int)floor(y);
	return (i+1-x) * (j+1-y) * mat(j,i)+
			(x-i) * (j+1-y) * mat(j,i+1)+
			(i+1-x) * (y-j) * mat(j+1,i)+
			(x-i) * (y-j) * mat(j+1,i+1);
}

void FluidSim :: initFluidBody(int FluidPos)
{
	int nX = this->sGrid->nX;
	int nY = this->sGrid->nY;
	double dx = this->sGrid->dx;
	int fluidPosBX,
	fluidPosBY,
	fluidPosTX,
	fluidPosTY;
	extern double zoomFactor;
	switch(FluidPos)
	{
	case 1: //dam_break- left right corner
		fluidPosBX = (nX-1-sb)-(nX-1-2*sb)/3;  //last one is boundary
		fluidPosBY = 1 ;//sb+1;
		fluidPosTX = nX-1-sb;
		fluidPosTY = ((nY-1-sb) - (nY-1-2*sb)/3) ;
		for (int i = 0+sb; i < (nX-sb)*(nY-sb)*16; ++i) {
			float xpos = randhashf(i * 2, 0.0+sb/nX, zoomFactor-sb/nX);
			float ypos = randhashf(i * 2 + 1,0.0+sb/nX, zoomFactor-sb/nX);

			if (xpos > fluidPosBX*dx && xpos < fluidPosTX*dx
					&& ypos > fluidPosBY*dx && ypos < fluidPosTY*dx ) {
				Particles* p = new Particles( xpos, ypos );
				if (p != NULL)
					sGrid->fluidParticles.push_back(p);
				else
					cout << "Error...Not Pushed back Particle";
			}
		}
		break;
	case 2: //dam_break- middle-middle
		fluidPosBY = 0+3*(nX-sb)/10;  //last one is boundary
		fluidPosBX = 0+3*(nY-sb)/10;
		fluidPosTY = fluidPosBY + 6*(nX-sb)/10;
		fluidPosTX = fluidPosBX + 4*(nY-sb)/10;

		for (int i = 0+sb; i < (nX-sb)*(nY-sb)*16; ++i) {
			float xpos = randhashf(i * 2, 0.0+sb/nX, zoomFactor-sb/nX);
			float ypos = randhashf(i * 2 + 1, 0.0+sb/nX, zoomFactor-sb/nX);


			if (xpos > fluidPosBX*dx && xpos < fluidPosTX*dx
					&& ypos > fluidPosBY*dx && ypos < fluidPosTY*dx ) {
				Particles* p = new Particles( xpos, ypos );
				if (p != NULL)
					sGrid->fluidParticles.push_back(p);
				else
					cout << "Error...Not Pushed back Particle";
			}
		}
		break;
	case 3: //Static Bed
		fluidPosBY = 0+1*(nX-sb)/20;  //last one is boundary
		fluidPosBX = 0+1*(nY-sb)/20;
		fluidPosTY = fluidPosBY + 6*(nX-sb)/20;
		fluidPosTX = 32 - fluidPosBX ;

		for (int i = 0+sb; i < (nX-sb)*(nY-sb)*16; ++i) {
			float xpos = randhashf(i * 2, 0.0+sb/nX, zoomFactor-sb/nX);
			float ypos = randhashf(i * 2 + 1, 0.0+sb/nX, zoomFactor-sb/nX);


			if (xpos > fluidPosBX*dx && xpos < fluidPosTX*dx
					&& ypos > fluidPosBY*dx && ypos < fluidPosTY*dx ) {
				Particles* p = new Particles( xpos, ypos );
				if (p != NULL)
					sGrid->fluidParticles.push_back(p);
				else
					cout << "Error...Not Pushed back Particle";
			}
		}
		break;

	case 4: //double Dam Break
		fluidPosBY = 0+1*(nX-sb)/20;  //last one is boundary
		fluidPosBX = 0+1*(nY-sb)/20;
		fluidPosTY = fluidPosBY + 12*(nX-sb)/20;
		fluidPosTX = fluidPosBX + 6*(nX-sb)/20 ;

		for (int i = 0+sb; i < (nX-sb)*(nY-sb)*16; ++i) {
			float xpos = randhashf(i * 2, 0.0+sb/nX, zoomFactor-sb/nX);
			float ypos = randhashf(i * 2 + 1, 0.0+sb/nX, zoomFactor-sb/nX);

			if (xpos > fluidPosBX*dx && xpos < fluidPosTX*dx
					&& ypos > fluidPosBY*dx && ypos < fluidPosTY*dx ) {
				Particles* p = new Particles( xpos, ypos );
				if (p != NULL)
					sGrid->fluidParticles.push_back(p);
				else
					cout << "Error...Not Pushed back Particle";
			}
		}
		for(int x = fluidPosBX; x<=fluidPosTX; x++)
			for(int y = fluidPosBY; y<fluidPosTY; y++){
				sGrid->cellType(y,x) = FLUID;
			}

		fluidPosBX = 0+19*(nX-sb)/20 - 5*(nX-sb)/20 ;  //last one is boundary
		fluidPosBY = 0+1*(nX-sb)/20;
		fluidPosTY = fluidPosBY + 12*(nY-sb)/20;
		fluidPosTX = 19*(nY-sb)/20 + 1;

		for (int i = 0+sb; i < (nX-sb)*(nY-sb)*16; ++i) {
			float xpos = randhashf(i * 2, 0.0+sb/nX, zoomFactor-sb/nX);
			float ypos = randhashf(i * 2 + 1, 0.0+sb/nX, zoomFactor-sb/nX);

			if (xpos > fluidPosBX*dx && xpos < fluidPosTX*dx
					&& ypos > fluidPosBY*dx && ypos < fluidPosTY*dx ) {
				Particles* p = new Particles( xpos, ypos );
				if (p != NULL)
					sGrid->fluidParticles.push_back(p);
				else
					cout << "Error...Not Pushed back Particle";
			}
		}
		break;
	}
	for(int x = fluidPosBX; x<=fluidPosTX; x++)
		for(int y = fluidPosBY; y<fluidPosTY; y++){
			sGrid->cellType(y,x) = FLUID;
		}
}

void FluidSim :: initFluidBody_Helper(int bx, int tx, int by, int ty, matrix<double>& mat, double val)
{
	for(int x = bx; x<=tx; x++)
		for(int y = by; y<ty; y++){
			mat(y,x) = val;
		}
}


matrix<double> FluidSim :: addGravity(matrix<double> ua, double dt)
{
	/* here we are not adding gravity factor to top line of grid i.e.(nY,0-nX) */
	matrix<double> ub = ua;
	setValidVelocity(1);
	for(int y=1; y < sGrid->nY-1; y++)
		for(int x=1; x < sGrid->nX-1; x++){
			if(vValid(y,x)){// && (fabs(ub(y,x)) > dt || fabs(ub(y,x))==0.0))
				ub(y,x) += dt*gravity;
			}
		}
	setValidVelocity(0);
	return ub;
}


void FluidSim :: initSolidBoundary(int choice)
{
	switch(choice)
	{
	case 1:
		//bottom boundary..
		for(int x = 0+sb; x < sGrid->nX-sb; x++){
			sGrid->cellType(0,x) = SOLID;
		}
		//Top boundary..
		for(int x = 0+sb; x < sGrid->nX-sb; x++){
			sGrid->cellType(sGrid->nY-1,x) = SOLID;
		}
		//Left boundary..
		for(int y = 0+sb; y < sGrid->nY-sb; y++){
			sGrid->cellType(y,0) = SOLID;
		}
		//Right boundary..
		for(int y = 0+sb; y < sGrid->nX-sb; y++){
			sGrid->cellType(y,sGrid->nX-1) = SOLID;
		}
		break;
	}
	this->applyBoundaryConditions(VELOCITY_BC2);
}


void FluidSim :: markFluidCells()
{
	double dx = sGrid->dx;

	sGrid->cellType.clear();
	sGrid->isFluidBoundary.clear();
	for (unsigned int i = 0; i < sGrid->fluidParticles.size(); i++) {
		int x = int(sGrid->fluidParticles.at(i)->x / dx);
		int y = int(sGrid->fluidParticles.at(i)->y / dx);
		sGrid->cellType(y, x) = FLUID; // means fluid cells
	}

	int eightNeighborCount = 0;
	int fourNeighborCount = 0;
    for (int y = 1; y < sGrid->nY - 1; y++)
        for (int x = 1; x < sGrid->nX - 1; x++) {
        	sGrid->isFluidBoundary(y, x) = 0;
        	if(sGrid->cellType(y,x) == FLUID){
        	eightNeighborCount =
                    int(sGrid->cellType(y + 1, x) == FLUID) +
                    int(sGrid->cellType(y - 1, x) == FLUID) +
                    int(sGrid->cellType(y, x - 1) == FLUID) +
                    int(sGrid->cellType(y, x + 1) == FLUID) +
                    int(sGrid->cellType(y + 1, x) == FLUID) +
                    int(sGrid->cellType(y - 1, x) == FLUID) +
                    int(sGrid->cellType(y, x - 1) == FLUID) +
                    int(sGrid->cellType(y, x + 1) == FLUID);
        	fourNeighborCount =
        			int(sGrid->cellType(y + 1, x + 1) == FLUID) +
                    int(sGrid->cellType(y - 1, x + 1) == FLUID) +
                    int(sGrid->cellType(y + 1, x - 1) == FLUID) +
                    int(sGrid->cellType(y - 1, x - 1) == FLUID);

        	if(eightNeighborCount>=7)
        		sGrid->cellType(y, x) = FLUID;
            if (fourNeighborCount < 4)
                sGrid->isFluidBoundary(y, x) = 1;
        	}
        }
	initSolidBoundary(1);
}


void FluidSim :: extrapolate2D(matrix<double> &grid, matrix<int> &valid)
{
	matrix<int> old_valid = valid;

	for(int layers = 0; layers < 2; ++layers) {
		old_valid = valid;
		matrix<double> temp_grid = grid;
		int nj = grid.size2()-1;
		int ni = grid.size1()-1;
		for( int j = 1; j <  (int)grid.size2()-1; ++j)
			for( int i = 1; i < (int)grid.size1()-1; ++i) {
			float sum = 0;
			int count = 0;

			if(!old_valid(i,j)) {

				if((i+1)<=ni && old_valid(i+1,j)) {
					sum += grid(i+1,j);\
					++count;
				}
				if((i-1)>=0 && old_valid(i-1,j)) {
					sum += grid(i-1,j);\
					++count;
				}
				if((j+1)<=nj && old_valid(i,j+1)) {
					sum += grid(i,j+1);\
					++count;
				}
				if((j-1)>=0 && old_valid(i,j-1)) {
					sum += grid(i,j-1);\
					++count;
				}

				//If any of neighbour cells were valid,
				//assign the cell their average value and tag it as valid
				if(count > 0) {
					temp_grid(i,j) = sum /(float)count;
					valid(i,j) = 1;
				}
			}
		}
		grid = temp_grid; //update with the new changes
	}
}

void FluidSim::solvePressureBridson(float dt) {

	unsigned int sys_size = sGrid->nX*sGrid->nY;
	if(rhs.size() != sys_size) {
		rhs.resize(sys_size);
		pressure.resize(sys_size);
		matrix1.resize(sys_size);
	}
	matrix1.zero();
	//Build the linear system for pressure
	double scale =  (dt / (double)sqr(sGrid->dx));
	matrix<double> cellType = sGrid->cellType;
	for(int j = 1; j < sGrid->nY-1; ++j) {
		for(int i = 1; i <  sGrid->nX-1; ++i) {
			int index = i + sGrid->nX*j;
			rhs[index] = 0;
			pressure[index] = 0;

			if(cellType(j,i)==FLUID)
			{
				int airCellCount = 0;
				if(cellType(j,i+1)==FLUID) {
					matrix1.add_to_element(index, index, scale);
					matrix1.add_to_element(index, index + 1, -scale);
				}
				else if(cellType(j,i+1)==AIR) {
					matrix1.add_to_element(index, index, scale);
					airCellCount++;
				}
				//left neighbour
				if(cellType(j,i-1)==FLUID) {
					matrix1.add_to_element(index, index, scale);
					matrix1.add_to_element(index, index - 1, -scale);
				}
				else if(cellType(j,i-1)==AIR) {
					matrix1.add_to_element(index, index, scale);
					airCellCount++;
				}

				//top neighbour
				if(cellType(j+1,i)==FLUID) {
					matrix1.add_to_element(index, index, scale);
					matrix1.add_to_element(index, index + sGrid->nX, -scale);
				}
				else  if(cellType(j+1,i)==AIR)
				{
					matrix1.add_to_element(index, index, scale);
					airCellCount++;
				}

				//bottom neighbour
				if(cellType(j-1,i)==FLUID) {
					matrix1.add_to_element(index, index, scale);
					matrix1.add_to_element(index, index - sGrid->nX, -scale);
				}
				else if(cellType(j-1,i)==AIR) {
					matrix1.add_to_element(index, index, scale);
					airCellCount++;
				}

				rhs[index]  = -((
						(sGrid->u(j,i+1) - sGrid->u(j,i) + sGrid->v(j+1,i) - sGrid->v(j,i))/(double)sGrid->dx ) - (float)airCellCount*0.00);;// /(float)sGrid->dx;
			}
		}
	}

	//Solve the system using Robert Bridson's incomplete Cholesky PCG solver
	double tolerance;
	int iterations;
	bool success = solver.solve(matrix1, rhs, pressure, tolerance, iterations);
	if(!success) {
		cout<<"WARNING: Pressure solve failed!************************************************\n";
	}
	else{
		for(int j = 1; j < sGrid->nY-1; ++j) {
			for(int i = 1; i < sGrid->nX-1; ++i) {
				int index = i + sGrid->nX*j;
				sGrid->p(j,i) = pressure[index];
			}
		}
#define FAC 1
		{
			double scale = dt / (1 * sGrid->dx); // book define rho value before use
			for(int y=1; y < sGrid->nY-1;y++)
				for(int x=1; x < sGrid->nX-1;x++){
					if(sGrid->cellType(y,x) == FLUID && sGrid->cellType(y,x-1) == AIR){
						sGrid->u(y,x) -= scale * ((sGrid->p(y,x)*FAC) - (sGrid->p(y,x-1)*FAC))/FAC;
					}

					if(sGrid->cellType(y,x) == FLUID && sGrid->cellType(y-1,x) == AIR){
						sGrid->v(y,x) -= scale * ((sGrid->p(y,x)*FAC) - (sGrid->p(y-1,x)*FAC))/FAC;
					}

					if(sGrid->cellType(y,x) == FLUID){
						sGrid->u(y,x+1) -= scale * ((sGrid->p(y,x+1)*FAC) - (sGrid->p(y,x)*FAC))/FAC;
						sGrid->v(y+1,x) -= scale * ((sGrid->p(y+1,x)*FAC) - (sGrid->p(y,x)*FAC))/FAC;
					}
				}
/*
			cout<<"Pressure Divergence : y "<<endl;
			for(int y=1; y < sGrid->nY-1;y++,cout<<endl)
				for(int x=1; x < sGrid->nX-1;x++){
					cout<<scale*((sGrid->p(y+1,x)*FAC) - (sGrid->p(y,x)*FAC))/FAC<<" ";
				}
			cout<<"Pressure Divergence : x "<<endl;
			for(int y=1; y < sGrid->nY-1;y++,cout<<endl)
				for(int x=1; x < sGrid->nX-1;x++){
					cout<<scale*((sGrid->p(y,x+1)*FAC) - (sGrid->p(y,x)*FAC))/FAC<<" ";
				}
*/
			for(int y=1; y < sGrid->nY-1;y++)
							for(int x=1; x <= sGrid->nX-1;x++)
								if(fabs(sGrid->u(y,x))<0.001)
									sGrid->u(y,x) = 0;

						for(int y=1; y <=sGrid->nY-1;y++)
							for(int x=1; x < sGrid->nX-1;x++)
								if(fabs(sGrid->v(y,x))<0.001)
									sGrid->v(y,x) = 0;

			setValidVelocity(1);
			for(int y=1; y < sGrid->nY-1;y++)
				for(int x=1; x <= sGrid->nX-1;x++)
					if(!uValid(y,x))
						sGrid->u(y,x) = 0;

			for(int y=1; y <=sGrid->nY-1;y++)
				for(int x=1; x < sGrid->nX-1;x++)
					if(!vValid(y,x))
						sGrid->v(y,x) = 0;
			setValidVelocity(0);
		}
	}
}

void FluidSim :: addViscosity(double visc,double dt)
{
	setValidVelocity(1);
	sGrid->u = addVisc_Helper(sGrid->u,visc,dt,1);
	sGrid->v = addVisc_Helper(sGrid->v,visc,dt,2);
	setValidVelocity(0);
}

#define forEach(i,start,end) for(int i=start;i<=end;i++)

matrix<double> FluidSim :: addVisc_Helper(matrix<double> mat, double visc, double dt, int comp){
	matrix<double> temp = mat ;
	temp.clear();
	double a = dt/(sGrid->dx * sGrid->dx) * visc ;
	//cout<<"a = "<<visc<<" "<<a<<endl ;
	forEach(k,0,500){
		forEach(i,1,(int)mat.size1()-2){
			forEach(j,1,(int)mat.size2()-2){
				double nval = 0 ; //neighbour val
				double ncnt = 0;//neighbour cnt
				if( (comp==1 && uValid(i,j)) || (comp==2 && vValid(i,j)))
				{
					ncnt = getNeighbours(temp,comp,i,j,nval);
					temp(i,j) = (mat(i,j) + a*nval)/(1+ncnt*a) ;
				}
			}
		}
	}
	return temp;
}

#define IsInBound(a,b,m,n) ( (a)>=1&&(b)>=1&&(a)<=(m-1)&&(b)<=(n-1)?true:false )

int FluidSim :: getNeighbours(matrix<double>mat, int comp, int i, int j, double &nval){
	int ncnt = 0;
	int m = mat.size1()-1;
	int n = mat.size2()-1;

	switch(comp){
		case 1 :
				if(IsInBound(i,j-1,m,n) && uValid(i,j-1)){
					nval += mat(i,j-1);
					ncnt++;
				}
				if(IsInBound(i+1,j,m,n) && uValid(i+1,j)){
					nval += mat(i+1,j);
					ncnt++;
				}
				if(IsInBound(i,j+1,m,n) && uValid(i,j+1)){
					nval += mat(i,j+1);
					ncnt++;
				}
				if(IsInBound(i-1,j,m,n) && uValid(i-1,j)){
					nval += mat(i-1,j);
					ncnt++;
				}
				break;
		case 2 :
				if(IsInBound(i,j-1,m,n) && vValid(i,j-1)){
					nval += mat(i,j-1);
					ncnt++;
				}
				if(IsInBound(i+1,j,m,n) && vValid(i+1,j)){
					nval += mat(i+1,j);
					ncnt++;
				}
				if(IsInBound(i,j+1,m,n) && vValid(i,j+1)){
					nval += mat(i,j+1);
					ncnt++;
				}
				if(IsInBound(i-1,j,m,n) && vValid(i-1,j)){
					nval += mat(i-1,j);
					ncnt++;
				}
				break;
		default:
				break;
	}
	return ncnt;
}
