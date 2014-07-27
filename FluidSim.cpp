/*
 * initFluidBody can be changed...
 *
 */

#include "FluidSim.h"
#include "Printer.h"
//#include "Vec.h"
#include <sys/time.h>
#include <math.h>
#include "pcgsolver/util.h"
#include <cassert>
using namespace std;

void FluidSim :: init(GridStag* sGrid)
{
	this->sGrid = sGrid;
	uValid.resize(sGrid->nY,sGrid->nX+1);
	vValid.resize(sGrid->nY+1,sGrid->nX);
	uValid.clear();
	vValid.clear();
}

void  FluidSim :: simulate(double timestep)
{
	struct timeval tt1, tt2;
	
	gettimeofday(&tt1, NULL);
	float cflVal = cfl();
	gettimeofday(&tt2, NULL);
	
	double cflTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;

	dt = cflVal > 0 && cflVal< 1 ? cflVal:0.1f ;
	{   //mark fluid cells
		gettimeofday(&tt1, NULL);		
		advectParticles(sGrid->fluidParticles, sGrid->u,sGrid->v, dt);
		gettimeofday(&tt2, NULL);		
		double advectParticlesTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;

		gettimeofday(&tt1, NULL);
		markFluidCells();
		gettimeofday(&tt2, NULL);
		double markFluidCellsTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;


		matrix<double > u = sGrid->u;
		matrix<double > v = sGrid->v;
		//advect velocity
		gettimeofday(&tt1, NULL);
		calculateLevelSetDistance();
		gettimeofday(&tt2, NULL);
		double calculateLevelSetDistanceTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;

		gettimeofday(&tt1, NULL);
		sGrid->u = advect2DSelf(sGrid->u,dt,u,v,1);
		sGrid->v = advect2DSelf(sGrid->v,dt,u,v,2);
		gettimeofday(&tt2, NULL);
		double advect2DSelfTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;

		gettimeofday(&tt1, NULL);
		applyBoundaryConditions(VELOCITY_BC2);
		gettimeofday(&tt2, NULL);
		double applyBoundaryConditionsTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;

		//add Gravity
		gettimeofday(&tt1, NULL);
		sGrid->v = addGravity(sGrid->v,dt);
		gettimeofday(&tt2, NULL);
		double addGravityTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;
		applyBoundaryConditions(VELOCITY_BC2);

		//add Viscosity
		//		addViscosity(0,dt);
		//		applyBoundaryConditions(VELOCITY_BC2);

		//apply Pressure
		gettimeofday(&tt1, NULL);
		solvePressureBridson((float)dt);
		gettimeofday(&tt2, NULL);
		double solvePressureBridsonTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;
		applyBoundaryConditions(VELOCITY_BC2);

		//extrapolation
		gettimeofday(&tt1, NULL);
		setValidVelocity(1);
		extrapolate2D(sGrid->u,uValid);
		extrapolate2D(sGrid->v,vValid);
		setValidVelocity(0);
		gettimeofday(&tt2, NULL);
		double extrapolate2DTime = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;
		applyBoundaryConditions(VELOCITY_BC2);

		//advect Particles
		//advectParticles(sGrid->fluidParticles, sGrid->u,sGrid->v, dt);
		//markFluidCells();
		cout<<"cflTime: "<<cflTime<<" advectParticlesTime: "<<advectParticlesTime<<" advect2DSelfTime: "
				<<advect2DSelfTime<<" applyBoundaryConditionsTime: "<<applyBoundaryConditionsTime
				<<" addGravityTime: "<<addGravityTime<<" solvePressureBridsonTime: "
				<<solvePressureBridsonTime<<" extrapolate2DTime: "<<extrapolate2DTime
				<<" calculateLevelSetDistanceTime : "<<calculateLevelSetDistanceTime<<endl;;
	}

}

double FluidSim ::  cfl() //keep
{
	double maxVel = 0.0;
	for(int y = 0;y < sGrid->nX; y++)
		for(int x = 0;x < sGrid->nY+1; x++){
			if(maxVel < fabs(sGrid->u(y,x)))
				maxVel = fabs(sGrid->u(y,x));
		}
	for(int y = 0;y < sGrid->nX+1; y++)
		for(int x = 0;x < sGrid->nY; x++){
			if(maxVel < fabs(sGrid->v(y,x)))
				maxVel = fabs(sGrid->v(y,x));
		}
	if ( !maxVel )
		return 0;
	return (sGrid->dx/maxVel);
}

void FluidSim :: advectParticles(std::vector <Particles*> & plist, matrix<double> u, matrix<double>v, double dt)//keep
{
	double dx = sGrid->dx;
	int part_rows = plist.size()/ NTHREADS;
//	omp_set_num_threads(NTHREADS);
	double posx,posy;
//	#pragma omp parallel for default(none) \
	                     shared(plist,u,v,dx,part_rows)\
	                     private(posx,posy)
	{	
//		#pragma omp for schedule(guided,part_rows)
	
	for ( unsigned i = 0; i < plist.size(); i++){
	
		posx = plist.at(i)->x /dx ;
		posy = plist.at(i)->y /dx ;
		RK2(posx, posy, u, v, dt);
		plist.at(i)->x = posx*dx;
		plist.at(i)->y = posy*dx;
	}
	}
//	omp_set_num_threads(1);
}
matrix<double> FluidSim :: advect2DSelf(matrix<double> q, double dt, matrix<double> u, matrix<double> v,int component)//keep
{
	//proper advection - RK2
	matrix<double> temp=q;
	temp.clear();
	int nX = this->sGrid->nX;
	int nY = this->sGrid->nY;
	double dx = sGrid->dx;
	//dt*=nX;
	double x,y, posx, posy;
	int part_rows = nX/ NTHREADS;
	int th_id;
//	omp_set_num_threads(NTHREADS);
	if(component==1){ //Horizontal Component
//		#pragma omp parallel shared(nX,nY,dx,temp,part_rows) private(x,y,posx,posy, th_id)
		{
			// th_id = omp_get_thread_num(); //th_id holds the thread number for each thread

//		#pragma omp for schedule(guided,part_rows)
		for(int i=1;i<=nY-2;i++){
			for(int j=1;j<=nX-1;j++){
		//		printf("Thread #%d is doing row %d.\n",th_id,i); //Uncomment this line to see which thread is doing each row
      
				x = (j)*dx;
				y = (i+0.5)*dx;
				posx = x/dx;
				posy = y/dx;
				RK2(posx,posy,sGrid->u,sGrid->v,-dt);
				temp(i,j) = 1 ;// getVelInterpolated(posx,posy-0.5,sGrid->u);
			}
		}
		}
	}
	else{ //component=2 i.e. Vertical component
//		#pragma omp parallel shared(nX,nY,dx,temp,part_rows) private(x,y,posx,posy)
		{
//		#pragma omp for schedule(guided,part_rows)
		for(int i=1;i<=nY-1;i++){
			for(int j=1;j<=nX-2;j++){
				x = (j+0.5)*dx;
				y = (i)*dx;
				posx = x/dx;
				posy = y/dx;
				RK2(posx,posy,sGrid->u,sGrid->v,-dt);
				temp(i,j) = getVelInterpolated(posx-0.5,posy,sGrid->v);
			}
		}
		}
	}
//	omp_set_num_threads(1);
	return temp;
}

matrix<double> FluidSim :: addGravity(matrix<double> ua, double dt) //keep
{
	/* here we are not adding gravity factor to top line of grid i.e.(nY,0-nX) */
	matrix<double> ub = ua;
	setValidVelocity(1);
	for(int y=1; y < sGrid->nY-1; y++)
		for(int x=1; x < sGrid->nX-1; x++){
			if(vValid(y,x)){// && (fabs(ub(y,x)) > dt || fabs(ub(y,x))==0.0))
				ub(y,x) += dt*GRAVITY;
			}
		}
	setValidVelocity(0);
	return ub;
}

//-----PRESSURE SOLVER----
void FluidSim::solvePressureBridson(float dt) { //keep

	unsigned int sys_size = sGrid->nX*sGrid->nY;
	if(rhs.size() != sys_size) {
		rhs.resize(sys_size);
		pressure.resize(sys_size);
		matrix1.resize(sys_size);
	}
	matrix1.zero();
	//Build the linear system for pressure
	double scale =  (dt / (double)(sGrid->dx * sGrid->dx));
	matrix<double> cellType = sGrid->cellType;
	for(int j = 1; j < sGrid->nY-1; ++j) {
		for(int i = 1; i <  sGrid->nX-1; ++i) {
			int index = i + sGrid->nX*j;
			rhs[index] = 0;
			pressure[index] = 0;

			if(cellType(j,i)==FLUID)
			{
				if(cellType(j,i+1)==FLUID) {
					matrix1.add_to_element(index, index, scale);
					matrix1.add_to_element(index, index + 1, -scale);
				}
				else if(cellType(j,i+1)==AIR) {
					matrix1.add_to_element(index, index, scale);
				}
				//left neighbour
				if(cellType(j,i-1)==FLUID) {
					matrix1.add_to_element(index, index, scale);
					matrix1.add_to_element(index, index - 1, -scale);
				}
				else if(cellType(j,i-1)==AIR) {
					matrix1.add_to_element(index, index, scale);
				}

				//top neighbour
				if(cellType(j+1,i)==FLUID) {
					matrix1.add_to_element(index, index, scale);
					matrix1.add_to_element(index, index + sGrid->nX, -scale);
				}
				else  if(cellType(j+1,i)==AIR)
				{
					matrix1.add_to_element(index, index, scale);
				}

				//bottom neighbour
				if(cellType(j-1,i)==FLUID) {
					matrix1.add_to_element(index, index, scale);
					matrix1.add_to_element(index, index - sGrid->nX, -scale);
				}
				else if(cellType(j-1,i)==AIR) {
					matrix1.add_to_element(index, index, scale);
				}

				rhs[index]  = -((
						(sGrid->u(j,i+1) - sGrid->u(j,i) + sGrid->v(j+1,i) - sGrid->v(j,i))/(double)sGrid->dx ) );;// /(float)sGrid->dx;
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
		double scale = dt / (1 * sGrid->dx); // book define rho value before use
		for(int y=1; y < sGrid->nY-1;y++)
			for(int x=1; x < sGrid->nX-1;x++){
				if(sGrid->cellType(y,x) == FLUID && sGrid->cellType(y,x-1) == AIR){
					sGrid->u(y,x) -= scale * ((sGrid->p(y,x)) - (sGrid->p(y,x-1)));
				}

				if(sGrid->cellType(y,x) == FLUID && sGrid->cellType(y-1,x) == AIR){
					sGrid->v(y,x) -= scale * ((sGrid->p(y,x)) - (sGrid->p(y-1,x)));
				}

				if(sGrid->cellType(y,x) == FLUID){
					sGrid->u(y,x+1) -= scale * ((sGrid->p(y,x+1)) - (sGrid->p(y,x)));
					sGrid->v(y+1,x) -= scale * ((sGrid->p(y+1,x)) - (sGrid->p(y,x)));
				}
			}

		/*for(int y=1; y < sGrid->nY-1;y++)
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
		 */		}

}

//-----Extrapolation........
void FluidSim :: extrapolate2D(matrix<double> &grid, matrix<int> &valid) //keep
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

void FluidSim :: calculateLevelSetDistance(){//keep
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

double FluidSim :: getMinDistance(double a, double b, double r){
	double d ;

	if (fabs(a-b) < 1 )
		d = (a + b + sqrt(2-(a-b)))/2.0;
	else
		d = (a>b?b:a) + 1;

	return(d<r?d:r);
}


void  FluidSim :: setValidVelocity(int val) //keep
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

void FluidSim :: applyBoundaryConditions(int bc)//boundary Condition //keep-BC2
{
	//bottom boundary
	cout<<"ok"<<endl;
	for(int x = 0; x < sGrid->nX; x++){
		sGrid->v(0,x) =  0;//-sGrid->v(2,x);
		sGrid->v(1,x) = 0;//sGrid->v(2,x);
	}
	//Top boundary..
	for(int x = 0; x < sGrid->nX; x++){  //allowing in x direction
		sGrid->v(sGrid->nY-1,x) = 0;//sGrid->v(sGrid->nY-2,x);
		sGrid->v(sGrid->nY,x) = 0;//-sGrid->v(sGrid->nY-2,x);
	}

	for(int y = 0; y < sGrid->nY; y++){
		sGrid->u(y,0) = 0;//-sGrid->u(y,2);
		sGrid->u(y,1) = 0;//sGrid->u(y,2);
	}
	//Right boundary..
	for(int y = 0; y < sGrid->nY; y++){
		sGrid->u(y,sGrid->nX-1) = 0;//-sGrid->u(y,sGrid->nX-2);
		sGrid->u(y,sGrid->nX) = 0;//sGrid->u(y,sGrid->nX-2);
	}

}

void FluidSim :: RK2(double &posx, double &posy,matrix<double> &u, matrix<double> &v, double dt){
	double dx = sGrid->dx;
	double x = posx*dx;
	double y = posy*dx;
	int nX = sGrid->nX;
	int nY = sGrid->nY;
	double velx = getVelInterpolated(posx,posy-0.5,u);
	double vely = getVelInterpolated(posx-0.5,posy,v);
	posx = x  + 0.5*dt*velx;
	posy = y  + 0.5*dt*vely;

	posx = posx < (1)*dx ? (1)*dx:posx;
	posy = posy < (1)*dx ? (1)*dx:posy;
	posx = posx > (nX-1)*dx ? (nX-1)*dx:posx;
	posy = posy > (nY-1)*dx ? (nY-1)*dx:posy;

	posx/=dx;
	posy/=dx;
	velx = getVelInterpolated(posx,posy-0.5,u);
	vely = getVelInterpolated(posx-0.5,posy,v);
	posx = x  + dt*velx;
	posy = y  + dt*vely;

	posx = posx < (1)*dx ? (1)*dx:posx;
	posy = posy < (1)*dx ? (1)*dx:posy;
	posx = posx > (nX-1)*dx ? (nX-1)*dx:posx;
	posy = posy > (nY-1)*dx ? (nY-1)*dx:posy;

	posx/=dx;
	posy/=dx;
}



matrix<double> FluidSim :: addForce(matrix<double> dest, double dt, matrix<double> src)
{
	for(int i=1;i < sGrid->nY-1;i++) //exclude the boundary cells
		for(int j=1;j<sGrid->nX-1;j++)
			dest( i, j ) = dest( i, j ) + dt* src(i,j);
	return dest;
}

double FluidSim :: getVelInterpolated(double x,double y, matrix<double> &mat)
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
		fluidPosBX = (nX-1)-(nX-1)/3;  //last one is boundary
		fluidPosBY = 1 ;//sb+1;
		fluidPosTX = nX-1;
		fluidPosTY = ((nY-1) - (nY-1)/3) ;
		for (int i = 0; i < (nX)*(nY)*PARTICLE_PER_CELL; ++i) {
			float xpos = randhashf(i * 2, 0.0, zoomFactor);
			float ypos = randhashf(i * 2 + 1,0.0, zoomFactor);

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
		fluidPosBY = 0+3*(nX)/10;  //last one is boundary
		fluidPosBX = 0+3*(nY)/10;
		fluidPosTY = fluidPosBY + 6*(nX)/10;
		fluidPosTX = fluidPosBX + 4*(nY)/10;

		for (int i = 0; i < (nX)*(nY)*PARTICLE_PER_CELL; ++i) {
			float xpos = randhashf(i * 2, 0.0, zoomFactor);
			float ypos = randhashf(i * 2 + 1, 0.0, zoomFactor);


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
		fluidPosBY = 0+1*(nX)/20;  //last one is boundary
		fluidPosBX = 0+1*(nY)/20;
		fluidPosTY = fluidPosBY + 6*(nX)/20;
		fluidPosTX = 32 - fluidPosBX ;

		for (int i = 0; i < (nX)*(nY)*PARTICLE_PER_CELL; ++i) {
			float xpos = randhashf(i * 2, 0.0, zoomFactor);
			float ypos = randhashf(i * 2 + 1, 0.0, zoomFactor);


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
		fluidPosBY = 0+1*(nX)/20;  //last one is boundary
		fluidPosBX = 0+1*(nY)/20;
		fluidPosTY = fluidPosBY + 12*(nX)/20;
		fluidPosTX = fluidPosBX + 6*(nX)/20 ;

		for (int i = 0; i < (nX)*(nY)*PARTICLE_PER_CELL; ++i) {
			float xpos = randhashf(i * 2, 0.0, zoomFactor);
			float ypos = randhashf(i * 2 + 1, 0.0, zoomFactor);

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

		fluidPosBX = 0+19*(nX)/20 - 5*(nX)/20 ;  //last one is boundary
		fluidPosBY = 0+1*(nX)/20;
		fluidPosTY = fluidPosBY + 12*(nY)/20;
		fluidPosTX = 19*(nY)/20 + 1;

		for (int i = 0; i < (nX)*(nY)*PARTICLE_PER_CELL; ++i) {
			float xpos = randhashf(i * 2, 0.0, zoomFactor);
			float ypos = randhashf(i * 2 + 1, 0.0, zoomFactor);

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

void FluidSim :: initSolidBoundary(int choice)//remove this choice factor..
{
	//bottom boundary..
	for(int x = 0; x < sGrid->nX; x++){
		sGrid->cellType(0,x) = SOLID;
	}
	//Top boundary..
	for(int x = 0; x < sGrid->nX; x++){
		sGrid->cellType(sGrid->nY-1,x) = SOLID;
	}
	//Left boundary..
	for(int y = 0; y < sGrid->nY; y++){
		sGrid->cellType(y,0) = SOLID;
	}
	//Right boundary..
	for(int y = 0; y < sGrid->nX; y++){
		sGrid->cellType(y,sGrid->nX-1) = SOLID;
	}
	this->applyBoundaryConditions(VELOCITY_BC2);
}

void FluidSim :: markFluidCells() //keep
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



//---VISCOSITY....Not That Perfect....
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
