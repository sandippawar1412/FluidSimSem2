#include "GridStag.h"
#include "Renderer.h"
#include "FluidSim.h"
#include "keyboard.h"
#include "Printer.h"
#include "Plotter.h"
#include "opengl.h"
#include <string>

#include <pthread.h>

GridStag* sGrid = new GridStag;
Renderer* render = new Renderer;
Printer * print = new Printer;
FluidSim* fluidSim = new FluidSim;
Plotter* plotr = new Plotter ;

int winId;
unsigned int winSizeX = 600;
unsigned int winSizeY = 600;
int  omx, omy, mouseButton[3];
double mx,my;

double timestep = 0.01;

void display ( void ) ;
void init ( void ) ;
void animate () ;


void mouseFun ( int button, int state, int x, int y ) ;
void motionFun ( int x, int y ) ;	
void preDisplay ( );
void postDisplay ( );
void openGlutWindow ( char* windowName) ;
void reshape ( int w, int h ) ;


unsigned int framenum=0;
unsigned char *pRGB;

int SCREEN_WIDTH=600;
int SCREEN_HEIGHT=600;

void capture_frame(unsigned int framenum)
{
  //global pointer float *pRGB
  pRGB = new unsigned char [3 * (SCREEN_WIDTH+1) * (SCREEN_HEIGHT + 1) ];
  // set the framebuffer to read
  //default for double buffered
  glReadBuffer(GL_BACK);

  glPixelStoref(GL_PACK_ALIGNMENT,1);//for word allignment

  glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, pRGB);

  char filename[200];
  sprintf(filename,"frame_%04d.ppm",framenum);
  std::ofstream out(filename, std::ios::out);
  out<<"P6"<<std::endl;
  out<<SCREEN_WIDTH<<" "<<SCREEN_HEIGHT<<" 255"<<std::endl;
  out.write(reinterpret_cast<char const *>(pRGB), (3 * (SCREEN_WIDTH+1) * (SCREEN_HEIGHT + 1)) * sizeof(int));
  out.close();

  //function to store pRGB in a file named count
  delete pRGB;
}


void Timer(int val) {
	if(framenum>100)
		fluidSim->simulate(timestep); //timestep is actually not used
    glutPostRedisplay();
    glutTimerFunc(FRAME_INTERVAL_ms, Timer, 1);
}

int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   init ();
   int i=pthread_getconcurrency();/*to resolve this error:Inconsistency detected by ld.so: dl-version.c: 224: _dl_check_map_versions: Assertion `needed != ((void *)0)' failed!
   ..Read this : https://bugs.launchpad.net/ubuntu/+source/nvidia-graphics-drivers-319/+bug/1248642?comments=all*/

   char windowName[]="   Liquid_Simulator-LevelSet+Surface" ;
   openGlutWindow(windowName);
 //  glutTimerFunc(FRAME_INTERVAL_ms,Timer,1);  //1s/40ms = 1000/40 -> 25 frames/sec
   glutMainLoop();
   return 0;
}

//---------glut and opengl.....Initializations
void animate()
{
	fluidSim->simulate(timestep);
#define MAX_ITERATION 1800
	/*

	extern int flagC;
	static int it = 0;

	if(flagC==2){
		it=0;
		flagC=10;
	}
	cout<< ++it <<endl ;
*/
} 

void init(void) 
{
   extern int flagC ;
   flagC = 0;
   sGrid->initGridStag();
   render->init(sGrid);
   print->init(sGrid);
   fluidSim->init(sGrid);
   fluidSim->initFluidBody(fBT);// 1: indicates dam break rhb
   plotr->init(sGrid,200,26,MAX_ITERATION);//number of intervals,noofFiles,maxIteration
}

void initParticles()
{
	sGrid->fluidParticles.clear();
}

void display(void){
	preDisplay();
	static bool flag[10]={false};
	switch(output1){

		case'~':
			flag[0] = !flag[0];
			//render->renderBoundary
			break;
		case '!':
			flag[1] = !flag[1];
			//render->renderGrid();
			break;
		case '@':
			flag[2] = !flag[2];
			//render->renderParticles();
			break;
		case '#':
			flag[3] = !flag[3];
			//render->renderSurfaceBoundary();
			break;
		case '$':
			flag[4] = !flag[4];
			//render->renderVector2D(sGrid->u,sGrid->v);
			break;
		case '%':
			flag[5] = !flag[5];
			//render->renderMat(sGrid->distanceLevelSet,1);
			break;
		/*case '^':
			flag[6] = !flag[6];
			//render->renderMat(sGrid->isFluidBoundary,1);
			break;*/
	}
	if(anyUpdation){
		cout<<"Flags :"<<" ~"<<flag[0]<<" !"<<flag[1]<<" @"<<flag[2]<<" #"<<flag[3]<<
			         " $"<<flag[4]<<" %"<<flag[5]/*<<" ^"<<flag[6]*/<<"   +"<<idleOn<<endl;
		anyUpdation = false;
	}
	if(flag[0])
		render->renderBoundary();
	if(flag[1])
		render->renderGrid();
	if(flag[2])
		render->renderParticles();
	if(flag[3])
		render->renderSurfaceBoundary();
	if(flag[4])
		render->renderVector2D(sGrid->u,sGrid->v);
	if(flag[5])
		render->renderMat(sGrid->distanceLevelSet,2);
	/*if(flag[6])
		render->renderMat(sGrid->isFluidBoundary,1);
*/
	output1 = ' ';
	postDisplay();


}



void display1(void)
{
   preDisplay();
   //print->matrices(sGrid->distanceLevelSet);
   switch(output)
   {
		case VEL :
			render->renderVector2D(sGrid->u,sGrid->v);
			break;
		case DEN :
			render->renderDen2D_Stam();	
			break;
		case PAR :
			render->renderParticles();
			break;
		case DUB :
			render->renderVector2D(sGrid->u,sGrid->v);
			render->renderParticles();
			break;
		case PRE :
			//render->renderMat(sGrid->p,2);
			//render->renderMat(sGrid->cellType,1);
			render->renderParticles();
			render->renderSurfaceBoundary();

//			render->renderMat(sGrid->distanceLevelSet,1);
			break;
		case LES :
				//render->renderMat(sGrid->p,2);
			render->renderMat(sGrid->distanceLevelSet,1);
			break;

		case BOU :
					//render->renderMat(sGrid->p,2);
		render->renderMat(sGrid->isFluidBoundary,1);

					break;
   }	   
	if(AdvectFlagCell || AdvectFlagParticle){
		render->renderParticle(traceAdvectCellsd[2],traceAdvectCellsd[3], 0.6,0,0,12);//R G B Size
		render->renderParticle(traceAdvectCellsd[0],traceAdvectCellsd[1], 0,1,0,8);//R G B Size
	}
//	capture_frame(framenum++);
//	cout<<framenum++<<endl;
	postDisplay();
}

void idleFun ( void )
{
	if(idleOn){
	 animate();
	 glutSetWindow ( winId );
	 glutPostRedisplay ( );}
}

void motionFun ( int x, int y )
{
	mx = x;
	my = y;
}

void processMouseInCoords(double mx, double my)
{
	cout<<endl<<"Clicked At     : "<<mx<<" "<<my<<endl;
	if(AdvectFlagCell){
		fluidSim->advect2DCell((int)mx,(int)my,fluidSim->dt,sGrid->u,sGrid->v,2);
	}
	else if(AdvectFlagParticle){
		fluidSim->advect2DParticle(mx, my,fluidSim->dt,sGrid->u,sGrid->v);
	}
}
void mouseFun ( int button, int state, int x, int y )
{
	mouseButton[button] = ( state == GLUT_DOWN ) ;	
	if(mouseButton[0]){
	mx = ( (x/(double)winSizeX) * (sGrid->nX) ) ;
	my = ( ((winSizeY - y) / (double)winSizeY) * ( sGrid->nY) ); //not covers exact cell..
	
	processMouseInCoords(mx,my);
	}
}

void reshape (int w, int h)
{
   glutSetWindow(winId);
   glutReshapeWindow(w ,h);
   winSizeX = w;
   winSizeY = h;
}

void preDisplay()
{
	glViewport ( 0, 0, winSizeX, winSizeY);
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho( 0, zoomFactor, 0, zoomFactor	, 0,1 ); //better to use ortho..
	glClearColor( 0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}

void postDisplay()
{
   glutSwapBuffers();
}	

void openGlutWindow(char* windowName)
{
   glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
   glutInitWindowSize (winSizeX, winSizeY); 
   glutInitWindowPosition (570, 100);
 
   winId = glutCreateWindow (&windowName[2]);

   glClearColor( 0.0f, 0.0f, 0.0f, 1.0f);
   glClear(GL_COLOR_BUFFER_BIT);
   glutSwapBuffers();
   glClear(GL_COLOR_BUFFER_BIT);
   glutSwapBuffers();
   createMenu();
   glutDisplayFunc(display); 
   glutSpecialFunc(&SpecialKeyPressed);
   glutKeyboardFunc(&KeyPressed);
   glutReshapeFunc(reshape);
   glutMouseFunc(mouseFun);
   glutMotionFunc(motionFun);
   glutIdleFunc(idleFun);
}


void menu(int num){
  if(num == 27){
    glutDestroyWindow(winId);
    exit(0);
  }
  switch(num)
  {
	case 0:	
	case 1:	
	case 2:
			ProjectFLAG = num  ;
			break;
	case 3:	
	case 4:	
	case 5:
	case 6:
			output = (renderFlag)num;
			break;
  }
  glutPostRedisplay();
}
