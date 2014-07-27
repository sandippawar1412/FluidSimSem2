#include "commonData.h"
#include "main.h"
#include <time.h>
#include <pthread.h>
#include <omp.h>
#include <sys/time.h>
/*
int main(int argc, char** argv)
{
	struct timeval tt1, tt2;
	init();
	int i=pthread_getconcurrency();
	int it=0;
	int ni= 0;
	while(ni++<NO_OF_ITERATIONS){
		gettimeofday(&tt1, NULL);
		animate();
		gettimeofday(&tt2, NULL);
		int milliSeconds = (tt2.tv_sec - tt1.tv_sec) * 1000 + (tt2.tv_usec - tt1.tv_usec)/1000;
		cout<<"Iteration "<<it<<" : "<<milliSeconds<<"ms"<<endl<<endl;
		it++;
	}
	return 0;
}
*/


void animate()
{
	fluidSim->simulate(timestep);
} 

void init(void) 
{
   sGrid->initGridStag();
   fluidSim->init(sGrid);
   fluidSim->initFluidBody(fBT);// 2: indicates dam break center
   print->init(sGrid);
   render->init(sGrid);
   
}

void initParticles()
{
	sGrid->fluidParticles.clear();
}

/*------------------Code Below this is used when rendered is ON---------*/
void display ( void ) ;
void preDisplay ( );
void postDisplay ( );
void openGlutWindow ( char* windowName) ;
void reshape ( int w, int h ) ;
void idleFun();


int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   init ();
   int i=pthread_getconcurrency();

   char windowName[]="   Liquid_Simulator-LevelSet+Surface" ;
   openGlutWindow(windowName);

   glutMainLoop();
   return 0;
}

void display(void){
	preDisplay();
	static bool flag[10]={true,true,true,true,true,false};

	char output1 = ' ';
	bool anyUpdation = false;

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
			         " $"<<flag[4]<<" %"<<flag[5]/*<<" ^"<<flag[6]*/<<"   +"<<endl;
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
void idleFun ( void )
{
	clock_t t1=clock(),t2;
	animate();
	t2=clock();
	double diff = t2-t1;
	cout<<"display "<<" : "<<diff/CLOCKS_PER_SEC*1000<<"ms"<<endl;
	glutSetWindow ( winId );
	glutPostRedisplay ( );
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
   glutDisplayFunc(display);
//   glutSpecialFunc(&SpecialKeyPressed);
//   glutKeyboardFunc(&KeyPressed);
   glutReshapeFunc(reshape);
   glutIdleFunc(idleFun);
}
