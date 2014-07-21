#ifndef _KEYBOARD_H
#define _KEYBOARD_H

#include "opengl.h"
#include "commonData.h"
/*
#define VEL 1	
#define DEN 2	
#define PAR 3	
#define PRE 4
*/
enum renderFlag{VEL=4, DEN, PAR, PRE, DUB, BOU, LES};//choose what to render..
renderFlag output=LES; //default is velocity
char output1='~';//extra.. to renderFlag.. it takes keys SHIFT + < ~ to ^ >  and ons offs alternatively different render
//enum ProjectionRoutines{ GSR,PCG,PCGBridson} ProjectFLAG ;
bool idleOn = true;// to run the simulator on idle
bool anyUpdation = true;

double nearFactor = 0;
double farFactor = zoomFactor; 

void animate();
void init();
/*--checkIt--*/

//bool GSrFlag = false; //flag to decide the project routine to be used..
bool AdvectFlag = true;
bool AdvectFlagP = true;
bool extrapolateFlag = true;

bool AdvectFlagCell = false;
bool AdvectFlagParticle = false;

#define GSRflag 0
#define PCGflag 1
#define PCGBridsonflag 2



int ProjectFLAG  = PCGBridsonflag;
	
void initParticles();	
extern int winId;


enum fluidBody{ DAM_BREAK=1,DAM_CENTER=2,STATIC_BED=3,DOUBLE_DAM=4} fBT=DAM_CENTER ; //fluidBodyType

void KeyPressed (unsigned char key, int x, int y)
{
	if (key==27)
		exit(0);
	switch( (key) )
	{
		default:
			output1 = key;
			anyUpdation = true;
			break;
		case '+':
			idleOn = !idleOn;
			anyUpdation = true;
			break;
		case 'n' :
				animate();
				//fluidSim->simulate();
				glutSetWindow ( winId );
				traceCellFlag = OFF;
				glutPostRedisplay();
				break;
				
		case 'r' :
				output = PRE;
				glutPostRedisplay();
				break;
				
		case 't' :
				output = VEL;
				glutPostRedisplay();
				break;
		/*case 'y':
				output = DEN	;
				glutPostRedisplay();	
				break;
		*/
		case 'u':
				output = DUB	;
				glutPostRedisplay();	
				break;
		case 'o':
				output = BOU	;
				glutPostRedisplay();
				break;
		case 'y':
				output = LES	;
				glutPostRedisplay();
				break;
		case 'p':
				output = PAR	;
				glutPostRedisplay();	
				break;
		case 'm':
				fBT= DAM_CENTER;
				init();
				flagC=2;
				glutPostRedisplay();
				break;	
		case 'l':
				fBT= DAM_BREAK;
				init();
				flagC=2;
				glutPostRedisplay();
				break;	
		case 's':
				fBT= fBT==STATIC_BED?DOUBLE_DAM:STATIC_BED;
				init();
				flagC=2;
				break;	
		case 'g':
				ProjectFLAG = (ProjectFLAG+1)%4;
				{char str[][15] = {"GSR","PCG","PCGBridson","NoPressure"};
				cout<<"\nProjectFlag : "<<str[ProjectFLAG]<<endl;}
				break;
				
		case 'q':
				AdvectFlagCell = !AdvectFlagCell;
				AdvectFlagParticle = false;
				break;
		case 'w':
				AdvectFlagCell = false;
				AdvectFlagParticle = !AdvectFlagParticle;
				break;
		case 'a':
				AdvectFlag = !AdvectFlag;
				{char str[][15] = {"RK2","Normal"};
				cout<<"\nAdvectFlag : "<<str[(int)AdvectFlag]<<endl;}
				break;	
		case 'e':
				extrapolateFlag = !extrapolateFlag;
				//glutPostRedisplay();
				{char str[][15] = {"Off","On"};
				cout<<"\nExtrapolation : "<<str[(int)extrapolateFlag]<<endl;}


				break;	
		case 'I':
				initParticles();
				break;

				
				
				
	}
}

void SpecialKeyPressed(int key, int x, int y) 
{
 double xpos = 0.0,
	    ypos = 0.0,
	    zpos = 0.0;

 switch (key)
      {   
       case GLUT_KEY_UP:
							ypos+=0.1; 
							break;
       case GLUT_KEY_DOWN:
							ypos-=0.1; 
							break;
       case GLUT_KEY_LEFT:
							xpos-=0.1;  
							break;
       case GLUT_KEY_RIGHT:
							xpos+=0.1;                
							break;
       case GLUT_KEY_PAGE_DOWN:
								zpos-=0.1;
	                    	    break;
       case GLUT_KEY_PAGE_UP:
							zpos+=0.1;
	                	    break;
       default:
				break;
      }	
	nearFactor+=xpos;
    zoomFactor+=ypos;
	farFactor+=zpos;
    glutPostRedisplay();
}
void menu(int);
void createMenu(void){     
	
	int submenu_id = glutCreateMenu(menu);
    glutAddMenuEntry("Velocity", VEL);
    glutAddMenuEntry("Pressure", PRE);
    glutAddMenuEntry("Particles", PAR);
    glutAddMenuEntry("Density", DEN); 
    glutAddMenuEntry("Mix", DEN); 
        
	int submenu2_id = glutCreateMenu(menu);
    glutAddMenuEntry("GSR", GSRflag);
    glutAddMenuEntry("PCG", PCGflag);
    glutAddMenuEntry("PCG_Bridson", PCGBridsonflag);
        
   glutCreateMenu(menu);
		glutAddSubMenu("ProjectRoutine", submenu2_id);
    	glutAddSubMenu("Render", submenu_id);
    glutAddMenuEntry("Quit", 27);   //ASCII of ESC is 27  
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
} 


#endif
