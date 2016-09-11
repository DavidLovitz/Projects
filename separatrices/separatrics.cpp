
#include "separatrics.h"

//complex double
typedef std::complex<double> dcomp;

//initialize graphics window and opengl settings
void initGL(void)
{
    //set up a square window of size WIN_SIZE
    //lower left corner initialized to (0,0)
	glViewport(0, 0, (GLsizei)WIN_SIZE, (GLsizei)WIN_SIZE);
    //set up camera view (could use GL_MODELVIEW?)
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();//resets a 'global' perspective matrix to the identity matrix?
	glTranslatef(-0.5, -0.5, 0.0);//moves the 'camera' position away from origin
    glScalef(1.0,1.0,1.0);//set scale
    //set textures, shading, blending etc.
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClear(GL_COLOR_BUFFER_BIT);
}

//interpolates a 2x2 tensor using bilinear interpolation and four input tensors at a point (x,y)
//(Assumes that the four tensors are arranged in a unit grid)
icMatrix2x2 bilinear_interpolate(icMatrix2x2 tensor[4], double x, double y){
    int i,j;
    icMatrix2x2 newTensor;
    for(i=0; i<2; i++){
        for(j=0; j<2; j++){
            newTensor.entry[i][j] = tensor[0].entry[i][j]*(1.0-x)*(1.0-y)+tensor[1].entry[i][j]*x*(1.0-y) +tensor[2].entry[i][j]*(1.0-x)*y + tensor[3].entry[i][j]*x*y;
        }
    }
    return newTensor;
}

//another way of finding singularities
//Algorithms computed with Mathematica
//Solves the system T11=T12 T12=0
//Gets the (x,y) coordinate of a singularity within the quad
//stores coordinates in array o[0]=x o[1]=y
//IN PROGRESS
void find_singularity(icMatrix2x2 tensor[4], double o[2])
{
    double T011 = tensor[0].entry[1][1];
    double T012 = tensor[0].entry[1][2];
    double T022 = tensor[0].entry[1][1];
    
    double T111 = tensor[1].entry[1][1];
    double T112 = tensor[1].entry[1][2];
    double T122 = tensor[1].entry[1][1];
    
    double T211 = tensor[2].entry[1][1];
    double T212 = tensor[2].entry[1][2];
    double T222 = tensor[2].entry[1][1];
    
    double T311 = tensor[3].entry[1][1];
    double T312 = tensor[3].entry[1][2];
    double T322 = tensor[3].entry[1][1];
    
    double x1 =
    -((T012*T211-2*T112*T211-T011*T212+T022*T212+2*T111*T212-2*T122*T212-T012*T222+2*T112*T222+T112*T311-T111*T312+T122*T312-T112*T322+sqrt(pow((-T011*T212+T022*T212+T012*(T211-T222)+(T111-T122)*(2*T212-T312)+T112*(-2*T211+2*T222+T311-T322)),2)-4*((-T111+T122)*T212+T112*(T211-T222))*((T011-T022-T111+ T122)*(T212-T312)-(T012-T112)*(T211-T222-T311+T322))))
      /(2*((T011-T022-T111+T122)*(T212-T312)-(T012-T112)*(T211-T222-T311+T322))));
    
    double y1 =(-T012*T211+T011*T212-T022*T212+T012*T222+T112*T311-2*T212*T311-T111*T312+T122*T312+2*T211*T312-2*T222*T312-T112*T322+2*T212*T322+sqrt(pow((-T011*T212+T022*T212+T012*(T211-T222)+(T111-T122)*(2*T212-T312)+T112*(-2*T211+2*T222+T311-T322)),2)-4*((-T111+T122)*T212+T112*(T211-T222))*((T011-T022-T111+T122)*(T212-T312)-(T012-T112)*(T211-T222-T311+T322))))
    /(2*((T111-T122-T211+T222)*(T012-T312)-(T112-T212)*(T011-T022-T311+T322)));
    
    double x2 =
    (-T012*T211+2*T112*T211+T011*T212-T022*T212-2*T111*T212+2*T122*T212+T012*T222-2*T112*T222-T112*T311+T111*T312-T122*T312+T112*T322+sqrt(pow((-T011*T212+T022*T212+T012*(T211-T222)+(T111-T122)*(2*T212-T312)+T112*(-2*T211+2*T222+T311-T322)),2)-4*((-T111+T122)*T212+T112*(T211-T222))*((T011-T022-T111+T122)*(T212-T312)-(T012-T112)*(T211-T222-T311+T322))))/(2*((T011-T022-T111+T122)*(T212-T312)-(T012-T112)*(T211-T222-T311+T322)));
    
    double y2 =
    -((T012*T211-T011*T212+T022*T212-T012*T222-T112*T311+2*T212*T311+T111*T312-T122*T312-2*T211*T312+2*T222*T312+T112*T322-2*T212*T322+sqrt(pow((-T011*T212+T022*T212+T012*(T211-T222)+(T111-T122)*(2*T212-T312)+T112*(-2*T211+2*T222+T311-T322)),2)-4*((-T111+T122)*T212+T112*(T211-T222))*((T011-T022-T111+T122)*(T212-T312)-(T012-T112)*(T211-T222-T311+T322))))/(2*((T111-T122-T211+T222)*(T012-T312)-(T112-T212)*(T011-T022-T311+T322))));
    printf("%f, %f\n",x1,y1);
    printf("%f, %f\n",x2,y2);
    
    if(x1<0 || x1>1 || y1<0 || y1>1){
        o[0] = x2;
        o[1] = y2;
        return;
    }
    o[0] = x1;
    o[1] = y1;
}

//get the (x,y) coordinates of a singularity within the quad
//stores coordinates in array o[0]=x o[1]=y
//basically solves t11=t22 & t12 =0
//written by Xiaofei Gao
void get_singularity(icMatrix2x2 tensor[4], double o[2])
{
    double a[4], b[4], c[4];
    
    for (int j = 0; j < 4; j++)
    {
        a[j] = 0.5*(tensor[j].entry[0][0] - tensor[j].entry[1][1]);
        b[j] = tensor[j].entry[0][1];
        c[j] = -0.5*(tensor[j].entry[0][0] - tensor[j].entry[1][1]);
    }
    double t[3];
    t[0] = a[0] * (b[1] - b[2] - b[3])
    + a[1] * (-b[0] + b[3])
    + a[2] * (b[0] - b[3])
    + a[3] * (b[0] - b[1] + b[2]);
    t[1] = a[0] * b[2]
    + a[1] * (-b[3])
    + a[2] * (-b[0] + 2 * b[3])
    + a[3] * (b[1] - 2 * b[2]);
    t[2] = a[3] * b[2] - a[2] * b[3];
    double delta, sol_x[2], sol_y[2];
    if ((t[1] * t[1] - 4 * t[0] * t[2]) >= 0)
        delta = sqrt(t[1] * t[1] - 4 * t[0] * t[2]);
    else
    {
        printf("check delta.\n");
        delta = 0;
    }
    if (t[0]==0)
        printf("t[0] = 0\n");

    sol_y[0] = (-t[1] + delta)* 0.5 / t[0];
    sol_y[1] = (-t[1] - delta)* 0.5 / t[0];
    //printf("original computing\n");

    for (int i = 0; i < 2; i++)
    {
        sol_x[i] = (-sol_y[i] * b[1] + sol_y[i] * b[2] + sol_y[i] * b[3] - b[2])
        / (sol_y[i] * b[0] - sol_y[i] * b[1] + sol_y[i] * b[2] - b[2] + b[3]);
        
    }
    for (int i = 0; i < 2; i++)
    {
        if ((sol_x[i] >= 0) && (sol_x[i] <= 1) && (sol_y[i] >= 0) && (sol_y[i] <= 1))
        {
            o[0] =  sol_x[i];
            o[1] =  sol_y[i];
            //printf("original = %lf, %lf \n", o[0],o[1]);
        }
        //printf("original incorrect.\n");
    }
}

//Computes the eigenvectors of a 2x2 tensor
void eigen_vector(icMatrix2x2 tensor, double evect[2])
{
    double a, b, c;
    if (tensor.entry[0][1] != tensor.entry[1][0])
    {
        printf("Not a symmetric tensor!\n");
    }
    a = tensor.entry[0][0];
    b = tensor.entry[1][0];
    c = tensor.entry[1][1];
    
    double lambda, theta;
    lambda = (a + c)/2.0;
    theta = atan(2 * b / (a - c))/2.0;
    evect[0] = lambda*cos(theta); // (for minor eigenvector add pi/2 to theta)
    evect[1] = lambda*sin(theta);
}

//even faster way?: http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
//An old way of computing eigenvectors
//void eigen_vector(icMatrix2x2 tensor, double ev[2]){
//    double a = tensor.entry[0][0];
//    double b = tensor.entry[0][1];
//    double c = tensor.entry[1][1];
//    
//    double d = a*a + 4.00*b*b -2.00*a*c + c*c;
//    
//    //change to -1.0 or +1.0 to switch between eigenvectors
//    double sign = 1.00;
//    
//    double x = -(-a+c+(sign*sqrt(d)))/(2.00*b);
//    double y = 1.00;
//    
//    //double eval = 0.5*(a+c+(sign*sqrt(d)));
//    
//    ev[0] = x/(sqrt(x*x + y*y));
//    ev[1] = y/(sqrt(x*x + y*y));
//    return;
//}

//draws the entire vector field
//input is the four tensors defined on the corners of a unit square
void draw_vector_field(icMatrix2x2 tensor[4])
{
    int numVects = 25;//sets a numVectsXnumVects grid need to adjust spacing if this changes
    double tx, ty;
    icMatrix2x2 temp;
    for (int i = 0; i <= numVects; i++)
    {
        for (int j = 0; j <= numVects; j++)
        {
            tx = i * 0.04;
            ty = j * 0.04;
            temp = bilinear_interpolate(tensor, tx, ty);
            draw_vector(temp, tx, ty);
        }
    }
}

//draws a single eigenvector at (x,y)
void draw_vector(icMatrix2x2 t, double x, double y)
{
    double ev[2],mag=0.00000004;
    eigen_vector(t, ev);
    glBegin(GL_LINE_STRIP);
    glColor3f(0.1, 0.3, 0.7);
    glVertex2f(x, y);
    glVertex2f(x+ev[0]*mag, y+ev[1]*mag);
    glEnd();
}

// computes the roots of a cubic polynomial of the form
// Ax^3 + Bx^2 + Cx + D = 0
// returns the singularity type (if there is one real root then it is a wedge otherwise it is a trisector)
// explanation of algorithm: http://mathworld.wolfram.com/CubicFormula.html
int cubic_roots(double A, double B, double C, double D, double root[3])
{
    int realRoots;
    
    B=B/A;
    C=C/A;
    D=D/A;
    dcomp x1(0.0,0.0);
    dcomp x2(0.0,0.0);
    dcomp x3(0.0,0.0);
    dcomp Q(((1.0/9.0)*(3.0*C-(B*B))),0.0);
    dcomp R((1.0/54.0)*(-2.0*(B*B*B) + 9.0*B*C-27.0*D),0.0);
    dcomp DD = Q*Q*Q + R*R;
    dcomp T = pow((R-sqrt(DD)),1.0/3.0);
    dcomp S = pow((sqrt(DD)+R),1.0/3.0);
    dcomp I(0.0,1.0);
    dcomp F,G,H;
    //double discriminant;
    //discriminant = B*B*C*C - 4*A*C*C*C - 4*B*B*B*D - 27*A*A*D*D + 18*A*B*C*D;
   
    //3 roots
    x1 = (S+T)-(B/3.0);
    x2 = -(B/3.0) + ((0.5)*(sqrt(3.0)*I*(S-T))) - ((S+T)/2.0);
    x3 = -(B/3.0) - ((0.5)*(sqrt(3.0)*I*(S-T))) - ((S+T)/2.0);
    
    //for tensor fields we only care about the real roots
    root[0] = x1.real();
    root[1] = x2.real();
    root[2] = x3.real();
    
    //account for floating point errors and count the number of real roots
    double precision = 0.00000001;
    if(x1.imag()<precision){
        realRoots++;
    }
    if(x2.imag()<precision){
        realRoots++;
    }
    if(x3.imag()<precision){
        realRoots++;
    }
    if(realRoots==1){
        return WEDGE;
    }
    return TRISECTOR;
}


//For defining tensor components with predefined functions
// matrix form:
// {{a b},
// {b c}}
double a(double x){
    return -2.0*x+1;
}
double b(double y){
    return y;
}
double c(double x){
    return x*x;
}


void display(void)
{
	initGL();
	double theta[3]; // separatrics angle 
	double root[3]; // the roots of the cubic equation from Delmarcelle
	double A, B, C, D, px, py;
	int type = 1;
    double singularity[2];
    icMatrix2x2 tensor[4];
    

   
//Symbolic input
//if the grid starts at (0,0) and goes to (1,1)
//    tensor[0].set(a(1.0),b(1.0),b(1.0),c(1.0));
//    tensor[1].set(a(0.0),b(1.0),b(1.0),c(0.0));
//    tensor[2].set(a(0.0),b(0.0),b(0.0),c(0.0));
//    tensor[3].set(a(1.0),b(0.0),b(0.0),c(1.0));

//Some numerical inputs (data from Abaqus):
    
    tensor[0].set(-486250,-13189,-13189,-524850);
    tensor[1].set(-501290,-2951,-2951,-466510);
    tensor[2].set(-511890,11300,11300,-472760);
    tensor[3].set(-495640,2649.20,2649.20,-533430);
    
//    tensor[0].set(-86962,-12453,-12453,-118940);
//    tensor[1].set(-72879,7083,7083,-99363);
//    tensor[2].set(-122820,297.65,297.65,-114860);
//    tensor[3].set(-135660,-19724,-19724,-136950);

    
//    tensor[0].set(-72879,-7083,-7083,-99363);
//    tensor[1].set(-86962,12453,12453,-118940);
//    tensor[2].set(-135660,19724,19724,19724);
//    tensor[3].set(-122820,-297.65,-297.65,-114860);

//for testing interpolation
//    double xtest=0.5, ytest=0.5;
//    icMatrix2x2 newTensor;
//    newTensor=bilinear_interpolate(tensor,xtest,ytest);

    get_singularity(tensor, singularity);
    
	glPushMatrix();
	glClearColor(1.0, 1.0, 1.0, 1.0);  // background for rendering color coding 
	glClear(GL_COLOR_BUFFER_BIT);

	//draws outer quad
	glBegin(GL_LINE_STRIP);
	glColor3f(1.0, 1.0, 1.0);
	glVertex2f(0.0, 0.0);
	glVertex2f(0.0, 1.0);
	glVertex2f(1.0, 1.0);
	glVertex2f(1.0, 0.0);
	glVertex2f(0.0, 0.0);
	glEnd();
	
    glEnable(GL_BLEND);
    draw_vector_field(tensor);
    glDisable(GL_BLEND);
	
    double x = singularity[0];
    double y = singularity[1];
    
    //These values come from Delmarcelle
    A= 0.5*(-tensor[2].entry[0][0]*(1-y)+tensor[3].entry[0][0]*(1-y)-tensor[1].entry[0][0]*y+tensor[0].entry[0][0]*y + tensor[2].entry[1][1]*(1-y)-tensor[3].entry[1][1]*(1-y)+tensor[1].entry[1][1]*y-tensor[0].entry[1][1]*y);
    
    B= 0.5*(-tensor[2].entry[0][0]*(1-x) - tensor[3].entry[0][0]*x + tensor[1].entry[0][0]*(1-x)+tensor[0].entry[0][0]*x +tensor[2].entry[1][1]*(1-x) + tensor[3].entry[1][1]*x + tensor[1].entry[1][1]*(1-x)-tensor[0].entry[1][1]*x);
    
    C= -tensor[2].entry[0][1]*(1-y) + tensor[3].entry[0][1]*(1-y)-tensor[1].entry[0][1]*y+tensor[0].entry[0][1]*y;
    
    D= -tensor[2].entry[0][1]*(1-x) - tensor[3].entry[0][1]*x + tensor[1].entry[0][1]*(1-x)+tensor[0].entry[0][1]*x;


    type = cubic_roots(D, (C+2*B), (2*A-D), -C, root);
    //type = cubic_roots(A,B,C,D,root);
    
    if (type == WEDGE) // only one root
	{
		//glEnable(GL_BLEND);
		
        theta[0] = atan(root[0]);
		px = 1.0*cos(theta[0]);
		py = 1.0*sin(theta[0]);

		glBegin(GL_LINE_STRIP);
		glColor3f(0.5, 0.5, 0);
		glVertex2f(singularity[0], singularity[1]);// all separatrices start from singularity point
        glVertex2f(singularity[0]+px, singularity[1]+py);
		glEnd();
		//glDisable(GL_BLEND);
	}
    
    //look at page 123 figure 4.9 of dissertation
	else if (type == TRISECTOR) // three roots
	{
        //trisectors have three separatrices
		theta[0] = atan(root[0]);
        theta[1] = atan(root[1]);
        theta[2] = atan(root[2]);
        
        
		px = 1.0*cos(theta[0]);
		py = 1.0*sin(theta[0]);

		glBegin(GL_LINE_STRIP);
		glColor3f(0.5, 0.5, 0);
		glVertex2f(singularity[0], singularity[1]);
		glVertex2f(singularity[0]+px, singularity[1]+py);
		glEnd();

		
		px = 1.0*cos(theta[1]);
		py = 1.0*sin(theta[1]);

		glBegin(GL_LINE_STRIP);
		glColor3f(0.5, 0.5, 0);
		glVertex2f(singularity[0], singularity[1]);
        glVertex2f(singularity[0]+px, singularity[1]+py);
		glEnd();
		
		px = 1.0*cos(theta[2]);
		py = 1.0*sin(theta[2]);

		glBegin(GL_LINE_STRIP);
		glColor3f(0.5, 0.5, 0);
		glVertex2f(singularity[0], singularity[1]);
        glVertex2f(singularity[0]+px, singularity[1]+py);
		glEnd();
        
        //left side
        px = 1.0*cos(theta[0]);
        py = 1.0*sin(theta[0]);
        
        glBegin(GL_LINE_STRIP);
        glColor3f(1.0, 0, 0);
        glVertex2f(singularity[0], singularity[1]);
        glVertex2f(singularity[0]-px, singularity[1]-py);
        glEnd();
        
        
        px = 1.0*cos(theta[1]);
        py = 1.0*sin(theta[1]);
        
        glBegin(GL_LINE_STRIP);
        glColor3f(1.0, 0, 0);
        glVertex2f(singularity[0], singularity[1]);
        glVertex2f(singularity[0]-px, singularity[1]-py);
        glEnd();
        
        px = 1.0*cos(theta[2]);
        py = 1.0*sin(theta[2]);
        
        glBegin(GL_LINE_STRIP);
        glColor3f(1.0, 0, 0);
        glVertex2f(singularity[0], singularity[1]);
        glVertex2f(singularity[0]-px, singularity[1]-py);
        glEnd();

	}
	glPopMatrix();
	glFlush();
	glutSwapBuffers();
	glFinish();

}

//not currently used
//for future keyboard input
void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 27: //esc
		exit(0);
		break;
	}//switch key
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(WIN_SIZE, WIN_SIZE);
	glutCreateWindow(argv[0]);
	glutDisplayFunc(display);
	glutIdleFunc(display);
    //glutkeyboardfunc(keyboard)
	glutMainLoop();
	return 1;
}