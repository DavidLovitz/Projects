//#include <glew.h>
//#include <math.h>
#include <GLUT/glut.h> //An easier to work with OpenGL API (OpenGL Utility Toolkit)
#include <stdlib.h>
#include <stdio.h>
#include <complex>//complex numbers
#include <cmath> //math on complex numbers
#include "icMatrix.h"
#include "icVector.h"



#define WIN_SIZE  512 //window widht&height in pixels

#define WEDGE			1
#define TRISECTOR		-1


void display(void);
void keyboard(unsigned char key, int x, int y);

int cubic_roots(double A, double B, double C, double D,double root[3]);

void eigen_vector(double tensor[2][2], double ev[2]);
void draw_vector(icMatrix2x2 t,double x,double y);
void draw_vector_field(icMatrix2x2 tensor[4]);
void get_singularity(icMatrix2x2 tensor[4], double o[2]);
void find_singularity(icMatrix2x2 tensor[4], double o[2]);
icMatrix2x2 interpolate(double tensor[4], double x, double y);

double a(double x);
double b(double y);
double c(double x);