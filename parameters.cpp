#include <cmath>
#include <math.h>

double x_start = 0.0;
double x_end = 1.0;
double y_start = 0.0;
double y_end = 1.0;

int LDC = 1;
int TG = 0;
int IBM = 1;

long expnt = 6;
long nx_int = pow(2,expnt);
long ny_int = pow(2,expnt);
long nx = nx_int;
long ny = ny_int;
long n = nx;

double radius = 0.15;
double ibx = 0.7;
double iby = 0.5;
double ELrat = 0.8;
long alpha = ceil(2*M_PI*radius*n/ELrat);

double Tfinal = 50.0;
double dx = (x_end-x_start)/nx_int;
double dy = (y_end-y_start)/ny_int;
double h = dx;
double s = (2.0*M_PI*radius)/alpha;

double CFL = 0.8;
double dt = h*CFL;
double tol = 10e-6;

double re = 1000.0;
double nu = 1.0/re;

char DirName[500] = "";