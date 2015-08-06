#include <stdio.h>
#include <math.h>

int main()
{
    printf("Navier-Stokes Solver\n");

    const int nx = 10;
    const int ny = 10;
    const double c = 1500;
    const double dt = 1;
    const double dx = 1;
    const double dy  =1;

    const int SIZE = nx*ny;



    double *rho, *u, *v, *rho_star, *u_star, *v_star, *rho_1, *u_1, *v_1;

    rho = (double*) malloc(SIZE*sizeof(double));
    u = (double*) malloc(SIZE*sizeof(double));
    v = (double*) malloc(SIZE*sizeof(double));
    rho_star = (double*) malloc(SIZE*sizeof(double));
    u_star = (double*) malloc(SIZE*sizeof(double));
    v_star = (double*) malloc(SIZE*sizeof(double));
    rho_1 = (double*) malloc(SIZE*sizeof(double));
    u_1 = (double*) malloc(SIZE*sizeof(double));
    v_1 = (double*) malloc(SIZE*sizeof(double));
    
    
    //getchar();
    return 0;
}
