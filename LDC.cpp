#include <cmath>
#include <iostream>
#include <omp.h>
#include <fstream>
#include <cstring>
#include <chrono>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "matrices.h"
#include "LinearAlgebra.h"
#include "parameters.h"
#include "simulation.h"
#include "multigrid.h"
#include "gridtransfer.h"
#include "output.h"

using namespace std;

int main(int argc, char* argv[])
{   
    //read name of directory and create it
    strcpy(DirName,argv[1]);
    mkdir(DirName);

    auto start = chrono::high_resolution_clock::now();

    // initialize solution arrays
    DoubleArray1D u(n*(n+1));
    DoubleArray1D v(n*(n+1));
    DoubleArray1D u_old(n*(n+1));
    DoubleArray1D v_old(n*(n+1));

    u.setToValue(0.0);
    v.setToValue(0.0);
    u_old = u;
    v_old = v;


    // initialize finite difference matrices

    SparseMatrix UV(n*(n+1),n*(n+1));
    SparseMatrix VU(n*(n+1),n*(n+1));

    UtoV(UV,n);
    VtoU(VU,n);

    SparseMatrix DxU(n*(n+1),n*(n+1));
    SparseMatrix DyU(n*(n+1),n*(n+1));

    GradxU(DxU,n);
    GradyU(DyU,n);

    SparseMatrix DxV(n*(n+1),n*(n+1));
    SparseMatrix DyV(n*(n+1),n*(n+1));

    GradxV(DxV,n);
    GradyV(DyV,n);

    SparseMatrix DxPU(n*(n+1),n*n);
    SparseMatrix DyPV(n*(n+1),n*n);

    GradxPtoU(DxPU,n);
    GradyPtoV(DyPV,n);

    SparseMatrix DxUP(n*n,n*(n+1));
    SparseMatrix DyVP(n*n,n*(n+1));

    GradxUtoP(DxUP,n);
    GradyVtoP(DyVP,n);

    SparseMatrix DDU(n*(n+1),n*(n+1));
    SparseMatrix DDV(n*(n+1),n*(n+1));

    LaplaceU(DDU,n);
    LaplaceV(DDV,n); 

    //Initialize distribution matrix for immersed boundary method
    SparseMatrix Dhu(alpha,n*(n+1));
    SparseMatrix Dhv(alpha,n*(n+1));
    double rx,ry;

    if (IBM == 1)
    {

        for (int k = 0;k<alpha;k++)
        {
            for(int i = 0;i<n+1;i++)
            {
                for(int j = 0;j<n;j++)
                {
                    ry = ((j+0.5)*h - (iby + radius*sin((2.0*M_PI*k)/alpha)))/h;
                    if (abs(ry) >= 2.0) continue;
                    rx = (i*h - (ibx + radius*cos((2.0*M_PI*k)/alpha)))/h;

                    if (abs(rx) >= 2.0) continue;
                    if (rx == 0.0 && ry == 0.0)
                    {
                        Dhu.resetRow(k);
                        Dhu.addEntry(k,i+(n+1)*j,1.0);
                        goto cont1;
                    }
                    Dhu.addEntry(k,i+(n+1)*j,(1.0/(h*h))*DistnF(rx)*DistnF(ry));
                }
            }
            cont1:;

            for(int i = 0;i<n;i++)
            {
                for(int j = 0;j<n+1;j++)
                {
                    ry = (j*h - (iby + radius*sin((2.0*M_PI*k)/alpha)))/h;
                    if (abs(ry) >= 2.0) continue;
                    rx = ((i+0.5)*h - (ibx + radius*cos((2.0*M_PI*k)/alpha)))/h;
                    if (abs(rx) >= 2.0) continue;
                    if (rx == 0.0 && ry == 0.0)
                    {
                        Dhv.resetRow(k);
                        Dhv.addEntry(k,i+n*j,1.0);
                        goto cont2;
                    }
                    Dhv.addEntry(k,i+n*j,(1.0/(h*h))*DistnF(rx)*DistnF(ry));
                }
            }
            cont2:;
        }
    }

    auto stop = chrono::high_resolution_clock::now();
    auto dur = chrono::duration_cast<chrono::milliseconds>(stop-start);
    double len = dur.count();
    printf("Initialize Arrays: %.0f ms\n",len);
    
    start = chrono::high_resolution_clock::now();

    Multigrid D_phi(n,expnt);
    Multigrid D_star_u(n,expnt);
    Multigrid D_star_v(n,expnt);
    
    D_phi.Setupp();
    D_star_u.Setupu();
    D_star_v.Setupv();

    stop = chrono::high_resolution_clock::now();
    dur = chrono::duration_cast<chrono::milliseconds>(stop-start);
    len = dur.count();
    printf("Initialize Grid: %.0f ms\n",len);

    simulateTime(u,v,D_star_u,D_star_v,D_phi,u_old,v_old,UV,VU,DxU,DyU,DxV,DyV,DDU,DDV,DxPU,DyPV,DxUP,DyVP,Dhu,Dhv);
    return 0;
}