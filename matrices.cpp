#include <cmath>
#include <omp.h>
#include "LinearAlgebra.h"
#include "parameters.h"

void UtoV(SparseMatrix& A,int n) //output on v grid
{
    //need to enforce boundary conditions for dirichlet problem after matrix multiplication
    int inr, i;
    for (i =0; i<n*(n+1);i++)
    {
        inr = (i/n)-1;
        if (i < n) //bottom
        {
            
        }
        else if (i >= n*n) //top
        {

        }
        else
        {
            A.addEntry(i,i+1+inr,0.25);
            A.addEntry(i,i+2+inr,0.25);
            A.addEntry(i,i-n+1+inr,0.25);
            A.addEntry(i,i-n+inr,0.25);
        }
        
    }
}

void VtoU(SparseMatrix& A,int n) //output on u grid
{
    //need to enforce boundary conditions for dirichlet problem after matrix multiplication
    int inr, i;
    for (i =0; i<n*(n+1);i++)
    {
        inr = i/(n+1);
        if (i%(n+1) == 0) //left
        {

        }
        else if (i%(n+1) == n) //right
        {
            
        }
        else
        {
            A.addEntry(i,i-inr,0.25);
            A.addEntry(i,i-inr-1,0.25);
            A.addEntry(i,i-inr+n,0.25);
            A.addEntry(i,i-inr-1+n,0.25);
        }
        
    }
}

void UtoP(SparseMatrix& A,int n) //output on p grid
{
    int inr;
    for (int i = 0;i<n*n;i++)
    {
        inr = i/n;
        A.addEntry(i,inr+i,0.5);
        A.addEntry(i,inr+i+1,0.5);
        
    }
}

void VtoP(SparseMatrix& A,int n) //output on p grid
{
    for (int i = 0;i<n*n;i++)
    {
        A.addEntry(i,i,0.5);
        A.addEntry(i,i+n,0.5);        
    }
}

void GradxU(SparseMatrix &A,int n) // output on u grid
{
    //need to enforce boundary conditions for dirichlet problem after matrix multiplication
    for (int i =0; i<n*(n+1);i++)
    {
        if (i%(n+1) == 0) //left
        {   
            
        }
        else if (i%(n+1) == n) //right
        {
            
        }
        else
        {
            A.addEntry(i,i+1,1.0);
            A.addEntry(i,i-1,-1.0);
        }
        
    }
    A = A/(2.0*h);
}

void GradyU(SparseMatrix& A,int n) // output on u grid
{
    //need to enforce boundary conditions for dirichlet problem after matrix multiplication
    for (int i =0; i<n*(n+1);i++)
    {
        if (i < n+1) //bottom
        {
            A.addEntry(i,i,2.0);
            A.addEntry(i,i+(n+1),2.0/3.0);            
        }
        else if (i >= (n-1)*(n+1)) //top
        {
            A.addEntry(i,i,-2.0);
            A.addEntry(i,i-(n+1),-2.0/3.0);
        }
        else
        {
            A.addEntry(i,i+(n+1),1.0);
            A.addEntry(i,i-(n+1),-1.0);
        }
        
    }
    A = A/(2.0*h);
}

void GradxV(SparseMatrix& A,int n) //output on v grid
{
    //need to enforce boundary conditions for dirichlet problem after matrix multiplication
    for (int i =0; i<n*(n+1);i++)
    {
        if (i%n == 0) //left
        {   
            A.addEntry(i,i,2.0);
            A.addEntry(i,i+1,2.0/3.0);
        }
        else if (i%n == n-1) //right
        {
            A.addEntry(i,i,-2.0);
            A.addEntry(i,i-1,-2.0/3.0);
        }
        else
        {
            A.addEntry(i,i+1,1.0);
            A.addEntry(i,i-1,-1.0);
            
        }
        
    }
    A = A/(2.0*h);

}
void GradyV(SparseMatrix& A,int n) // output on v grid
{
    //need to enforce boundary conditions for dirichlet problem after matrix multiplication
    for (int i =0; i<n*(n+1);i++)
    {
        if (i < n) //bottom
        {
            
        }
        else if (i >= n*n) //top
        {
            
        }
        else
        {
            A.addEntry(i,i+n,1.0);
            A.addEntry(i,i-n,-1.0);
        }
        
    }
    A = A/(2.0*h);

}

void GradxPtoU(SparseMatrix& A,int n) //output on u grid
{
    //need to enforce boundary conditions for dirichlet problem after matrix multiplication
    int inr;
    for (int i =0; i<n*(n+1);i++)
    {
        inr = i/(n+1);
        if (i%(n+1) == 0) //left
        {
            
        }
        else if (i%(n+1) == n) //right
        {

        }
        else
        {
            A.addEntry(i,i-inr-1,-1.0);
            A.addEntry(i,i-inr,1.0);
        }
        
    }
    A = A/h;
}

void GradyPtoV(SparseMatrix& A,int n) //output on v grid
{
    for (int i =0; i<n*(n+1);i++)
    {
        if (i < n) //bottom
        {
        //     if (Nmpt == 0)
        //     {
        //         A.addEntry(i,i,-2.0);
        //     }
            
        }
        else if (i >= n*n) //top
        {
            // if (Nmpb == 0)
            // {
            //     A.addEntry(i,i-n,2.0);
            // }
        }
        else
        {
            A.addEntry(i,i,1.0);
            A.addEntry(i,i-n,-1.0);
        }
        
    }
    A = A/h;
}

void GradxUtoP(SparseMatrix& A, int n) //output on p grid
{
    int inr;
    for (int i = 0;i<n*n;i++)
    {
        inr = i/n;
        A.addEntry(i,i+inr,-1.0);
        A.addEntry(i,i+inr+1,1.0);
        
    }
    A = A/h;
}

void GradyVtoP(SparseMatrix& A, int n) //output on p grid
{
    for (int i = 0;i<n*n;i++)
    {
        A.addEntry(i,i,-1.0);
        A.addEntry(i,i+n,1.0);
        
    } 
    A = A/h;
}

void LaplaceU(SparseMatrix& A, int n)
{
    for (int i = 0;i<n*(n+1);i++)
    {
        if (i%(n+1) == 0 || i%(n+1) == n) //left or right
        {

        }
        else
        {
            A.addEntry(i,i-1,1.0);
            A.addEntry(i,i+1,1.0);
            if (i > 0 && i < n) //interior bottom
            {
                A.addEntry(i,i,-6.0);
                A.addEntry(i,i+n+1,4.0/3.0);
            }
            else if (i>(n-1)*(n+1) && i < (n*(n+1) - 1)) //interior top
            {
                A.addEntry(i,i,-6.0);
                A.addEntry(i,i-n-1,4.0/3.0);
            }
            else
            {
                A.addEntry(i,i,-4.0);
                A.addEntry(i,i+n+1,1.0);
                A.addEntry(i,i-n-1,1.0);
            }
        }
    }
    A = A/pow(h,2.0);

}

void LaplaceV(SparseMatrix& A, int n)
{
    for (int i = 0;i<n*(n+1);i++)
    {
        if (i < n || i > n*n-1) //bottom or top
        {

        }
        else
        {
            A.addEntry(i,i-n,1.0);
            A.addEntry(i,i+n,1.0);
            if (i != 0 && i != n*n && i%n == 0) //interior left
            {
                A.addEntry(i,i,-6.0);
                A.addEntry(i,i+1,4.0/3.0);
            }
            else if (i != n-1 && i != n*(n+1)-1 && i%n == n-1) //interior right
            {
                A.addEntry(i,i,-6.0);
                A.addEntry(i,i-1,4.0/3.0);
            }
            else
            {
                A.addEntry(i,i,-4.0);
                A.addEntry(i,i+1,1.0);
                A.addEntry(i,i-1,1.0);
            }
        }
    }
    A = A/pow(h,2.0);
    
}

void createDstaru(SparseMatrix& A, int n)
{
    double hi = (x_end-x_start)/n;
    for (int i = 0;i<n*(n+1);i++)
    {
        if (i%(n+1) == 0 || i%(n+1) == n) //left or right
        {
            A.addEntry(i,i,1.0);
        }
        else
        {
            A.addEntry(i,i-1,1.0);
            A.addEntry(i,i+1,1.0);
            if (i > 0 && i < n) // interior bottom
            {
                A.addEntry(i,i,-6.0 - 2.0*pow(hi,2.0)/(dt*nu));
                A.addEntry(i,i+n+1,4.0/3.0);
            }
            else if (i>(n-1)*(n+1)-1 && i < (n*(n+1) - 1)) //interior top
            {
                A.addEntry(i,i,-6.0 - 2.0*pow(hi,2.0)/(dt*nu));
                A.addEntry(i,i-n-1,4.0/3.0);
            }
            else
            {
                A.addEntry(i,i,-4.0 - 2.0*pow(hi,2.0)/(dt*nu));
                A.addEntry(i,i+n+1,1.0);
                A.addEntry(i,i-n-1,1.0);
            }
        }        
    }
    A = A/pow(hi,2.0);
}

void createDstarv(SparseMatrix& A, int n)
{
    double hi = (x_end-x_start)/n;
    for (int i = 0;i<n*(n+1);i++)
    {
        if (i < n || i > n*n-1) //bottom or top
        {
            A.addEntry(i,i,1.0);
        }
        else
        {
            A.addEntry(i,i-n,1.0);
            A.addEntry(i,i+n,1.0);
            if (i%n == 0) //interior left
            {
                A.addEntry(i,i,-6.0 - 2.0*pow(hi,2.0)/(dt*nu));
                A.addEntry(i,i+1,4.0/3.0);
            }
            else if (i%n == n-1) //interior right
            {
                A.addEntry(i,i,-6.0 - 2.0*pow(hi,2.0)/(dt*nu));
                A.addEntry(i,i-1,4.0/3.0);
            }
            else
            {
                A.addEntry(i,i,-4.0 - 2.0*pow(hi,2.0)/(dt*nu));
                A.addEntry(i,i+1,1.0);
                A.addEntry(i,i-1,1.0);
            }
        }

    }
    A = A/pow(hi,2.0);
}

void createDphi(SparseMatrix& A, int n)
{
    double hi = (x_end-x_start)/n;
    // nx and ny are the number of nodes not intervals
    for (int i = 0;i<n*n;i++){
        A.addEntry(i,i,-4.0);
        if (i < n*n-1)
        {
            A.addEntry(i,i+1,1.0);
        }
        if (i > 0)
        {
            A.addEntry(i,i-1,1.0);
        }
        if (i < n*n-n)
        {
            A.addEntry(i,i+n,1.0);
        }
        if (i > n-1)
        {
            A.addEntry(i,i-n,1.0);
        }

        if (i%n == 0)
        {//left
            A.addEntry(i,i,-3.0);
            if (i == 0)
            {
                A.addEntry(i,i+1,1.0);
            }
            else
            {
                A.addEntry(i,i-1,0.0);
                A.addEntry(i,i+1,1.0);
            }
        }
        if (i%n == n-1)
        {//right
            A.addEntry(i,i,-3.0);
            if (i == n*n-1)
            {
                A.addEntry(i,i-1,1.0);
            }
            else
            {
                A.addEntry(i,i+1,0.0);
                A.addEntry(i,i-1,1.0);
            }
        }
        if (i >= n*n-n)
        {//top
            A.addEntry(i,i,-3.0);
            A.addEntry(i,i-n,1.0);
        }
        if (i < n)
        {//bottom
            A.addEntry(i,i,-3.0);
            A.addEntry(i,i+n,1.0);
        }
        
        if (i == 0 || i == n-1 || i == (n-1)*n || i == n*n-1)
        {
            A.addEntry(i,i,-2.0);            
        }
    }
    A = A/pow(hi,2.0);

}

void GradxC(SparseMatrix& A, int n)
{
    double hi = (x_end-x_start)/n;
    for (int i = 0;i < n*n;i++)
    {
        
    }
    A = A/(2.0*hi);
}

void GradyC(SparseMatrix& A, int n)
{
    double hi = (x_end-x_start)/n;
    for (int i = 0;i < n*n;i++)
    {


        
    }
    A = A/(2.0*hi);
}

void LaplaceC(SparseMatrix& A, int n)
{
    double hi = (x_end-x_start)/n;
    for (int i = 0;i < n*n;i++)
    {
        
    }
    A = A/(2.0*hi);
}

void createCm(SparseMatrix& A, int n)
{
    double hi = (x_end-x_start)/n;
    for (int i = 0;i < n*n;i++)
    {
        
    }
    A = A/pow(hi,2.0);
}

double f1(double r)
{
    return (3.0 - 2.0*abs(r) + sqrt(1.0+4.0*abs(r)-4.0*pow(abs(r),2.0)))/8.0;
}

double DistnF(double r)
{
    if (abs(r) > 1)
    {
        return 0.5 - f1(2.0-abs(r));
    }
    else
    {
        return f1(r);
    }
}