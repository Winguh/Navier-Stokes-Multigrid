#include <cmath>
#include <omp.h>
#include "LinearAlgebra.h"
// #include "Grid.h"
#include "multigrid.h"
#include "matrices.h"
#include "gridtransfer.h"

Multigrid::Multigrid(long n,long l)
{
    max_lvls = l;
    initialize(n);
};

void Multigrid::initialize(long nn)
{
    Grid.n = new long[max_lvls+1];
    Grid.R = new SparseMatrix[max_lvls];
    Grid.I = new SparseMatrix[max_lvls];
    Grid.M = new SparseMatrix[max_lvls];
    Grid.n[0] = nn;
};

void Multigrid::Setupu()
{
    for (int i=0;i<max_lvls;i++)
    {   
        Grid.n[i+1] = Grid.n[i]/2;
        SparseMatrix mu(Grid.n[i]*(Grid.n[i]+1),Grid.n[i]*(Grid.n[i]+1));
        createDstaru(mu,Grid.n[i]);
        Grid.M[i] = mu;
        SparseMatrix ru(Grid.n[i+1]*(Grid.n[i+1]+1),Grid.n[i]*(Grid.n[i]+1));
        SparseMatrix iu(Grid.n[i]*(Grid.n[i]+1),Grid.n[i+1]*(Grid.n[i+1]+1));
        RestrictU(ru,Grid.n[i]);
        Grid.R[i] = ru;
        ProlongateU(iu,Grid.n[i+1]);
        Grid.I[i] = iu;
    }
};

void Multigrid::Setupv()
{
    for (int i=0;i<max_lvls;i++)
    {   
        Grid.n[i+1] = Grid.n[i]/2;
        SparseMatrix mv(Grid.n[i]*(Grid.n[i]+1),Grid.n[i]*(Grid.n[i]+1));
        createDstarv(mv,Grid.n[i]);
        Grid.M[i] = mv;
        SparseMatrix rv(Grid.n[i+1]*(Grid.n[i+1]+1),Grid.n[i]*(Grid.n[i]+1));
        SparseMatrix iv(Grid.n[i]*(Grid.n[i]+1),Grid.n[i+1]*(Grid.n[i+1]+1));
        RestrictV(rv,Grid.n[i]);
        Grid.R[i] = rv;
        ProlongateV(iv,Grid.n[i+1]);
        Grid.I[i] = iv;
    }
};

void Multigrid::Setupp()
{
    for (int i=0;i<max_lvls;i++)
    {   
        Grid.n[i+1] = Grid.n[i]/2;
        SparseMatrix mm(Grid.n[i]*Grid.n[i],Grid.n[i]*Grid.n[i]);
        createDphi(mm,Grid.n[i]);
        Grid.M[i] = mm;
        SparseMatrix rr(Grid.n[i+1]*Grid.n[i+1],Grid.n[i]*Grid.n[i]);
        SparseMatrix ii(Grid.n[i]*Grid.n[i],Grid.n[i+1]*Grid.n[i+1]);
        RestrictP(rr,Grid.n[i]);
        Grid.R[i] = rr;
        ProlongateP(ii,Grid.n[i+1]);
        Grid.I[i] = ii;
    }
};

void Multigrid::solve( DoubleArray1D f, DoubleArray1D& v)
{
    FMG(f,v,0);
}

void Multigrid::vCyc(DoubleArray1D f, DoubleArray1D& v)
{
    vCycle(f,v,0);
}

void Multigrid::GSm(SparseMatrix A,DoubleArray1D f,DoubleArray1D& v,long iter)
{
    for (int j = 0; j < iter; j++)
    {
        for (int i = 0; i < A.getSize1(); i++)
        {
            // v(i) = (0.25)*v(i) + (0.75)*(f(i)-A.getRow(i).dot(v)+A(i,i)*v(i))/A(i,i);
            v(i) = (f(i)-A.getRow(i).dot(v)+A(i,i)*v(i))/A(i,i);
        }
    }
}

DoubleArray1D Multigrid::findresidue(DoubleArray1D f,DoubleArray1D v,long m)
{
    return f-Grid.M[m]*v;
}

void Multigrid::vCycle(DoubleArray1D& f,DoubleArray1D& v,long m)
{
    if (m == max_lvls-1)
    {   
        GSm(Grid.M[m],f,v,10);
    }
    else
    {
        GSm(Grid.M[m],f,v,5);
        DoubleArray1D f2 = Grid.R[m]*findresidue(f,v,m);
        DoubleArray1D v2(f2.getIndex1Size());
        vCycle(f2,v2,m+1);
        v = v + Grid.I[m]*v2;
        GSm(Grid.M[m],f,v,5);
    }
}

void Multigrid::FMG(DoubleArray1D& f,DoubleArray1D& v,long m)
{
    if (m == max_lvls-1)
    {
        v.setToValue(0.0);
        vCycle(f,v,m);
    }
    else
    {
        DoubleArray1D f2(Grid.R[m]*f);
        DoubleArray1D v2(f2.getIndex1Size());
        FMG(f2,v2,m+1);

        v = Grid.I[m]*v2;

        vCycle(f,v,m);
    }
}
