#ifndef __Multigrid_H__
#define __Multigrid_H__

#include "LinearAlgebra.h"

struct LvlData
{
    SparseMatrix* R;
    SparseMatrix* I;
    SparseMatrix* M;
    long* n;
};

class Multigrid
{
    public:    
    Multigrid(long n,long l);
    void initialize(long nn);
    void solve(DoubleArray1D f,DoubleArray1D& v);
    void vCyc(DoubleArray1D f, DoubleArray1D& v);
    void Setupu();
    void Setupv();
    void Setupp();


    LvlData Grid;

    protected:
    //solution functions
    void GSm(SparseMatrix A,DoubleArray1D f,DoubleArray1D& v,long iter);
    void vCycle(DoubleArray1D& f,DoubleArray1D& v,long m);
    DoubleArray1D findresidue(DoubleArray1D f,DoubleArray1D v,long m);
    void FMG(DoubleArray1D& f,DoubleArray1D& v,long m);

    //restriction/interpolation matrix construction
    SparseMatrix restrictionB(long n);
    void restrictionK(SparseMatrix& A,long n);

    long max_lvls;
};

#endif