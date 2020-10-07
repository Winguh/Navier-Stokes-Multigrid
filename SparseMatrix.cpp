#include <iostream>
#include <iomanip>
#include "LinearAlgebra.h"
using namespace std;

//
//###################################################################
//                 Constructors/Initialization
//###################################################################
//
SparseMatrix::SparseMatrix()
{
    size1 = 0;
    size2 = 0;
    rows = NULL;
}

SparseMatrix::SparseMatrix(long m, long n)
{
    initialize(m,n);
}

SparseMatrix::SparseMatrix(const SparseMatrix& D)
{
    size1 = D.size1;
    size2 = D.size2;
    rows = new SparseRow[size1];
    for (int i = 0; i < size1; i++)
    {
        rows[i] = D.rows[i];
    }
}

SparseMatrix::~SparseMatrix()
{
    if( rows != NULL)
    {
        delete [] rows;
    }
}

void SparseMatrix::initialize(long m, long n)
{
    size1 = m;
    size2 = n;
    rows = new SparseRow[m];
    for (int i = 0; i < m; i++)
    {
        rows[i] = SparseRow(n);
    }
}


//
//###################################################################
//                  Element Access
//###################################################################
//
    
double SparseMatrix::operator()(long i1, long i2)
{
    return rows[i1](i2);
}

SparseRow&  SparseMatrix::operator()(long i1)
{
    return rows[i1];
}

const SparseRow&  SparseMatrix::operator()(long i1) const
{
    return rows[i1];
}

//
//###################################################################
//                Array Structure Access Functions
//###################################################################
//

void SparseMatrix::addEntry(long i, long j, double v)
{
    rows[i].addEntry(j,v);
}
    
//
//###################################################################
//                     Array Operators
//###################################################################
//

SparseMatrix SparseMatrix::operator+(const SparseMatrix& D)
{
    SparseMatrix M(*this);
    int i;
    for (i = 0; i < D.size1; i++)
    {
        M.rows[i] += D.rows[i];
    }
    return M;
}

SparseMatrix SparseMatrix::operator-(const SparseMatrix& D)
{
    SparseMatrix M(*this);
    int i;
    for (i = 0; i < D.size1; i++)
    {
        M.rows[i] -= D.rows[i];
    }
    return M;
}

SparseMatrix SparseMatrix::operator*(double alpha)
{
    SparseMatrix M(*this);
    int i;
    for (i = 0; i < size1; i++)
    {
        M.rows[i] *= alpha;
    }
    return M;
}

SparseMatrix operator*(double alpha, const SparseMatrix& D)
{
    SparseMatrix M(D);
    M *= alpha;
    return M;
}

SparseMatrix SparseMatrix::operator/(double alpha)
{
    SparseMatrix M(*this);
    M *= 1.0/alpha;
    return M;
}

void SparseMatrix::operator=(const SparseMatrix& D)
{
    if (size1 != D.size1 || size2 != D.size2)
    {
        size1 = D.size1;
        size2 = D.size2;
        if( rows != NULL)
        {
            delete [] rows;
        }
    }
    rows = new SparseRow[size1];
    int i;
    for (i = 0; i < size1; i++)
    {
        rows[i] = D.rows[i];
    }
}

void SparseMatrix::operator*=(double alpha)
{
    int i;
    for (i = 0; i < size1; i++)
    {
        rows[i] *= alpha;
    }
}

void SparseMatrix::operator+=(const SparseMatrix& D)
{
    int i;
    for (i = 0; i < D.size1; i++)
    {
        rows[i] += D.rows[i];
    }
}

void SparseMatrix::operator-=(const SparseMatrix& D)
{
    int i;
    for (i = 0; i < D.size1; i++)
    {
        rows[i] -= D.rows[i];
    }
}

DoubleArray1D SparseMatrix::operator*(const DoubleArray1D& b) const
{
    DoubleArray1D v(size1);
    int i;
    for (i = 0; i < size1; i++)
    {
        v(i) = rows[i].dot(b);
    }
    return v;
}

SparseMatrix SparseMatrix::operator*(const SparseMatrix M) const
{
    SparseMatrix A(size1,M.getSize2());
    double* val;
    long* ind;
    long nnz;
    double inter;
    for (int i = 0;i<size1;i++)
    {   
        for (int j = 0;j<M.getSize2();j++)
        {
            val = rows[i].getValue();
            ind = rows[i].getIndex();
            nnz = rows[i].getNnz();
            inter = 0.0;
            for (int k = 0;k<nnz;k++)
            {
                if (M.rows[ind[k]](j) != 0.0)
                {
                    inter += val[k]*M.rows[ind[k]](j);           
                }
            }
            if (inter != 0.0)
            {
                A.addEntry(i,j,inter);
            }
        }
    }
    return A;
}

void SparseMatrix::resetRow(long i)
{
    rows[i] = SparseRow(size2);
}

DoubleArray1D SparseMatrix::transMultiply(const DoubleArray1D& b) const
{
    DoubleArray1D v(size2);
    double* value;
    long* index;
    long nnz;
    for (int i = 0; i < size1; i++)
    {
        value = rows[i].getValue();
        index = rows[i].getIndex();
        nnz = rows[i].getNnz();
        for (int j = 0; j < nnz; j++)
        {
            v(index[j]) += value[j]*b(i);
        }
    }
    return v;
}

SparseMatrix SparseMatrix::transpose() const
{
    SparseMatrix B(size2,size1);
    double* val;
    long* ind;
    long nnz;
    for (int i=0;i<size1;i++)
    {
        val = rows[i].getValue();
        ind = rows[i].getIndex();
        nnz = rows[i].getNnz();
        for (int j=0;j<nnz;j++)
        {   
            B.addEntry(ind[j],i,val[j]);
        }
    }
    return B;
}

//
//  Output
//
ostream& operator<<(ostream& outStream, const SparseMatrix& V)
{
    for (int i = 0; i < V.size1; i++)
    {
        if (V.rows[i].getNnz()>0)
            outStream<<"Row "<<i<<" : "<< V.rows[i];
            // outStream<< V.rows[i];
    }
    outStream<<endl;
    return outStream;
}

