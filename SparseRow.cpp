#include <iostream>
#include <iomanip>
#include <cmath>
#include "LinearAlgebra.h"
using namespace std;

#ifdef  _DEBUG
#include <stdio.h>
#include <assert.h>
#endif


//
//###################################################################
//                 Constructors/Initialization
//###################################################################
//
SparseRow::SparseRow()
{
    initialize(0);
}

SparseRow::SparseRow(long m)
{
    initialize(m);
}

SparseRow::SparseRow(const SparseRow& D)
{
    size = D.size;
    nnz = D.nnz;
    index = new long[nnz];
    value = new double[nnz];
    for (long i = 0; i < nnz; i++)
    {
        index[i] = D.index[i];
        value[i] = D.value[i];
    }
}

SparseRow::~SparseRow()
{
    if( index != NULL) delete [] index;
    if( value != NULL) delete [] value;
}

void SparseRow::initialize(long m)
{
    size = m;
    nnz = 0;
    index = NULL;
    value = NULL;
}

    
//
//###################################################################
//                  Element Access
//###################################################################
//
    
double SparseRow::operator()(long i1)
{
    for (long i = 0; i < nnz; i++)
    {
        if (i1 == index[i])
        {
            return value[i];
        }
    }
    return 0.0;
}

//
//###################################################################
//                Array Structure Access Functions
//###################################################################
//
void SparseRow::addEntry(long i1, double v1)
{
    if (nnz==0)
    {
        if( index != NULL) delete [] index;
        if( value != NULL) delete [] value;
        index = new long[nnz+1];
        value = new double[nnz+1];
        index[0] = i1;
        value[0] = v1;
        nnz = 1;
    }
    else
    {
        bool found = 0;
        long j = 0;
        long* newIndex = new long[nnz+1];
        double* newValue = new double[nnz+1];
        for (long i = 0; i < nnz; i++)
        {
            if (i1 > index[i])
            {
                newValue[j] = value[i];
                newIndex[j] = index[i];
            }
            else
            {
                if (i1==index[i])
                {
                    newValue[j] = v1;
                    newIndex[j] = i1;
                }
                else
                {
                    if (found == 0)
                    {
                        newValue[j] = v1;
                        newIndex[j] = i1;
                        j++;
                    }
                    newValue[j] = value[i];
                    newIndex[j] = index[i];
                }
                found = 1;
            }
            j++;
        }
        if (found ==0)
        {
            newValue[j] = v1;
            newIndex[j] = i1;
            j++;
        }
        nnz = j;
        delete [] index;
        delete [] value;
        index = newIndex;
        value = newValue;
    }
}


//
//###################################################################
//                     Array Operators
//###################################################################
//

SparseRow SparseRow::operator+(const SparseRow& D)
{
#ifdef _DEBUG
    assert(size == D.size);
#endif
    SparseRow A(size);
    if (nnz==0)
    {
        A = SparseRow(D);
    }
    else
    {
        long iB = 0, iD = 0, i = 0;
        A.index = new long[nnz+D.nnz];
        A.value = new double[nnz+D.nnz];
        
        while (iB<nnz || iD<D.nnz)
        {
            if (iB<nnz && iD<D.nnz)
            {
                if (index[iB]<D.index[iD])
                {
                    A.index[i] = index[iB];
                    A.value[i] = value[iB];
                    iB++;
                }
                else if (index[iB]>D.index[iD])
                {
                    A.index[i] = D.index[iD];
                    A.value[i] = D.value[iD];
                    iD++;
                }
                else
                {
                    A.index[i] = D.index[iD];
                    A.value[i] = D.value[iD]+value[iB];
                    iD++;
                    iB++;
                }
            }
            else if (iB==nnz)
            {
                A.index[i] = D.index[iD];
                A.value[i] = D.value[iD];
                iD++;
            }
            else if (iD==D.nnz)
            {
                A.index[i] = index[iB];
                A.value[i] = value[iB];
                iB++;
            }
            i++;
        }
        A.nnz = i;
    }
    return A;
}

SparseRow SparseRow::operator-(const SparseRow& D)
{
    return (*this)+((-1.0)*D);
}

SparseRow SparseRow::operator*(double alpha)
{
    SparseRow A(*this);
    for (int i = 0; i < nnz; i++)
    {
        A.value[i] = value[i]*alpha;
    }
    return A;
}

SparseRow operator*(double alpha, const SparseRow& D)
{
    SparseRow A(D);
    for (int i = 0; i < A.nnz; i++)
    {
        A.value[i] = A.value[i]*alpha;
    }
    return A;
}

SparseRow SparseRow::operator/(double alpha)
{
#ifdef _DEBUG
    assert(alpha != 0);
#endif
    return (1.0/alpha)*(*this);
}

void SparseRow::operator=(const SparseRow& D)
{
    if( index != NULL) delete [] index;
    if( value != NULL) delete [] value;
    index = new long[D.nnz];
    value = new double[D.nnz];
    nnz = D.nnz;
    size = D.size;
    for (long i = 0; i < nnz; i++)
    {
        index[i] = D.index[i];
        value[i] = D.value[i];
    }
}

void SparseRow::operator*=(double alpha)
{
    for (int i = 0; i < nnz; i++)
    {
        value[i] = value[i]*alpha;
    }
}

void SparseRow::operator/=(double alpha)
{
    for (int i = 0; i < nnz; i++)
    {
        value[i] = value[i]*(1.0/alpha);
    }
}

void SparseRow::operator+=(const SparseRow& D)
{
#ifdef _DEBUG
    assert(size == D.size);
#endif
    if (nnz==0)
    {
        nnz = D.nnz;
        for (long i = 0; i < nnz; i++)
        {
            index[i] = D.index[i];
            value[i] = D.value[i];
        }
    }
    else
    {
        long iB = 0, iD = 0, i = 0;
        long* newIndex = new long[nnz+D.nnz];
        double* newValue = new double[nnz+D.nnz];
        
        while (iB<nnz || iD<D.nnz)
        {
            if (iB<nnz && iD<D.nnz)
            {
                if (index[iB]<D.index[iD])
                {
                    newIndex[i] = index[iB];
                    newValue[i] = value[iB];
                    iB++;
                }
                else if (index[iB]>D.index[iD])
                {
                    newIndex[i] = D.index[iD];
                    newValue[i] = D.value[iD];
                    iD++;
                }
                else
                {
                    newIndex[i] = D.index[iD];
                    newValue[i] = D.value[iD]+value[iB];
                    iD++;
                    iB++;
                }
            }
            else if (iB==nnz)
            {
                newIndex[i] = D.index[iD];
                newValue[i] = D.value[iD];
                iD++;
            }
            else if (iD==D.nnz)
            {
                newIndex[i] = index[iB];
                newValue[i] = value[iB];
                iB++;
            }
            i++;
        }
        nnz = i;
        delete [] index;
        delete [] value;
        index = newIndex;
        value = newValue;
    }
}

void SparseRow::operator-=(const SparseRow& D)
{
    (*this)+=((-1)*D);
}

double SparseRow::dot(const DoubleArray1D& B) const
{
    double sum = 0;
    for (int i = 0; i < nnz; i++)
    {
        sum += value[i]*B(index[i]);
    }
    return sum;
}

double SparseRow::norm(double p) const
{
    double sum = 0;
    for (int i = 0; i < nnz; i++)
    {
        sum += pow(value[i],p);
    }
    return pow(sum,1.0/p);
}

//
//  Output
//
ostream& operator<<(ostream& outStream, const SparseRow& V)
{
    
    long i;
    for(i = 0; i <  V.nnz; i++)
    {
        outStream << "(" << V.index[i] <<", "<<  setw(5) << V.value[i] << ")," << "     ";
    }
    outStream << endl;
    return outStream;
}

