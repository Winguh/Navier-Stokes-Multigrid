#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "LinearAlgebra.h"

#ifdef  _DEBUG 
#include <stdio.h>
#include <assert.h>
#endif

DoubleArray1D::DoubleArray1D()
{
    dataPtr       = 0;
    index1Size    = 0;
}

DoubleArray1D::DoubleArray1D(long m)
{
     dataPtr       = 0;
    index1Size    = 0;
    initialize(m);
};

DoubleArray1D::DoubleArray1D(double* d, long m)
{
    dataPtr       = 0;
    index1Size    = 0;
    initialize(d,m);
};

DoubleArray1D::DoubleArray1D(const DoubleArray1D& D)
{
    index1Size    = D.index1Size;
    dataPtr       = new double[index1Size];
    long i;
    for(i = 0; i < index1Size; i++){dataPtr[i] = D.dataPtr[i];}
};

DoubleArray1D::~DoubleArray1D()
{
    if( dataPtr != 0) delete [] dataPtr;
}

void DoubleArray1D::initialize(long m)
{

    if(index1Size != m)
    {
        delete [] dataPtr;
        dataPtr = new double[m];
    }
    index1Size    = m;

    long i;
    for(i = 0; i < index1Size; i++)
    {
        dataPtr[i] = 0.0;
    }
};

void DoubleArray1D::initialize(double* d, long m)
{

    initialize(m);

    long i;
    for(i = 0; i < index1Size; i++)
    {
        dataPtr[i] = d[i];
    }
};

//
//###################################################################
//                  Element Access 
//###################################################################
//

#ifdef _DEBUG 

double& DoubleArray1D::operator()(long i1)
{
    boundsCheck(i1, 0, index1Size-1,1);
    return *(dataPtr + i1);
};

const double& DoubleArray1D::operator()(long i1) const
{
    boundsCheck(i1, 0, index1Size-1,1);
    return *(dataPtr + i1);
};

#else

double& DoubleArray1D::operator()(long i1)
{
    return *(dataPtr + i1);
};

const double& DoubleArray1D::operator()(long i1) const
{
    return *(dataPtr + i1);
};

#endif

//
//###################################################################
//                Array Structure Access Functions
//###################################################################
//

//
// Resizes array to exactly newSize
//

void DoubleArray1D::set(long i, double v)
{
    dataPtr[i] = v;
}

void DoubleArray1D::resize(long newSize)
{
    long i;
    double*  newDataPtr = new double[newSize];
    double*  tmpDataPtr;

    if(newSize > index1Size) 
    {
        for(i = 0; i < index1Size; i++)
            newDataPtr[i] = dataPtr[i];
    }
    else
    {
        for(i = 0; i < newSize; i++)
            newDataPtr[i] = dataPtr[i];
    }

    index1Size = newSize;
    tmpDataPtr = dataPtr;
    dataPtr    = newDataPtr;

    if(tmpDataPtr != 0)
        delete [] tmpDataPtr;

}

//
//###################################################################
//                     Array Operators
//###################################################################
//
DoubleArray1D DoubleArray1D::operator+(const DoubleArray1D& D)
{
#ifdef _DEBUG
    sizeCheck(this->index1Size,D.index1Size);
#endif
    DoubleArray1D R(*this);
    long i;
    for(i = 0; i < index1Size; i++)
    {
        R.dataPtr[i] += D.dataPtr[i];
    }
    return R;
}

DoubleArray1D DoubleArray1D::operator-(const DoubleArray1D& D)
{
#ifdef _DEBUG
    sizeCheck(this->index1Size,D.index1Size);
#endif
    DoubleArray1D R(*this);
    long i;
    for(i = 0; i < index1Size; i++)
    {
        R.dataPtr[i] -= D.dataPtr[i];
    }
    return R;
}

DoubleArray1D DoubleArray1D::operator^(const DoubleArray1D& D)
{
    DoubleArray1D R(*this);
    long i;
    for (i = 0;i < index1Size;i++)
    {
        R.dataPtr[i] *= D.dataPtr[i];
    }
    return R;
}

DoubleArray1D DoubleArray1D::operator*(double alpha)
{
    DoubleArray1D R(*this);
    long i;
    for(i = 0; i < index1Size; i++)
    {
        R.dataPtr[i] *= alpha;
    }
    return R;
}

DoubleArray1D operator*(double alpha, const DoubleArray1D& D)
{
    DoubleArray1D R(D);
    long i;
    for(i = 0; i < D.index1Size; i++)
    {
        R.dataPtr[i] *= alpha;
    }
    return R;
}
                                                                ///<p>
DoubleArray1D DoubleArray1D::operator/(double alpha)
{
    DoubleArray1D R(*this);
    long i;
    for(i = 0; i < index1Size; i++)
    {
        R.dataPtr[i] /= alpha;
    }
    return R;
}
//
//###################################################################
//          
//###################################################################
//
void DoubleArray1D::operator=(const DoubleArray1D& D)
{
#ifdef _DEBUG
    if(index1Size != 0)
    {
        sizeCheck(this->index1Size,D.index1Size);
    }
#endif
    if(index1Size == 0)
    {
        initialize(D.index1Size);
        index1Size    = D.index1Size;
    }
//
//  copy over the data
//
    long i;
    for(i = 0; i < D.index1Size; i++)
    {
        dataPtr[i] = D.dataPtr[i];
    }
}



void DoubleArray1D::operator*=(double alpha)
{
    long i;
    for(i = 0; i < index1Size; i++)
    {
        dataPtr[i] *= alpha;
    }
}

void DoubleArray1D::operator/=(double alpha)
{
    long i;
    for(i = 0; i < index1Size; i++)
    {
        dataPtr[i] *= (1.0/alpha);
    }
}

void DoubleArray1D::operator+=(const DoubleArray1D& D)
{
#ifdef _DEBUG
    if(index1Size != 0)
    {
        sizeCheck(this->index1Size,D.index1Size);
    }
#endif

    long i;
    for(i = 0; i < D.index1Size; i++)
    {
        dataPtr[i] += D.dataPtr[i];
    }
}

void DoubleArray1D::operator-=(const DoubleArray1D& D)
{
#ifdef _DEBUG
    if(index1Size != 0)
    {
        sizeCheck(this->index1Size,D.index1Size);
    }
#endif

    long i;
    for(i = 0; i < D.index1Size; i++)
    {
        dataPtr[i] -= D.dataPtr[i];
    }
}

void DoubleArray1D::operator-=(double alpha)
{
    long i;
    for(i = 0; i < index1Size; i++)
    {
        dataPtr[i] -= alpha;
    }
}

void DoubleArray1D::setToValue(double d)
{
    long i;
    for(i = 0; i < index1Size; i++)
    {
        dataPtr[i] = d;
    }
}

double DoubleArray1D::dot(const DoubleArray1D& B) const
{
#ifdef _DEBUG
    sizeCheck(this->index1Size,B.index1Size);
#endif

    double R;
    R  = 0;
    long i;
    for(i = 0; i < index1Size; i++)
    {
        R += dataPtr[i]*B.dataPtr[i];
    }
    return R;
}

double DoubleArray1D::norm(double p) const
{
    double n=0;
    for (int i = 0; i<index1Size; i++)
        {
        n += pow(dataPtr[i],p);
        }
    n = pow(n, 1.0/p);
    return n;
}

double DoubleArray1D::mean() const
{
    double sm = 0.0;
    for (int i =0; i<index1Size;i++)
    {
        sm += dataPtr[i];
    }
    sm = sm/index1Size;
    return sm;
}

//
//  Output
//
ostream& operator<<(ostream& outStream, const DoubleArray1D& V)
{

    long i; 
    for(i = 0; i <  V.index1Size; i++)
    { 
        outStream <<  setw(5) << V.dataPtr[i] << " ";
        outStream << endl;
    }
    return outStream;
}

//
//###################################################################
//                      Bounds Checking
//###################################################################
//

#ifdef _DEBUG 
void DoubleArray1D::boundsCheck(long i, long begin, long end, int coordinate)
{
    if((i < begin)||(i  > end))
    {
        printf("Array index %d out of bounds \n",coordinate);
        printf("Offending index value %ld : Acceptable Range [%ld, %ld] \n",i, begin, end);
        assert(0);
    }
}
#else
void DoubleArray1D::boundsCheck(long, long, long, int){}
#endif

#ifdef _DEBUG 
void DoubleArray1D::sizeCheck(long size1, long size2)
{
    if(size1 != size2)
    {
        printf("Array Sizes Are Incompatable %ld != %ld \n", size1, size2);
        assert(0);
    }
}
#else
void DoubleArray1D::sizeCheck(long, long){}
#endif

