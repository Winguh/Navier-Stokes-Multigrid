#ifndef __LinearAlgebra_H__
#define __LinearAlgebra_H__

#include <iostream>
#include <iomanip>
using namespace std;

class DoubleArray1D 
{
public :
//
//###################################################################
//                 Constructors/Initialization
//###################################################################
//
    DoubleArray1D();
    DoubleArray1D(long m);
    DoubleArray1D(double* d, long m);
    DoubleArray1D(const DoubleArray1D& D);
    ~DoubleArray1D();
    void initialize(long m);
    void initialize(double* d, long m);

//
//###################################################################
//                  Element Access 
//###################################################################
//

#ifdef _DEBUG 
    double&  operator()(long i1);
    const double&  operator()(long i1) const;
#else
    double&  operator()(long i1);
    const double&  operator()(long i1) const;
#endif

//
//###################################################################
//                Array Structure Access Functions
//###################################################################
//

    double* getDataPointer(){return dataPtr;};
    long getIndex1Size()  const {return index1Size;}
    void set(long i, double v);

//
//  Get/Set specifically for one dimensional arrays
//
    long getSize()       const {return index1Size;}

    //
    // Resizes array to exactly newSize
    //
    void resize(long newSize);

//
//###################################################################
//                     Array Operators
//###################################################################
//
    DoubleArray1D operator+(const DoubleArray1D& D);
    DoubleArray1D operator-(const DoubleArray1D& D);
    DoubleArray1D operator*(double alpha);
    friend DoubleArray1D operator*(double alpha, const DoubleArray1D& D);
    DoubleArray1D operator/(double alpha);
    DoubleArray1D operator^(const DoubleArray1D& D);
    
//
//###################################################################
//          
//###################################################################
//
    void operator=(const DoubleArray1D& D);
    void operator*=(double alpha);
    void operator/=(double alpha);
    void operator+=(const DoubleArray1D& D);
    void operator-=(const DoubleArray1D& D);
    void operator-=(double alpha);
    void setToValue(double d);
    double dot(const DoubleArray1D& B) const;
    double norm(double p) const;
    double mean() const;

//
//  Output
//
    friend ostream& operator<<(ostream& outStream, const DoubleArray1D& V);
    
//
//###################################################################
//                      Class Data Members
//###################################################################
//
    protected :

    double*      dataPtr;     // data pointer
    long    index1Size;       // coordinate 1 size

//
//###################################################################
//                      Bounds Checking
//###################################################################
//

    static void boundsCheck(long, long, long, int);
    static void sizeCheck(long, long);

    
};

// class DoubleArray2D
// {
    
//     public :
//     //
//     //###################################################################
//     //                 Constructors/Initialization
//     //###################################################################
//     //
//     DoubleArray2D();
//     DoubleArray2D(long m, long n);
//     DoubleArray2D(double* d, long m, long n);
//     DoubleArray2D(const DoubleArray2D& D);
//     ~DoubleArray2D();
//     void initialize(long m, long n);
//     void initialize(double* d, long m, long n);
//     //
//     //###################################################################
//     //                  Element Access
//     //###################################################################
//     //
//     double&  operator()(long i1, long i2);
//     const double&  operator()(long i1, long i2) const;
//     DoubleArray1D row(long i);
//     void setrow(long i, const DoubleArray1D& x);
//     void set(long i, long j, double m);
//     DoubleArray1D  column(long j);
//     void setcolumn(long j, const DoubleArray1D& x);
    
    
//     //
//     //###################################################################
//     //                Array Structure Access Functions
//     //###################################################################
//     //
//     double* getDataPointer(){return dataPtr;};
//     long getIndex1Size()  const {return index1Size;}
//     long getIndex2Size()  const {return index2Size;}
    
//     //
//     //###################################################################
//     //                     Array Operators
//     //###################################################################
//     //
    
//     DoubleArray2D operator+(const DoubleArray2D& D);
//     DoubleArray2D operator-(const DoubleArray2D& D);
//     DoubleArray2D operator*(double alpha);
//     DoubleArray1D operator*(const DoubleArray1D& x);
//     DoubleArray2D operator*(const DoubleArray2D& D);
//     friend DoubleArray2D operator*(double alpha, const DoubleArray2D& D);
//     DoubleArray2D operator/(double alpha);
//     void operator=(const DoubleArray2D& D);
//     void operator*=(double alpha);
//     void operator+=(const DoubleArray2D& D);
//     void operator-=(const DoubleArray2D& D);
//     void setToValue(double d);
//     double dot(const DoubleArray2D& D) const;
    
//     //  Input/Output
//     friend ostream&  operator <<(ostream& outStream, const DoubleArray2D& A);
    
//     //
//     //###################################################################
//     //                      Class Data Members
//     //###################################################################
//     //
//     protected :
    
//     double*      dataPtr;     // data pointer
//     long      index1Size;     // coordinate 1 size
//     long      index2Size;     // coordinate 2 size
    
//     //
//     //###################################################################
//     //                      Bounds Checking
//     //###################################################################
//     //
    
//     static void boundsCheck(long, long, long, int);
    
//     static void sizeCheck(long, long, long, long);
    
// };

class SparseRow
{
    public:
    //
    //###################################################################
    //                 Constructors/Initialization
    //###################################################################
    //
    SparseRow();
    SparseRow(long m);
    SparseRow(const SparseRow& D);
    ~SparseRow();
    void initialize(long m);
    
    //
    //###################################################################
    //                  Element Access
    //###################################################################
    //
    
    double  operator()(long i1);
    
    //
    //###################################################################
    //                Array Structure Access Functions
    //###################################################################
    //
    
    double* getValue(){return value;};
    long* getIndex(){return index;};
    long getNnz() const {return nnz;};
    long getIndex1Size()  const {return size;};
    void addEntry(long i, double v);
    
    //
    //###################################################################
    //                     Array Operators
    //###################################################################
    //
    SparseRow operator+(const SparseRow& D);
    SparseRow operator-(const SparseRow& D);
    SparseRow operator*(double alpha);
    friend SparseRow operator*(double alpha, const SparseRow& D);
    SparseRow operator/(double alpha);
    
    void operator=(const SparseRow& D);
    void operator*=(double alpha);
    void operator/=(double alpha);
    void operator+=(const SparseRow& D);
    void operator-=(const SparseRow& D);
    double dot(const DoubleArray1D& B) const;
    double norm(double p) const;
    
    //
    //  Output
    //
    friend ostream& operator<<(ostream& outStream, const SparseRow& V);
    
    //
    //###################################################################
    //                      Class Data Members
    //###################################################################
    //
protected:
    long size;
    long nnz;
    long* index;
    double* value;
};

class SparseMatrix
{
    public :
    //
    //###################################################################
    //                 Constructors/Initialization
    //###################################################################
    //
    SparseMatrix();
    SparseMatrix(long m, long n);
    SparseMatrix(const SparseMatrix& D);
    ~SparseMatrix();
    void initialize(long m, long n);
    
    //
    //###################################################################
    //                  Element Access
    //###################################################################
    //
    
    double  operator()(long i1, long i2);
    SparseRow&  operator()(long i1);
    const SparseRow&  operator()(long i1) const;
    
    //
    //###################################################################
    //                Array Structure Access Functions
    //###################################################################
    //
    
    const SparseRow& getRow(long i){return rows[i];};
    long getSize1() const {return size1;};
    long getSize2() const {return size2;};
    void addEntry(long i, long j, double v);
    
    //
    //###################################################################
    //                     Array Operators
    //###################################################################
    //
    SparseMatrix operator+(const SparseMatrix& D);
    SparseMatrix operator-(const SparseMatrix& D);
    SparseMatrix operator*(double alpha);
    friend SparseMatrix operator*(double alpha, const SparseMatrix& D);
    SparseMatrix operator/(double alpha);
    SparseMatrix operator*(const SparseMatrix D) const;
    // SparseMatrix operator*(SparseMatrix& D);
    
    void operator=(const SparseMatrix& D);
    void operator*=(double alpha);
    void operator+=(const SparseMatrix& D);
    void operator-=(const SparseMatrix& D);
    DoubleArray1D operator*(const DoubleArray1D& B) const;
    DoubleArray1D transMultiply(const DoubleArray1D& B) const;
    SparseMatrix transpose() const;
    void resetRow(long i);
    
    //
    //  Output
    //
    friend ostream& operator<<(ostream& outStream, const SparseMatrix& V);
    
    //
    //###################################################################
    //                      Class Data Members
    //###################################################################
    //
protected:
    long size1;
    long size2;
    SparseRow* rows;
};


#endif


 
