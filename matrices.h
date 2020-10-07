#ifndef __Matrices_H__
#define __Matrices_H__
#include <cmath>
#include <omp.h>
#include "LinearAlgebra.h"

void UtoV(SparseMatrix& A,int n);
void VtoU(SparseMatrix& A,int n);

void UtoP(SparseMatrix& A,int n);
void VtoP(SparseMatrix& A,int n);

void GradxU(SparseMatrix &A,int n);
void GradyU(SparseMatrix &A,int n);

void GradxV(SparseMatrix& A,int n);
void GradyV(SparseMatrix& A,int n);

void GradxPtoU(SparseMatrix& A,int n);
void GradyPtoV(SparseMatrix& A,int n);

void GradxUtoP(SparseMatrix& A, int n);
void GradyVtoP(SparseMatrix& A, int n);

void LaplaceU(SparseMatrix& A, int n);
void LaplaceV(SparseMatrix& A, int n);

void GradxC(SparseMatrix& A, int n);
void GradyC(SparseMatrix& A, int n);
void LaplaceC(SparseMatrix& A, int n);
void createCm(SparseMatrix& A, int n);

void createDstaru(SparseMatrix& A, int n);
void createDstarv(SparseMatrix& A, int n);
void createDphi(SparseMatrix& A, int n);

double f1(double r);
double DistnF(double r);
#endif