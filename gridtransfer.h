#ifndef __gridtransfer_H__
#define __gridtransfer_H__
#include <cmath>
#include <omp.h>
#include "LinearAlgebra.h"

void RestrictU(SparseMatrix& R, int n);
void RestrictV(SparseMatrix& R, int n);
void RestrictP(SparseMatrix& R, int n);

void ProlongateU(SparseMatrix& P, int n);
void ProlongateV(SparseMatrix& P, int n);
void ProlongateP(SparseMatrix& P, int n);

#endif