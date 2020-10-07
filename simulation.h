#ifndef __simulation_H__
#define __simulation_H__

#include "LinearAlgebra.h"
#include "multigrid.h"

// void simulateTime(Multigrid& D_phi,Multigrid& D_star,DoubleArray1D& u,DoubleArray1D& v,DoubleArray1D& u_old,DoubleArray1D& v_old,DoubleArray1D& dp_dx,DoubleArray1D& dp_dy,SparseMatrix Dx,SparseMatrix Dy,SparseMatrix DD,SparseMatrix Dh);
void simulateTime(DoubleArray1D& u,DoubleArray1D& v,Multigrid D_star_u,Multigrid D_star_v,Multigrid D_phi,DoubleArray1D& u_old,DoubleArray1D& v_old,SparseMatrix UV,SparseMatrix VU,SparseMatrix DxU,SparseMatrix DyU,SparseMatrix DxV,SparseMatrix DyV,SparseMatrix DDU,SparseMatrix DDV,SparseMatrix DxPU,SparseMatrix DyPV,SparseMatrix DxUP,SparseMatrix DyVP,SparseMatrix Dhu, SparseMatrix Dhv);

void GSm(SparseMatrix A,DoubleArray1D f,DoubleArray1D& v,long iter);

#endif