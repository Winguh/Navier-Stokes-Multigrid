#ifndef __OUTPUT_H__
#define __OUTPUT_H__
#include "LinearAlgebra.h"

void outputVelocityx(int outnum,FILE* outer,DoubleArray1D& u);
void outputVelocityy(int outnum,FILE* outer,DoubleArray1D& u);
void outputPressure(int outnum,FILE* outer,DoubleArray1D& u);
void outputForcex(int outnum,FILE* outer,DoubleArray1D& u);
void outputForcey(int outnum,FILE* outer,DoubleArray1D& u);
// void outputVelocity(FILE* outer,DoubleArray1D& u,DoubleArray1D& v);

#endif