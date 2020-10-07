#include <iostream>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>
#include <direct.h>
#include <cstring>
#include "LinearAlgebra.h"
// #include "Grid.h"
#include "parameters.h"

FILE* outer;

void outputVelocityx(int outnum,FILE* outer,DoubleArray1D& u)
{
    char outname[500];
    strcpy(outname,DirName);
    strcat(outname,"\\u");
    // mkdir(outname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(outname);
    sprintf(outname,"%s/%d",outname,outnum);
	strcat(outname,".txt");
    outer = fopen(outname, "w");
    rewind(outer);
    for (int i = 0;i<u.getIndex1Size();i++)
    {
        fprintf(outer,"%6E\n",u(i));
    }
    fflush(outer);
    fclose(outer);
}

void outputVelocityy(int outnum,FILE* outer,DoubleArray1D& u)
{
    char outname[500];
    strcpy(outname,DirName);
    strcat(outname,"\\v");
    // mkdir(outname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(outname);
    sprintf(outname,"%s/%d",outname,outnum);
	strcat(outname,".txt");
    outer = fopen(outname, "w");
    rewind(outer);
    for (int i = 0;i<u.getIndex1Size();i++)
    {
        fprintf(outer,"%6E\n",u(i));
    }
    fflush(outer);
    fclose(outer);
}

void outputPressure(int outnum,FILE* outer,DoubleArray1D& u)
{
    char outname[500];
    strcpy(outname,DirName);
    strcat(outname,"\\p");
    // mkdir(outname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(outname);
    sprintf(outname,"%s/%d",outname,outnum);
	strcat(outname,".txt");
    outer = fopen(outname, "w");
    rewind(outer);
    for (int i = 0;i<u.getIndex1Size();i++)
    {
        fprintf(outer,"%6E\n",u(i));
    }
    fflush(outer);
    fclose(outer);
}

void outputForcex(int outnum,FILE* outer,DoubleArray1D& u)
{
    char outname[500];
    strcpy(outname,DirName);
    strcat(outname,"\\fx");
    // mkdir(outname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(outname);
    sprintf(outname,"%s/%d",outname,outnum);
	strcat(outname,".txt");
    outer = fopen(outname, "w");
    rewind(outer);
    for (int i = 0;i<u.getIndex1Size();i++)
    {
        fprintf(outer,"%6E\n",u(i));
    }
    fflush(outer);
    fclose(outer);
}

void outputForcey(int outnum,FILE* outer,DoubleArray1D& u)
{
    char outname[500];
    strcpy(outname,DirName);
    strcat(outname,"\\fy");
    // mkdir(outname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(outname);
    sprintf(outname,"%s/%d",outname,outnum);
	strcat(outname,".txt");
    outer = fopen(outname, "w");
    rewind(outer);
    for (int i = 0;i<u.getIndex1Size();i++)
    {
        fprintf(outer,"%6E\n",u(i));
    }
    fflush(outer);
    fclose(outer);
}