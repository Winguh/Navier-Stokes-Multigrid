#include <cmath>
#include <iostream>
#include "LinearAlgebra.h"

using namespace std;

void RestrictU(SparseMatrix& R, int n)
{
    int nr = n/2;
    int inr,ii;
    for (int i = 0;i<nr*(nr+1);i++)
    {
        inr = i/(nr+1);
        ii = i%(nr+1);
        if (i%(nr+1) == 0) //left
        {
            // R.addEntry(i,inr*2*(n+1)+ii*2,2.0);
            // R.addEntry(i,inr*2*(n+1)+ii*2+1,1.0);

            // R.addEntry(i,(inr*2+1)*(n+1)+ii*2,2.0);
            // R.addEntry(i,(inr*2+1)*(n+1)+ii*2+1,1.0);
        }
        else if (i%(nr+1) == nr) //right
        {
            // R.addEntry(i,inr*2*(n+1)+ii*2,2.0);
            // R.addEntry(i,inr*2*(n+1)+ii*2-1,1.0);

            // R.addEntry(i,(inr*2+1)*(n+1)+ii*2,2.0);
            // R.addEntry(i,(inr*2+1)*(n+1)+ii*2-1,1.0);
        }
        else
        {
            R.addEntry(i,inr*2*(n+1)+ii*2,2.0);
            R.addEntry(i,inr*2*(n+1)+ii*2+1,1.0);
            R.addEntry(i,inr*2*(n+1)+ii*2-1,1.0);

            R.addEntry(i,(inr*2+1)*(n+1)+ii*2,2.0);
            R.addEntry(i,(inr*2+1)*(n+1)+ii*2+1,1.0);
            R.addEntry(i,(inr*2+1)*(n+1)+ii*2-1,1.0);
        }
    }
    R = R/8.0;
}

void RestrictV(SparseMatrix& R, int n)
{
    int nr = n/2;
    int inr,ii;
    for (int i = 0;i<nr*(nr+1);i++)
    {
        inr = i/nr;
        ii = i%nr;
        if (i < nr) //bottom
        {
            // R.addEntry(i,inr*2*n+ii*2,2.0);
            // R.addEntry(i,inr*2*n+ii*2+1,2.0);

            // R.addEntry(i,inr*2*n+ii*2+n,1.0);
            // R.addEntry(i,inr*2*n+ii*2+1+n,1.0);
        }
        else if (i >= nr*nr) //top
        {
            // R.addEntry(i,inr*2*n+ii*2,2.0);
            // R.addEntry(i,inr*2*n+ii*2+1,2.0);

            // R.addEntry(i,inr*2*n+ii*2-n,1.0);
            // R.addEntry(i,inr*2*n+ii*2+1-n,1.0);
        }
        else
        {
            R.addEntry(i,inr*2*n+ii*2,2.0);
            R.addEntry(i,inr*2*n+ii*2+1,2.0);

            R.addEntry(i,inr*2*n+ii*2-n,1.0);
            R.addEntry(i,inr*2*n+ii*2+1-n,1.0);

            R.addEntry(i,inr*2*n+ii*2+n,1.0);
            R.addEntry(i,inr*2*n+ii*2+1+n,1.0);
        }
    }
    R = R/8.0;
}

void RestrictP(SparseMatrix& R, int n)
{
    int nr = n/2;
    int inr,ii;
    for (int i = 0;i<nr*nr;i++)
    {
        inr = i/nr;
        ii = i%nr;
        R.addEntry(i,inr*2*n+ii*2,1.0);
        R.addEntry(i,inr*2*n+ii*2+1,1.0);
        R.addEntry(i,inr*2*n+ii*2+n,1.0);
        R.addEntry(i,inr*2*n+ii*2+1+n,1.0);
    }
    R = R/4.0;
}

void ProlongateU(SparseMatrix& P, int n)
{
    int np = 2*n;
    int inr,inri,ii,iieo;
    for (int i = 0;i<np*(np+1);i++)
    {
        inr = i/(np+1);
        inri = inr/2;
        ii = i%(np+1);
        iieo = ii%2;
        if (ii == 0 || ii == np) //left or right
        {
            if (i<np)
            {
                P.addEntry(i,inri*(n+1)+ii/2,6.0);       
            }
            else if (i>np*np-2)
            {
                P.addEntry(i,inri*(n+1)+ii/2,6.0);                
            }
            else
            {
                if (inr%2 == 0)
                {
                    P.addEntry(i,inri*(n+1)+ii/2,6.0);
                    P.addEntry(i,(inri-1)*(n+1)+ii/2,2.0);         
                } 
                else
                {
                    P.addEntry(i,inri*(n+1)+ii/2,6.0);
                    P.addEntry(i,(inri+1)*(n+1)+ii/2,2.0);            
                }
            }

        }
        else
        {
            if (i<np)
            {
                if (iieo == 0)
                {
                    P.addEntry(i,inri*(n+1)+ii/2,6.0);
                }
                else
                {
                    P.addEntry(i,inri*(n+1)+ii/2,3.0);
                    P.addEntry(i,inri*(n+1)+ii/2+1,3.0);
                }                
            }
            else if (i>np*np-2)
            {
                if (iieo == 0)
                {
                    P.addEntry(i,inri*(n+1)+ii/2,6.0);
                }
                else
                {
                    P.addEntry(i,inri*(n+1)+ii/2,3.0);
                    P.addEntry(i,inri*(n+1)+ii/2+1,3.0);
                }         
            }
            else
            {
                if (inr%2 == 0)
                {
                    if (iieo == 0)
                    {
                        P.addEntry(i,inri*(n+1)+ii/2,6.0);
                        P.addEntry(i,(inri-1)*(n+1)+ii/2,2.0);
                    }
                    else
                    {
                        P.addEntry(i,inri*(n+1)+ii/2,3.0);
                        P.addEntry(i,inri*(n+1)+ii/2+1,3.0);
                        P.addEntry(i,(inri-1)*(n+1)+ii/2,1.0);
                        P.addEntry(i,(inri-1)*(n+1)+ii/2+1,1.0);
                    }                    
                } 
                else
                {
                    if (iieo == 0)
                    {
                        P.addEntry(i,inri*(n+1)+ii/2,6.0);
                        P.addEntry(i,(inri+1)*(n+1)+ii/2,2.0);
                    }
                    else
                    {
                        P.addEntry(i,inri*(n+1)+ii/2,3.0);
                        P.addEntry(i,inri*(n+1)+ii/2+1,3.0);
                        P.addEntry(i,(inri+1)*(n+1)+ii/2,1.0);
                        P.addEntry(i,(inri+1)*(n+1)+ii/2+1,1.0);
                    }                    
                }
            }
        }
    }
    P = P/8.0;
}

void ProlongateV(SparseMatrix& P, int n) 
{
    int np = 2*n;
    int inr,inri,ii;
    for (int i = 0;i<np*(np+1);i++)
    {
        inr = i/np;
        inri = inr/2;
        ii = i%np;
        if (i>(np-1) && i<np*np)
        {
            if (inr%2 == 0)
            {
                if (i%np == 0 || i%np == np-1)
                {
                    P.addEntry(i,inri*n+ii/2,6.0);
                }
                else if (ii%2 == 0)
                {
                    P.addEntry(i,inri*n+ii/2,6.0);
                    P.addEntry(i,inri*n+ii/2-1,2.0);
                }
                else
                {
                    P.addEntry(i,inri*n+ii/2,6.0);
                    P.addEntry(i,inri*n+ii/2+1,2.0);
                }
            }
            else
            {
                if (i%np == 0 || i%np == np-1)
                {
                    P.addEntry(i,inri*n+ii/2,3.0);
                    P.addEntry(i,inri*n+ii/2+n,3.0);
                }
                else if (ii%2 == 0)
                {
                    P.addEntry(i,inri*n+ii/2,3.0);
                    P.addEntry(i,inri*n+ii/2+n,3.0);
                    P.addEntry(i,inri*n+ii/2-1,1.0);
                    P.addEntry(i,inri*n+ii/2+n-1,1.0);
                }
                else
                {
                    P.addEntry(i,inri*n+ii/2,3.0);
                    P.addEntry(i,inri*n+ii/2+n,3.0);
                    P.addEntry(i,inri*n+ii/2+1,1.0);
                    P.addEntry(i,inri*n+ii/2+n+1,1.0);
                }
            }
        }

    }
    P = P/8.0;
}

void ProlongateP(SparseMatrix& P, int n)
{
    int np = 2*n;
    int inr,inri,ii;
    for (int i = 0;i<np*np;i++)
    {
        inr = i/np;
        inri = inr/2;
        ii = i%np;
        if (inr == 0)
        {
            P.addEntry(i,inri*n+ii/2,3.0);
            if (i%np != 0 && i%np != np-1)
            {
                if (ii%2 == 0)
                {
                    P.addEntry(i,inri*n+ii/2-1,1.0); 
                }
                else
                {
                    P.addEntry(i,inri*n+ii/2+1,1.0);
                }

            }
        }
        else if (inr == np-1)
        {
            P.addEntry(i,inri*n+ii/2,3.0);
            if (i%np != 0 && i%np != np-1)
            {
                if (ii%2 == 0)
                {
                    P.addEntry(i,inri*n+ii/2-1,1.0); 
                }
                else
                {
                    P.addEntry(i,inri*n+ii/2+1,1.0);
                }

            }

        }
        else
        {
            P.addEntry(i,inri*n+ii/2,3.0);
            if (inr%2 != 0)
            {
                P.addEntry(i,inri*n+ii/2+n,1.0);
                if (i%np != 0 && i%np != np-1)
                {
                    if (ii%2 == 0)
                    {
                        P.addEntry(i,inri*n+ii/2-1,1.0); 
                        P.addEntry(i,inri*n+ii/2-1+n,1.0); 
                    }
                    else
                    {
                        P.addEntry(i,inri*n+ii/2+1,1.0);
                        P.addEntry(i,inri*n+ii/2+n+1,1.0);
                    }
                }
            }
            else
            {
                P.addEntry(i,inri*n+ii/2-n,1.0);
                if (i%np != 0 && i%np != np-1)
                {
                    if (ii%2 == 0)
                    {
                        P.addEntry(i,inri*n+ii/2-1,1.0); 
                        P.addEntry(i,inri*n+ii/2-n-1,1.0); 
                    }
                    else
                    {
                        P.addEntry(i,inri*n+ii/2+1,1.0);
                        P.addEntry(i,inri*n+ii/2-n+1,1.0);
                    }
                }
                
            }


        }
    }
    P = P/6.0;
}
