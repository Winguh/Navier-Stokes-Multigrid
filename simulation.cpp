#include <cmath>
#include <chrono>
#include "simulation.h"
#include "LinearAlgebra.h"
#include "matrices.h"
#include "parameters.h"
#include "output.h"
#include "multigrid.h"

void simulateTime(DoubleArray1D& u,DoubleArray1D& v,Multigrid D_star_u,Multigrid D_star_v,Multigrid D_phi,DoubleArray1D& u_old,DoubleArray1D& v_old,SparseMatrix UV,SparseMatrix VU,SparseMatrix DxU,SparseMatrix DyU,SparseMatrix DxV,SparseMatrix DyV,SparseMatrix DDU,SparseMatrix DDV,SparseMatrix DxPU,SparseMatrix DyPV,SparseMatrix DxUP,SparseMatrix DyVP,SparseMatrix Dhu, SparseMatrix Dhv)
{    
    int sizet = Tfinal/dt + 2;
    int sizeM = n*(n+1);
    int sizeP = n*n;

    DoubleArray1D T(sizet);
    for (int i = 0;i<sizet;i++)
    {
        T(i) = i*dt;
    }

    SparseMatrix UP(sizeP,sizeP);
    SparseMatrix VP(sizeP,sizeP);
    DoubleArray1D ut(sizeP);
    DoubleArray1D vt(sizeP);
    DoubleArray1D fxt(alpha);
    DoubleArray1D fyt(alpha);

    UtoP(UP,n);
    VtoP(VP,n);

    //initialize vectors for calculations
    DoubleArray1D du_dx(sizeM);
    DoubleArray1D du_dy(sizeM);
    DoubleArray1D dv_dx(sizeM);
    DoubleArray1D dv_dy(sizeM);

    DoubleArray1D du_old_dx(sizeM);
    DoubleArray1D du_old_dy(sizeM);
    DoubleArray1D dv_old_dx(sizeM);
    DoubleArray1D dv_old_dy(sizeM);

    DoubleArray1D uv(sizeM);
    DoubleArray1D vu(sizeM);
    DoubleArray1D uv_old(sizeM);
    DoubleArray1D vu_old(sizeM);

    DoubleArray1D dp_dx(sizeM);
    DoubleArray1D dp_dy(sizeM);

    DoubleArray1D ddu(sizeM);
    DoubleArray1D ddv(sizeM);
    DoubleArray1D u_star(sizeM);
    DoubleArray1D v_star(sizeM);
    DoubleArray1D u_rhs(sizeM);
    DoubleArray1D v_rhs(sizeM);

    DoubleArray1D du_star_dx(sizeP);
    DoubleArray1D dv_star_dy(sizeP);
    DoubleArray1D phi_rhs(sizeP);
    DoubleArray1D phi(sizeP);

    DoubleArray1D dphi_dx(sizeM);
    DoubleArray1D dphi_dy(sizeM);
    DoubleArray1D u_new(sizeM);
    DoubleArray1D v_new(sizeM);
    DoubleArray1D d_phi_rhs_dx(sizeM);
    DoubleArray1D d_phi_rhs_dy(sizeM);

    DoubleArray1D pres(sizeP);

    DoubleArray1D phiR(sizeP);
    DoubleArray1D phiup(sizeP);
    DoubleArray1D uR(sizeM);
    DoubleArray1D uup(sizeM);
    DoubleArray1D vR(sizeM);
    DoubleArray1D vup(sizeM);    

    DoubleArray1D dtU(nx*ny);
    DoubleArray1D dtV(nx*ny);
    DoubleArray1D fKu(alpha);
    DoubleArray1D fKv(alpha);

    DoubleArray1D FKu(nx*(nx+1));
    DoubleArray1D FKv(nx*(nx+1));

    DoubleArray1D ubdry(sizeM);
    DoubleArray1D vbdry(sizeM);
    DoubleArray1D ubdrynew(sizeM);
    DoubleArray1D vbdrynew(sizeM);
    DoubleArray1D ubdryold(sizeM);
    DoubleArray1D vbdryold(sizeM);
    
    DoubleArray1D usoln(sizeM);
    DoubleArray1D vsoln(sizeM);
    DoubleArray1D psoln(sizeP);
    DoubleArray1D pxsoln(sizeP);
    DoubleArray1D pysoln(sizeP);

    DoubleArray1D gradtestx(sizeM);
    DoubleArray1D gradtesty(sizeM);
    
    if (TG == 1)
    {
        //apply initial condition for velocity, pressure gradient
        for (int i = 0;i<n+1;i++)
        {
            for (int j = 0;j<n;j++)
            {
                u(i+j*(n+1)) = cos(i*h)*sin((0.5+j)*h);
                dp_dx(i+j*(n+1)) = 0.5*sin(2.0*i*h);
            }
        }
        for (int i = 0;i<n;i++)
        {
            for (int j = 0;j<n+1;j++)
            {
                v(i+j*n) = -sin((i+0.5)*h)*cos(j*h);
                dp_dy(i+j*n) = 0.5*sin(2.0*j*h);
            }
        }
        for (int i  = 0;i<n;i++)
        {
            for (int j = 0;j<n;j++)
            {
            pres(i+j*n) = -0.25*(cos(2.0*(i+0.5)*h) + cos(2.0*(j+0.5)*h));
            }
        }
    }


    if (LDC == 1)
    {
        u.setToValue(0.0);
        v.setToValue(0.0);
        u_old = u;
        v_old = v;
        //init u boundary
        for(int i = 1;i<n;i++) //interior bottom
        {
            ubdry(i) = 0.0;
            ubdrynew(i) = 0.0;
            ubdryold(i) = 0.0;
        }
        for(int i = n*n;i<n*(n+1)-1;i++) //interior top
        {
            ubdry(i) = 1.0;
            ubdrynew(i) = 1.0;
            ubdryold(i) = 1.0;
        }
        for (int i=0;i<n*n;i+=n+1)//left
        {
            ubdry(i) = 0.0;
            ubdrynew(i) = 0.0;
            ubdryold(i) = 0.0;
        }
        for (int i=n;i<n*(n+1);i+=n+1)//right
        {
            ubdry(i) = 0.0;
            ubdrynew(i) = 0.0;
            ubdryold(i) = 0.0;
        }

        //init v boundary
        for (int i=0;i<n;i++)//bottom
        {
            vbdry(i) = 0.0;
            vbdrynew(i) = 0.0;
            vbdryold(i) = 0.0;
        }
        for (int i=n*n;i<n*(n+1);i++)//top
        {
            vbdry(i) = 0.0;
            vbdrynew(i) = 0.0;
            vbdryold(i) = 0.0;
        }
        for (int i = n;i<n*n;i+=n) //interior left
        {
            vbdry(i) = 0.0;
            vbdrynew(i) = 0.0;
            vbdryold(i) = 0.0;
        }
        for (int i = 2*n-1;i<n*(n+1)-1;i+=n) //interior right
        {
            vbdry(i) = 0.0;
            vbdrynew(i) = 0.0;
            vbdryold(i) = 0.0;
        }
    }

    DoubleArray1D uerr(sizeM);
    DoubleArray1D verr(sizeM);
    DoubleArray1D perr(sizeM);
    DoubleArray1D pxerr(sizeM);
    DoubleArray1D pyerr(sizeM);

    double respx = 0.0;
    double respy = 0.0;
    double relu = 0.0;
    double relv = 0.0;
    double relp = 0.0;

    long outnum = 0;
    FILE* outer;
    double resu = 0.0,resv = 0.0,resp = 0.0;
    double len1,len2,len3;

    DoubleArray1D dvgnc(sizeP);

    int its;
    for (int k=0;k<sizet;k++)
    {
        if (TG == 1)
        {
            // update analytic solution for error analysis
            for (int i  = 0;i<n;i++)
            {
                for (int j = 0;j<n;j++)
                {
                    if (k == 0)
                    {
                        psoln(i+j*n) = -0.25*(cos(2.0*(i+0.5)*h) + cos(2.0*(j+0.5)*h));
                    }
                    else
                    {
                        psoln(i+j*n) = -0.25*(cos(2.0*(i+0.5)*h) + cos(2.0*(j+0.5)*h))*exp(-4.0*nu*(dt*(k-0.5)));

                    }
                }
            }
            for (int i = 0;i<n+1;i++)
            {
                for (int j = 0;j<n;j++)
                {
                    
                    usoln(i+j*(n+1)) = cos(i*h)*sin((0.5+j)*h)*exp(-2.0*nu*T(k));
                    // if (k == 0)
                    // {
                    //     pxsoln(i+j*(n+1)) = 0.0;
                    // }
                    // else
                    // {
                    //     pxsoln(i+j*(n+1)) = 0.5*sin(2.0*i*h)*exp(-4.0*nu*(dt*(k-0.5)));
                    // }
                }
            }
            for (int i = 0;i<n;i++)
            {
                for (int j = 0;j<n+1;j++)
                {
                    vsoln(i+j*n) = -sin((i+0.5)*h)*cos(j*h)*exp(-2.0*nu*T(k));
                    // if (k == 0)
                    // {
                    //     pysoln(i+j*n) = 0.0;
                    // }
                    // else
                    // {
                    //     pysoln(i+j*n) = 0.5*sin(2.0*j*h)*exp(-4.0*nu*(dt*(k-0.5)));
                    // }
                }
            }
        
            uerr = u-usoln;
            verr = v-vsoln;
            perr = pres-psoln;
            // pxerr = dp_dx-pxsoln;
            // pyerr = dp_dy-pysoln;

            //update u boundary
            for(int i = 1;i<n;i++) //interior bottom
            {
                ubdry(i) = 0.0;
                ubdrynew(i) = 0.0;
                ubdryold(i) = 0.0;
            }
            for(int i = n*n;i<n*(n+1)-1;i++) //interior top
            {
                ubdry(i) = 0.0;
                ubdrynew(i) = 0.0;
                ubdryold(i) = 0.0;
            }
            for (int i=0;i<n*n;i+=n+1)//left
            {
                ubdrynew(i) = sin(((i/(n+1))+0.5)*h)*exp(-2.0*nu*T(k+1));
            }
            for (int i=n;i<n*(n+1);i+=n+1)//right
            {
                ubdrynew(i) = sin(((i/(n+1))+0.5)*h)*exp(-2.0*nu*T(k+1));
            }

            //update v boundary
            for (int i = n;i<n*n;i+=n) //interior left
            {
                vbdry(i) = 0.0;
                vbdrynew(i) = 0.0;
                vbdryold(i) = 0.0;
            }
            for (int i = 2*n-1;i<n*(n+1)-1;i+=n) //interior right
            {
                vbdry(i) = 0.0;
                vbdrynew(i) = 0.0;
                vbdryold(i) = 0.0;
            }
            for (int i=0;i<n;i++)//bottom
            {
                vbdrynew(i) = -sin((i+0.5)*h)*exp(-2.0*nu*T(k+1));
            }
            for (int i=n*n;i<n*(n+1);i++)//top
            {
                vbdrynew(i) = -sin(((i%n)+0.5)*h)*exp(-2.0*nu*T(k+1));
            }
        }
        
        //start section 2
        auto start = chrono::high_resolution_clock::now();

        //center difference first derviative of x velocity
        du_dx = DxU*u;
        du_dy = DyU*u;

        //center difference first derivative of previous x velocity
        du_old_dx = DxU*u_old;
        du_old_dy = DyU*u_old;

        //laplacian of x velocity
        ddu = DDU*u;

        vu = VU*v;//v on u grid
        vu_old = VU*v_old;

        //add part from Dirichlet interpolation to inner bottom/top
        for(int i = 1;i<n;i++) //interior bottom
        {
            du_dy(i) -= (8.0/3.0)*ubdry(i)/(2.0*h);
            du_old_dy(i) -= (8.0/3.0)*ubdryold(i)/(2.0*h);
            ddu(i) += (8.0/3.0)*(ubdry(i) + ubdrynew(i))/pow(h,2.0);
        }
        for(int i = n*n;i<n*(n+1)-1;i++) //interior top
        {
            du_dy(i) += (8.0/3.0)*ubdry(i)/(2.0*h);
            du_old_dy(i) += (8.0/3.0)*ubdryold(i)/(2.0*h);
            ddu(i) += (8.0/3.0)*(ubdry(i) + ubdrynew(i))/pow(h,2.0);
        }

        //center difference first derviative of y velocity
        dv_dx = DxV*v;
        dv_dy = DyV*v;

        //center difference first derivative of previous y velocity
        dv_old_dx = DxV*v_old;
        dv_old_dy = DyV*v_old;

        //laplacian of y velocity
        ddv = DDV*v;

        uv = UV*u;//u on v grid
        uv_old = UV*u_old;

        for (int i = n;i<n*n;i+=n) //interior left
        {
            dv_dx(i) -= (8.0/3.0)*vbdry(i)/(2.0*h);
            dv_old_dx(i) -= (8.0/3.0)*vbdryold(i)/(2.0*h);
            ddv(i) += (8.0/3.0)*(vbdry(i) + vbdrynew(i))/pow(h,2.0);
        }
        for (int i = 2*n-1;i<n*(n+1)-1;i+=n) //interior right
        {
            dv_dx(i) += (8.0/3.0)*vbdry(i)/(2.0*h);
            dv_old_dx(i) += (8.0/3.0)*vbdryold(i)/(2.0*h);
            ddv(i) += (8.0/3.0)*(vbdry(i) + vbdrynew(i))/pow(h,2.0);
        }

        if (IBM == 1)
        {
            dtU = ((-1.0)*u)/dt;
            dtV = ((-1.0)*v)/dt;

            fKu = Dhu*(dtU + (1.5)*((u^du_dx) + (vu^du_dy)) - (0.5)*((u_old^du_old_dx) + (vu_old^du_old_dy)) - (nu/2.0)*ddu + dp_dx)*pow(h,2);
            fKv = Dhv*(dtV + (1.5)*((uv^dv_dx) + (v^dv_dy)) - (0.5)*((uv_old^dv_old_dx) + (v_old^dv_old_dy)) - (nu/2.0)*ddv + dp_dy)*pow(h,2);

            FKu = Dhu.transMultiply(fKu)*pow(s,2);
            FKv = Dhv.transMultiply(fKv)*pow(s,2);
        }
        else
        {
            FKu.setToValue(0.0);
            FKv.setToValue(0.0);
        }

        //forms right hand side for systems to solve u_star, v_star
        u_rhs = -(1.5)*((u^du_dx) + (vu^du_dy)) + (0.5)*((u_old^du_old_dx) + (vu_old^du_old_dy)) - dp_dx + (nu/2.0)*ddu + FKu;
        u_rhs = u + u_rhs*dt;
        u_rhs = -2.0/(dt*nu)*u_rhs;

        v_rhs = -(1.5)*((uv^dv_dx) + (v^dv_dy)) + (0.5)*((uv_old^dv_old_dx) + (v_old^dv_old_dy)) - dp_dy + (nu/2.0)*ddv + FKv;
        v_rhs = v + v_rhs*dt;
        v_rhs = -2.0/(dt*nu)*v_rhs;

        //apply boundary conditions
        for (int i=0;i<n*n;i+=n+1)//left
        {
            u_rhs(i) = ubdrynew(i)/pow(h,2.0);
        }
        for (int i=n;i<n*(n+1);i+=n+1)//right
        {
            u_rhs(i) = ubdrynew(i)/pow(h,2.0);
        }
        for (int i=0;i<n;i++)//bottom
        {
            v_rhs(i) = vbdrynew(i)/pow(h,2.0);
        }
        for (int i=n*n;i<n*(n+1);i++)//top
        {
            v_rhs(i) = vbdrynew(i)/pow(h,2.0);
        }

        // cout << "u\n";
        start = chrono::high_resolution_clock::now();
        D_star_u.solve(u_rhs,u_star);
        auto stop = chrono::high_resolution_clock::now();
        auto dur = chrono::duration_cast<chrono::microseconds>(stop-start);
        uR = u_rhs - D_star_u.Grid.M[0]*u_star; 
        resu = uR.norm(2);
        // cout << dur.count() << " " << resu << "\n";
        its = 0;
        while (resu > tol && its < 10)
        {
            start = chrono::high_resolution_clock::now();
            D_star_u.solve(uR,uup);
            stop = chrono::high_resolution_clock::now();
            dur = chrono::duration_cast<chrono::microseconds>(stop-start);

            u_star = u_star + uup;
            uR = u_rhs - D_star_u.Grid.M[0]*u_star;    
            resu = uR.norm(2);
            // cout << dur.count() << " " << resu << "\n";
            its++;
        }
        
        // cout << "v\n";
        start = chrono::high_resolution_clock::now();
        D_star_v.solve(v_rhs,v_star);
        stop = chrono::high_resolution_clock::now();
        dur = chrono::duration_cast<chrono::microseconds>(stop-start);
        vR = v_rhs - D_star_v.Grid.M[0]*v_star;
        resv = vR.norm(2);
        // cout << dur.count() << " " << resv << "\n";
        its = 0;
        while (resv > tol && its < 10)
        {
            start = chrono::high_resolution_clock::now();
            D_star_v.solve(vR,vup);
            stop = chrono::high_resolution_clock::now();
            dur = chrono::duration_cast<chrono::microseconds>(stop-start);
            v_star = v_star + vup;
            vR = v_rhs - D_star_v.Grid.M[0]*v_star;
            resv = vR.norm(2);
            // cout << dur.count() << " " << resv << "\n";
            its++;
        }  

        // end section 2
        stop = chrono::high_resolution_clock::now();
        dur = chrono::duration_cast<chrono::milliseconds>(stop-start);
        len1 = dur.count();

        // start section 3
        start = chrono::high_resolution_clock::now();

        du_star_dx = DxUP*u_star;
        dv_star_dy = DyVP*v_star;

        //forms rhs for systems to solve phi
        phi_rhs = du_star_dx + dv_star_dy;
        phi_rhs /= dt;
        
        phi_rhs -= phi_rhs.mean();

        // cout << "p\n";
        start = chrono::high_resolution_clock::now();
        D_phi.solve(phi_rhs,phi);
        stop = chrono::high_resolution_clock::now();
        dur = chrono::duration_cast<chrono::microseconds>(stop-start);
        phiR = phi_rhs - D_phi.Grid.M[0]*phi;
        resp = phiR.norm(2);
        // cout << dur.count() << " " << resp << "\n";
        its = 0;
        while (resp > tol && its < 10)
        {
            start = chrono::high_resolution_clock::now();
            D_phi.solve(phiR,phiup);
            stop = chrono::high_resolution_clock::now();
            dur = chrono::duration_cast<chrono::microseconds>(stop-start);
            phi = phi + phiup;
            phiR = phi_rhs - D_phi.Grid.M[0]*phi;
            resp = phiR.norm(2);
            // cout << dur.count() << " " << resp << "\n";
            its++;
        }
        
        //Gram-Schmidt step to enforce convergence criterion
        phi -= phi.mean();

        //end section 3
        stop = chrono::high_resolution_clock::now();
        auto dur2 = chrono::duration_cast<chrono::milliseconds>(stop-start);
        len2 = dur2.count();

        //start section 4
        start = chrono::high_resolution_clock::now();

        dphi_dx = DxPU*phi;
        dphi_dy = DyPV*phi;

        pres += phi - phi_rhs*(nu*dt/2.0);

        u_new = u_star - dphi_dx*dt;
        v_new = v_star - dphi_dy*dt;

        d_phi_rhs_dx = DxPU*phi_rhs;
        d_phi_rhs_dy = DyPV*phi_rhs;

        dp_dx += dphi_dx - (d_phi_rhs_dx)*(nu*dt/2.0);
        dp_dy += dphi_dy - (d_phi_rhs_dy)*(nu*dt/2.0);

        u_old = u;
        u = u_new;

        v_old = v;
        v = v_new;

        //end section 4
        stop = chrono::high_resolution_clock::now();
        auto dur3 = chrono::duration_cast<chrono::microseconds>(stop-start);
        len3 = dur3.count();

        // relu = uerr.norm(2.0)/usoln.norm(2.0);
        // relv = verr.norm(2.0)/vsoln.norm(2.0);

        // resu = uerr.norm(2.0)/n;
        // resv = verr.norm(2.0)/n;
        // resp = perr.norm(2.0)/n;
        // respx = pxerr.norm(2.0)/n;
        // respy = pyerr.norm(2.0)/n;

        // relp = perr.norm(2.0)/psoln.norm(2.0);

        //output data to txt files
        if (k%1 == 0)
        {   
            ut = UP*u;
            vt = VP*v;
            fxt = Dhu*u*pow(h,2);
            fyt = Dhv*v*pow(h,2);

            printf("%ld\t: %f | %f | %f | %f || Time: %.0f | %.0f | %.0f | %.0f ms  ||  %f\n",outnum,resu,resv,resp,(DxUP*u+DyVP*v).norm(2),len1,len2,len3,(len1+len2+len3*0.001),dt*(k+1));
            // printf("%ld\t %6E  %6E  %6E  %6E  %6E  %6E  %6E  %6E  %f\n",outnum,resu,relu,resv,relv,resp,relp,respx,respy,dt*k);
            // printf("%6E  %6E\n",resu,resv);
            outputVelocityx(outnum,outer,ut);
            outputVelocityy(outnum,outer,vt);
            outputPressure(outnum,outer,pres);
            // outputForcex(outnum,outer,FKu);
            // outputForcey(outnum,outer,FKv);
            outnum++;
        }
    }      
}