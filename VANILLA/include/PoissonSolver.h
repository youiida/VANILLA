#ifndef _INC_POISSON_SOLVER_H_
#define _INC_POISSON_SOLVER_H_

#include "Field.h"
#ifdef _OPENMP
#include <omp.h>
#endif

class PoissonSolver {
protected:
    string    name;
    int       size;
    double    d2inv;
    double    convergence;
public:
    PoissonSolver() {name="base_class";};
    ~PoissonSolver() {};
    virtual bool solve(ScalarField & p, ScalarField & ds, ScalarField & dv, int & counter) {
        caution("base solver!");
        return false;
    };
};

class Jacobi : public PoissonSolver {
private:
    ScalarField    temp_p;//Jacobi
public:
    Jacobi() {name="Jacobi";};
    Jacobi(int in_size,
           double in_d2inv,
           Mesh in_mesh,
           double in_convergence) {
        name        = "Jacobi";
        size        = in_size;
        d2inv       = in_d2inv;
        temp_p.resize(in_mesh);
        convergence = in_convergence;
    };
    ~Jacobi() {};
    //20160810 youiida
    //solve pressure using MAC(JACOBI)-method
    bool solve(ScalarField & p,
               ScalarField & ds,
               ScalarField & dv,
               int & counter) {
        double maxp, cerror, merror, prev;
        counter = 0;
        maxp    = -1.0;


        while (1) {
            //update pressure
            merror = -1.0;
            //OPENMP parallel for...
            for (int ix = 0; ix < size; ix++) {
                prev       = p[ix];
                temp_p[ix] = d2inv * (-p.macValues(ix) + ds[ix] + dv[ix]);//Jacobi
                if (maxp < fabs(temp_p[ix])) {
                    maxp = fabs(temp_p[ix]);
                }
                cerror = fabs(temp_p[ix] - prev) / maxp;
                if (merror < cerror) merror = cerror;
            }

            for (int ix = 0; ix < size; ix++) {
                p[ix] = temp_p[ix];
            }

            counter++;
            // not converge
            if (counter == INT_MAX) {
                return true;
                //break;
            }
            // converge
            if (merror < convergence) {
                return false;
                //break;
            }
        }
        return true;
    };

};

class GaussSeidel : public PoissonSolver {
public:
    GaussSeidel() {name="Gauss-Seidel";};
    GaussSeidel(int in_size,
        double in_d2inv,
        Mesh in_mesh,
        double in_convergence) {
        name        = "Gauss-Seidel";
        size        = in_size;
        d2inv       = in_d2inv;
        convergence = in_convergence;
    };
    ~GaussSeidel() {};
    //same as MAC(SOR) omega=1
    //20160810 youiida
    //solve pressure using MAC(GS)-method
    bool solve(ScalarField & p,
        ScalarField & ds,
        ScalarField & dv,
        int & counter) {
        double maxp, cerror, merror, prev;
        counter = 0;
        maxp    = -1.0;
        while (1) {
            //update pressure
            merror = -1.0;
            for (int ix = 0; ix < size; ix++) {
                prev = p[ix];
                p[ix] = d2inv * (-p.macValues(ix) + ds[ix] + dv[ix]);//Gauss-Seidel
                if (maxp < fabs(p[ix])) {
                    maxp = fabs(p[ix]);
                }
                cerror = fabs(p[ix] - prev) / maxp;
                if (merror < cerror) merror = cerror;
            }
            counter++;
            // not converge
            if (counter == INT_MAX) {
                return true;
                //break;
            }
            // converge
            if (merror < convergence) {
                return false;
                //break;
            }
        }
        return true;
    };

};

//Successive Over Relaxation
class SOR : public PoissonSolver {
private:
    double omega_SOR;
public:
    SOR() {};
    SOR(int in_size,
        double in_d2inv,
        Mesh in_mesh,
        double in_convergence,
        double in_omega_SOR) {
        name        = "SOR";
        size        = in_size;
        d2inv       = in_d2inv;
        convergence = in_convergence;
        omega_SOR   = in_omega_SOR;
        if (omega_SOR >= 2.0 || omega_SOR < 0.0) {
            warning("0.0 < omega_SOR < 2.0");
        }
    };
    ~SOR() {};
    //solve pressure using MAC(SOR,GS)-method
    bool solve(ScalarField & p,
               ScalarField & ds,
               ScalarField & dv,
               int & counter) {
        double  maxp, cerror, merror, prev;
        counter = 0;
        maxp    = -1.0;

        while (1) {
            //update pressure
            merror = -1.0;
            for (int ix = 0; ix < size; ix++) {
                prev = p[ix];
                p[ix] = p[ix] + omega_SOR * (d2inv * (-p.macValues(ix) + ds[ix] + dv[ix]) - p[ix]);//SOR
                if (maxp < fabs(p[ix])) {
                    maxp = fabs(p[ix]);
                }
                cerror = fabs(p[ix] - prev) / maxp;
                if (merror < cerror) merror = cerror;
            }
            counter++;

            // do not converge
            if (counter == INT_MAX) {
                return true;
                //break;
            }

            // converged
            if (merror < convergence) {
                return false;
                //break;
            }
        }
        return true;
    };


};

//conjugate gradient 
class CGNR : public PoissonSolver {
private:
    //CGNR-method
    vector<double>  cg_p;//working array
    vector<double>  cg_q;//working array
    vector<double>  cg_r;//working array
    vector<double>  cg_b;//NOT Positive-defined vector
    //vector<double>            cg_a;//NOT Positive-defined matrix
    vector<double>  cg_ta;//to convert positive-defined matrix
    vector<double>  cg_m;//size*size matrix
    double ip(ScalarField & v1,
        ScalarField & v2) {
        int size   = v1.size();
        double val = 0.0;
        for (int i = 0; i<size; i++)val += v1[i] * v2[i];
        return val;
    }

    vector<double> mv(ScalarField & m, ScalarField & v) {
        int range = v.size();
        vector<double> val(range, 0.0);
        for (int i = 0; i < range; i++) {
            for (int j = 0; j < range; j++) {
                val[i] += m[i*range + j] * v[j];
            }
        }
        return val;
    }
    double ip(vector<double> & v1, vector<double> & v2) {
        int size = (int)v1.size();
        double val = 0.0;
        for (int i = 0; i<size; i++)val += v1[i] * v2[i];
        return val;
    }

    vector<double> mv(vector< double > & m, vector<double> & v) {
        int size = (int)v.size();
        vector<double> val(size, 0.0);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                val[i] += m[i*size + j] * v[j];
            }
        }
        return val;
    }

public:
    CGNR() {};
    CGNR(int in_size,
         double in_d2inv,
         Mesh in_mesh,
         double in_convergence) {
        name        = "CGNR";
        size        = in_size;
        d2inv       = in_d2inv;
        convergence = in_convergence;

        ////solve pressure, CGNR-method
        //cg_p.resize(size);
        //cg_q.resize(size);
        ////CGNR working array
        //vector<double>        cg_a;
        //cg_b.resize(size);
        //cg_a.resize(size*size, 0.0);
        //cg_ta.resize(size*size);
        //
        ////set a
        //for (int i = 0; i < size; i++) {
        //        cg_a[i] = ;
        //        cg_a[mx(i)] = ;
        //        cg_a[my(i)] = ;
        //        cg_a[mz(i)] = ;
        //        cg_a[px(i)] = ;
        //        cg_a[py(i)] = ;
        //        cg_a[pz(i)] = ;
        //}
        ////set a^t
        //for (int i = 0; i < size; i++) {
        //    for (int j = 0; j < size; j++) {
        //        cg_ta[i*size+j] = cg_a[j*size+i];
        //    }
        //}
        ////set m
        //for (int i = 0; i < size; i++) {
        //    for (int j = 0; j < size; j++) {
        //        cg_m[i*size + j] = 0.0;
        //        for (int k = 0; k < size; k++) {
        //            cg_m[i*size+j] = cg_m[i*size + j] + cg_ta[i*size+k]*cg_a[k*size+j];
        //        }
        //    }
        //}
    };
    ~CGNR() {};
    //solve pressure using CG: Conjugate Gradient method
    //                   + NR: Convert normal equations to positive-(semi)definite equations(DEFINED ON CONSTRACTOR)
    // A*x = b'
    // A^t*A*x = A^t*b' -> M*x = b
    // A:discrete-laplacian operator -> A and A^t is defined on constractor
    // x:pressure
    // b:divs, div_vtot
    bool solve(ScalarField & p,
        ScalarField & ds,
        ScalarField & dv,
        int & counter) {
        //cg-method
        double    rho, rho_m, rho_p;//rho_i-1, rho_i-2, rho_i+1
        double    beta, alpha, inner_product;
        bool      flag = false;


        //CGNR-method
        //make b = A^t*temp_b
        for (int i = 0; i < size; i++) {
            cg_b[i] = ds[i] + dv[i];
        }
        cg_r = mv(cg_ta, cg_b);
        //for (int i = 0; i < size; i++) {
        //    for (int j = 0; j < size; j++) {
        //        cg_r[i] = cg_ta[i*size+j] * cg_b[j];
        //    }
        //}

        //compute r0
        //cg_r = cg_b;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                cg_r[i] -= cg_m[i*size + j] * p[j];
            }
        }
        cg_p = cg_r;
        rho  = ip(cg_r, cg_r);

        counter = 0;
        while (1) {

            cg_q          = mv(cg_m, cg_p);
            inner_product = ip(cg_p, cg_q);
            alpha         = rho / inner_product;

            rho_p = 0.0;
            for (int i = 0; i < size; i++) {
                p[i] = p[i] + alpha*cg_p[i];
                cg_r[i] =  cg_r[i] - alpha*cg_q[i];
                rho_p   += cg_r[i] * cg_r[i];
            }
            if (rho_p < convergence) break;

            rho_m = rho;
            rho   = rho_p;
            //rho = ip(cg_r, cg_r);
            beta  = rho / rho_m;
            for (int i = 0; i < size; i++) {
                cg_p[i] = cg_r[i] + beta*cg_p[i];
            }

            counter++;
        }

        return true;
    };

};

//solve using Oseen tensor
class Oseen : public PoissonSolver {
private:
    //solve real_vtot using MAC-method, do not use oseen tensor
    //StressField        real_oseen;
    //StressField        oseen;
    //CVectorField    imag_vtot, imag_force;//to calculate vtot
    //VectorField    real_force;//to calculate vtot

    // {imag_vtot} = [oseen]{imag_force}
    // [oseen] = -(I-kk/k^2)/eta*k^2
    // {imag_force} = fft({real_force})
    // {real_force} = sum(phi*grad(mu)) - nabla dot sigma
public:
    Oseen() {};
    ~Oseen() {};
    //Oseen() {
        //real_force.resize(in_mesh);
        //imag_force.resize(in_mesh);
        //imag_vtot.resize(in_mesh);
        //real_force.fftSetup();
        //imag_force.fftSetup();
        //real_vtot.fftSetup();
        //imag_vtot.fftSetup();
        //real_oseen.resize(in_mesh);
        //oseen.resize(in_mesh);
        //real_oseen.fftSetup();


        //setup oseen tensor
        //Vector3d    r;
        //int            ix, iy, iz, c;


        //double        dqx = PI2 / (nx*dx);
        //double        dqy = PI2 / (ny*dy);
        //double        dqz = PI2 / (nz*dz);

        //double        dou_x, dou_y, dou_z;

        //double        s2, cfx, cfy, cfz;
        //int            plx, ply, plz;
        //for (iz = 0; iz < nz; iz++) {
        //    dou_z = iz*dqz;
        //    plz = iz;
        //    if (iz >(nz / 2)) {
        //        dou_z -= dqz*nz;
        //        plz -= nz;
        //    }
        //    cfz = sin(PI*double(plz) / double(nz)) / (PI*double(plz));
        //    if (iz == 0) {
        //        cfz = double(nz);
        //    }

        //    for (iy = 0; iy < ny; iy++) {
        //        dou_y = iy*dqy;
        //        ply = iy;
        //        if (iy >(ny / 2)) {
        //            dou_y -= dqy*ny;
        //            ply -= ny;
        //        }
        //        cfy = sin(PI*double(ply) / double(ny)) / (PI*double(ply));
        //        if (iy == 0) {
        //            cfy = double(ny);
        //        }

        //        for (ix = 0; ix < nx; ix++) {
        //            dou_x = ix*dqx;
        //            plx = ix;
        //            if (ix >(nx / 2)) {
        //                dou_x -= dqx*nx;
        //                plx -= nx;
        //            }
        //            cfx = sin(PI*double(plx) / double(nx)) / (PI*double(plx));
        //            if (ix == 0) {
        //                cfx = double(nx);
        //            }
        //            r.set(0.5*dou_x, 0.5*dou_y, 0.5*dou_z);
        //            //s2 = r.length2();
        //            s2 = sin(r.x)*sin(r.x) + sin(r.y)*sin(r.y) + sin(r.z)*sin(r.z);
        //            //s2 = sin(r.x) * sin(r.x);
        //            //s = r.length();

        //            c = ix + iy*nx + iz*nx*ny;

        //            //real_oseen[c].xx = (1.0 - (r.x*r.x) / (s*s)) / (4.0*PI2*s);
        //            //real_oseen[c].xy = (0.0 - (r.x*r.y) / (s*s)) / (4.0*PI2*s);
        //            //real_oseen[c].xz = (0.0 - (r.x*r.z) / (s*s)) / (4.0*PI2*s);

        //            //real_oseen[c].yy = (1.0 - (r.y*r.y) / (s*s)) / (4.0*PI2*s);
        //            //real_oseen[c].yz = (0.0 - (r.y*r.z) / (s*s)) / (4.0*PI2*s);

        //            //real_oseen[c].zz = (1.0 - (r.z*r.z) / (s*s)) / (4.0*PI2*s);

        //            //real_oseen[c].yx = real_oseen[c].xy;
        //            //real_oseen[c].zx = real_oseen[c].xz;
        //            //real_oseen[c].zy = real_oseen[c].yz;

        //            //oseen[c].xx = cfx*cfy*cfz*(1.0 - r.x*r.x / s2) / s2;
        //            //oseen[c].xy = cfx*cfy*cfz*(0.0 - r.x*r.y / s2) / s2;
        //            //oseen[c].xz = cfx*cfy*cfz*(0.0 - r.x*r.z / s2) / s2;

        //            //oseen[c].yx = oseen[c].xy;
        //            //oseen[c].yy = cfx*cfy*cfz*(1.0 - r.y*r.y / s2) / s2;
        //            //oseen[c].yz = cfx*cfy*cfz*(0.0 - r.y*r.z / s2) / s2;

        //            //oseen[c].zx = oseen[c].xz;
        //            //oseen[c].zy = oseen[c].yz;
        //            //oseen[c].zz = cfx*cfy*cfz*(1.0 - r.z*r.z / s2) / s2;

        //            oseen[c].xx = (1.0 - sin(r.x)*sin(r.x)/s2) / s2;
        //            oseen[c].xy = (    - sin(r.x)*sin(r.y)/s2) / s2;
        //            oseen[c].xz = (    - sin(r.x)*sin(r.z)/s2) / s2;

        //            oseen[c].yx = oseen[c].xy;
        //            oseen[c].yy = (1.0 - sin(r.y)*sin(r.y)/s2) / s2;
        //            oseen[c].yz = (    - sin(r.y)*sin(r.z)/s2) / s2;

        //            oseen[c].zx = oseen[c].xz;
        //            oseen[c].zy = oseen[c].yz;
        //            oseen[c].zz = (1.0 - sin(r.z)*sin(r.z)/s2) / s2;
        //        }
        //    }
        //}
        //oseen[0] = zero_stress;//index-0 is q=0
        //real_oseen.fft_scale(oseen);
        //end setup oseen tensor


    //}
};

#endif