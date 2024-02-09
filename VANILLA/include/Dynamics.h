#ifndef _INC_DYNAMICS_H_
#define _INC_DYNAMICS_H_

#include "Component.h"
#include "Field.h"
#include "PoissonSolver.h"
#ifdef _OPENMP
#include <omp.h>
#endif


class Dynamics {
protected:
    string                name, error_message;
    double                init_dt, dt, dx, dy, dz;
    int                   M;  // number of components
    int                   nx, ny, nz, steps;
    int                   size;  // nx*ny*nz
    vector<VectorField>   flow;
    vector<VectorField>   phi_vec;
    fstream               lobj;

    virtual void updateVelocity(void) { return; }

    void updatePhiVec(vector<ScalarField> & phi) {
        for (int p = 0; p < M; p++) {
            for (int ix = 0; ix < size; ix++) {
                phi_vec[p][ix] = phi[p].stag(ix);
            }
        }
        return;
    }

    void updateDensity(vector<ScalarField> & phi) {
        int p, ix;
        // update flow
        for (p = 0; p < M; p++) {
            for (ix = 0; ix < size; ix++) {
                flow[p][ix] = phi_vec[p][ix].mp(velocity[p][ix]);
            }
        }

        // update density
        for (p = 0; p < M; p++) {
            for (ix = 0; ix < size; ix++) {
                phi[p][ix] -= dt * flow[p].div(ix);
            }
        }
        return;
    }

    void updateDensityWithChangeDt(vector<ScalarField> & phi) {
        int     p, ix;
        double  val;

        dt = init_dt;

        // update flow
        for (p = 0; p < M; p++) {
            for (ix = 0; ix < size; ix++) {
                flow[p][ix] = phi_vec[p][ix].mp(velocity[p][ix]);
            }
        }

        // pre-update density to find 'dt'
        for (p = 0; p < M; p++) {
            for (ix = 0; ix < size; ix++) {
                val = phi[p][ix] - dt * flow[p].div(ix);
                if (val <= 0.0e0) {
                    dt = 0.9 * phi[p][ix] / flow[p].div(ix);
                }
            }
        }

        // update density
        for (p = 0; p < M; p++) {
            for (ix = 0; ix < size; ix++) {
                phi[p][ix] -= dt * flow[p].div(ix);
            }
        }

        return;
    }

    void baseSetup(string in_name,
                   Mesh in_mesh,
                   double _dt,
                   int _M) {
        name    = in_name;
        size    = in_mesh.nx*in_mesh.ny*in_mesh.nz;
        nx      = in_mesh.nx;
        ny      = in_mesh.ny;
        nz      = in_mesh.nz;
        dx      = in_mesh.dx;
        dy      = in_mesh.dy;
        dz      = in_mesh.dz;
        dt      = _dt;
        init_dt = dt;
        M       = _M;

        flow.resize(_M);
        phi_vec.resize(_M);
        velocity.resize(_M);
        for (int p = 0; p < M; p++) {
            flow[p].resize(in_mesh);
            phi_vec[p].resize(in_mesh);
            velocity[p].resize(in_mesh);
        }
        return;
    }

public:
    explicit Dynamics(string in_name = "NULL") { name = in_name; }
    ~Dynamics() {}
    string    getName(void) { return name; }

    virtual double    evolve(vector<ScalarField> & phi,
                             vector<ScalarField> & mu,
                             vector<Component> & comps) {
        return dt;
    }
    virtual double    evolveWithChangeDt(vector<ScalarField> & phi,
                                         vector<ScalarField> & mu,
                                         vector<Component> & comps) {
        return dt;
    }
    virtual void updateIsPolymer(vector<Component> & comps) {
        return;
    }
    vector<VectorField> velocity;
    vector<StressField> shear;
    vector<ScalarField> bulk;

    virtual void    logSetting(const string name) {}
    virtual bool    logger(int record, vector<ScalarField> & phi) { return false; }
    string          passError() { return error_message; }
    double          get_dt() { return dt; }
};

class TDGL : public Dynamics {
private:
public:
    TDGL() { name = "Time Depend Ginzburg-Landau"; }
    TDGL(Mesh in_mesh,
         double in_dt,
         int in_M) {
        baseSetup("Time Depend Ginzburg-Landau", in_mesh, in_dt, in_M);
        steps = 0;
    }
    ~TDGL() {}

    void updateVelocity(vector<ScalarField> & mu) {
        Vector3d    lambda, denom;
        Vector3d    zero(0.0, 0.0, 0.0);
        int         ix, p;

        for (ix = 0; ix < size; ix++) {
            lambda = zero;
            denom = zero;
            for (p = 0; p < M; p++) {
                velocity[p][ix] = -mu[p].grad(ix);
                lambda += phi_vec[p][ix].mp(velocity[p][ix]);
                denom  += phi_vec[p][ix];
            }
            for (p = 0; p < M; p++) {
                velocity[p][ix] -= lambda.dv(denom);
            }
        }
        return;
    }

    double    evolve(vector<ScalarField> & phi,
                     vector<ScalarField> & mu,
                     vector<Component> & comps) {
        updatePhiVec(phi);
        updateVelocity(mu);
        // updateFlow();  // merged updateDensity 20160526 youiida
        updateDensity(phi);

        return dt;
    }
    double    evolveWithChangeDt(vector<ScalarField> & phi,
                                 vector<ScalarField> & mu,
                                 vector<Component> & comps) {
        // same as Dynamics class
        updatePhiVec(phi);

        updateVelocity(mu);

        // same as Dynamics class
        updateDensityWithChangeDt(phi);
        steps++;
        return dt;
    }

    bool logger(int record, vector<ScalarField> & phi) {
        double    max = -1.0;
        double    min = 1.0;
        bool      flag_isnum = true;
        for (int p = 0; p < M; p++) {
            for (int ix = 0; ix < size; ix++) {
                if (phi[p][ix] > max)     max       = phi[p][ix];
                if (phi[p][ix] < min)     min       = phi[p][ix];
                if (isnan(phi[p][ix]))    flag_isnum = false;
            }
        }
        lobj << record << "," << steps << "," << dt << "," << max << "," << min << endl;

        if (!flag_isnum) {
            error_message = "density field has Inf/NaN";
        }
        steps = 0;

        return !flag_isnum;
    }
};

class ViscoElastic : public Dynamics {
private:
    // parameters
    double          d2inv, dx2inv, dy2inv, dz2inv, dtinv, ne, ne3inv, phic;
    double          dim_inv;  // dimension inverse(0.33333...)
    vector<bool>    is_polymer;

    // calc. stress working array
    StressField temp_shear;
    ScalarField temp_bulk;

    // solve poisson-eq for real_vtot
    VectorField     real_vtot;
    ScalarField     pressure;
    ScalarField     divs;
    VectorField     sorce;
    ScalarField     div_vtot;

    PoissonSolver*  solver;

    // parameters for poisson-solver iteration
    int     max_counter, min_counter, average_counter, counter1loop;
    double  init_average_counter;
    bool    average_counter_flag;
    bool    iter_flag;


    vector<Vector3d>    coeff_mat, coeff_inv;  // M*M-matrix to calculate velocity
    vector<Vector3d>    vec_vector;  // M-vector to calculate velocity
    vector<double>      x;  // working M-array
    vector<int>         indx;  // working M-array
    // {vec_vector}_i = [coeff_mat]_ij{velocity}_j
    // [coeff_inv]_ij{vec_vector}_j = [coeff_inv]_ij[coeff_mat]_jk{velocity}_k
    // [coeff_inv]_ij{vec_vector}_j = delta_ik{velocity}_k
    //
    // original equations(index 0, a represent components)
    // LHS  = \frac{\nabla \cdot \sigma_a}{\phi_a} - \frac{\nabla \cdot \sigma_0}{\phi_0}
    //       - \grad \mu_a + \grad \mu_0
    //
    // RHS  = X_a \vecv_a - X_0 \vecv_0 + \sum_j(-Y_aj + Y_bj)\vecv_j
    // X_a  = \frac{Z^m + N\zeta_a^m + \zeta_l^p}{\phi_a}
    // Y_aj = \frac{\zeta_j^m + \zeta_a^m + Z^{-1}\zeta_a^p\zeta_j^p}{\phi_a}


    // velocity matrix(RHS 1st or 2nd term) coefficient
    Vector3d    getX(vector<Component> & comps,
                     int p,
                     int ix,
                     Vector3d zm) {
        return (zm + M * getMz(p, ix, comps) + getPz(p, ix, comps) ).dv(phi_vec[p][ix]);
    }

    void    setZ(Vector3d & inv,
                 Vector3d & zm,
                 vector<Component> & comps,
                 int ix) {
        Vector3d zero(0.0, 0.0, 0.0);
        Vector3d zp;

        zp = zero;  // sum zeta^p = Z -> to get Z^{-1}
        zm = zero;  // sum zeta^m = Z^m
        for (int p = 0; p < M; p++) {
            zp += phi_vec[p][ix] * comps[p].pz();
            zm += phi_vec[p][ix] * comps[p].mz();
        }
        if (zp.length2() > 1e-20) {
            inv = Vector3d(1.0, 1.0, 1.0).dv(zp);
        } else {
            inv = zero;
        }
        return;
    }

    // \zeta_{index-p}^p
    Vector3d    getPz(int p,
                      int ix,
                      vector<Component> & comps) {
        return comps[p].pz() * phi_vec[p][ix];
    }

    // \zeta_{index-p}^m
    Vector3d    getMz(int p,
                      int ix,
                      vector<Component> & comps) {
        return comps[p].mz() * phi_vec[p][ix];
    }

    // Y_p1p2
    Vector3d    getY(int p1,
                     int p2,
                     int ix,
                     vector<Component> & comps,
                     Vector3d inv) {
        return (getMz(p2, ix, comps)
                +getMz(p1, ix, comps)
                +inv.mp(getPz(p1, ix, comps), getPz(p2, ix, comps)) ).dv(phi_vec[p1][ix]);
    }


    Stress getPhipTensor(vector<ScalarField> & phi,
                         int ix) {
        Stress phip = Stress(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0.0);
        for (int p = 0; p < M; p++) {
            if (is_polymer[p]) {
                phip += phi[p].toTensor(ix);
            }
        }
        return phip;
    }

    // tensor calc sigma*exp(-dt/tau)
    Stress decay(Stress s, Stress tau) {
        return Stress(s.xx * exp(-dt / tau.xx),
                      s.xy * exp(-dt / tau.xy),
                      s.xz * exp(-dt / tau.xz),
                      s.yx * exp(-dt / tau.yx),
                      s.yy * exp(-dt / tau.yy),
                      s.yz * exp(-dt / tau.yz),
                      s.zx * exp(-dt / tau.zx),
                      s.zy * exp(-dt / tau.zy),
                      s.zz * exp(-dt / tau.zz) );
    }

    // tensor calc sigma( nabla(v) + nabla(v)^T )
    Stress grow(Stress s, Stress nv) {
        return Stress(s.xx * (nv.xx + nv.xx),
                      s.xy * (nv.xy + nv.yx),
                      s.xz * (nv.xz + nv.zx),
                      s.yx * (nv.yx + nv.xy),
                      s.yy * (nv.yy + nv.yy),
                      s.yz * (nv.yz + nv.zy),
                      s.zx * (nv.zx + nv.xz),
                      s.zy * (nv.zy + nv.yz),
                      s.zz * (nv.zz + nv.zz) );
    }


    // to modify stress model change this function
    void    updateStress(vector<ScalarField> & phi,
                         vector<Component> & comps) {
        // update stress-parameters
        double    gb, tb, phip, nphi3, n3;
        Stress    unit(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        Stress    gs, ts, ppt, nppt3, nv;  // G_s, tau_s, phi_polymer_tensor, n2_phi_polymer_tensor, nabla v_pol

        int       p1, ix;
        // update stress
        for (p1 = 0; p1 < M; p1++) {
            if (is_polymer[p1]) {
                for (ix = 0; ix < size; ix++) {
                    ppt   = getPhipTensor(phi, ix);
                    phip  = ppt.xx;
                    // change N, phi dependence 20161019 youiida
                    n3    = comps[p1].N * comps[p1].N * comps[p1].N * ne3inv;
                    nphi3 = n3 * phip * phip;
                    nppt3 = n3 * ppt.square();

                    gs = comps[p1].get_gs() * ppt.square();
                    ts = comps[p1].get_ts() * nppt3;
                    gb = comps[p1].get_gb() * step(phip, phic);
                    tb = comps[p1].get_tb() * nphi3;

                    nv = velocity[p1].nabla(ix);  // nabla(v)

                    temp_shear[ix]  = decay(shear[p1][ix], ts)
                        + (- vectorDotNablaStress(velocity[p1], shear[p1], ix)
                           + stressDotNablaVectorSymmetry(shear[p1], velocity[p1], ix)
                           + grow(gs, nv) )*dt;
                    temp_shear[ix] -= temp_shear[ix].trace() * dim_inv * unit;

                    temp_bulk[ix]   = bulk[p1][ix] * exp(-dt / tb)
                        + (- vectorDotNablaScalar(velocity[p1], bulk[p1], ix)
                           + gb * velocity[p1].div(ix) )*dt;
                }
                // Bug-fix using temp_shear, temp_bulk
                // 20160802 youiida
                for (ix = 0; ix < size; ix++) {
                    shear[p1][ix] = temp_shear[ix];
                    bulk[p1][ix]  = temp_bulk[ix];
                }
            }
        }  // end update stress
        return;
    }

    void    updateVelocity(vector<ScalarField> & phi,
                           vector<ScalarField> & mu,
                           vector<Component> & comps) {
        int         ix, p1, p2;
        Vector3d    zero(0.0, 0.0, 0.0);
        Vector3d    zinv, zm, x0, mu0, force, divs0;

        // update stress
        // to modify stress model change this function
        updateStress(phi, comps);


        // update sorce, divsorce, div_vtot
        for (ix = 0; ix < size; ix++) {
            force    = zero;
            for (p1 = 0; p1 < M; p1++) {
                force +=  phi_vec[p1][ix].mp(mu[p1].grad(ix))  // phi grad mu
                        - nablaDotStrss(shear[p1], bulk[p1], ix);
            }
            sorce[ix] = - force;
        }

        for (ix = 0; ix < size; ix++) {
            divs[ix]     = sorce.div(ix);
            div_vtot[ix] = real_vtot.div(ix) * dtinv;
        }
        // end update sorce, divsorce, div_vtot

        // solve pressure
        iter_flag = solver->solve(pressure, divs, div_vtot, counter1loop);
        if (iter_flag) {
            warning("cannot solve pressure!");
        }
        // end solve pressure

        // iteration counter infomation
        if (max_counter < counter1loop) {
            max_counter = counter1loop;
        }
        if (min_counter > counter1loop) {
            min_counter = counter1loop;
        }
        average_counter += counter1loop;


        // update real_vtot
        for (ix = 0; ix < size; ix++) {
            real_vtot[ix] += dt * (sorce[ix] - pressure.grad(ix));
        }



        // update velocity field <main>
        for (ix = 0; ix < size; ix++) {
            // set coeff_mat and vec_vector

            // p1=0 -> vtot elements, 1st row
            for (p2 = 0; p2 < M; p2++) {
                coeff_mat[p2] = phi_vec[p2][ix];  // p1 = 0
            }
            vec_vector[0] = real_vtot[ix];

            // re-use parameters
            setZ(zinv, zm, comps, ix);
            x0    = getX(comps, 0, ix, zm);
            mu0   = mu[0].grad(ix);
            divs0 = nablaDotStrss(shear[0], bulk[0], ix).dv(phi_vec[0][ix]);

            // p1=1-M -> A-matrix, B-vector, 2nd-M-th column
            for (p1 = 1; p1 < M; p1++) {
                // LHS, vec_vector
                vec_vector[p1] =   nablaDotStrss(shear[p1], bulk[p1], ix).dv(phi_vec[p1][ix])
                                 - divs0
                                 - mu[p1].grad(ix)
                                 + mu0;

                // -Y_aj+Y_bj
                for (p2 = 1; p2 < M; p2++) {
                    coeff_mat[p1*M + p2] = - getY(p1, p2, ix, comps, zinv)
                                           + getY(0,  p2, ix, comps, zinv);
                }

                // -X_b
                coeff_mat[p1*M]       = - x0;  // p2=0

                // X_a
                coeff_mat[p1*M + p1] += getX(comps, p1, ix, zm);
            }

            // solve coeff_inv
            setInvMatrix(coeff_mat, coeff_inv, x, indx, M);

            // fix constraint using lagrange multiplyer(lambda->vtot,denom) 20160229 youiida
            // -> makes bug?? 20160426 youiida
            // lagrange multiplyer conflicts to MAC-method and causes double correction
            // 20160614 youiida
            // update velocity fields
            for (p1 = 0; p1 < M; p1++) {
                velocity[p1][ix] = zero;
                for (p2 = 0; p2 < M; p2++) {
                    velocity[p1][ix] += coeff_inv[M*p1 + p2].mp(vec_vector[p2]);
                }
                // denom += phi_vec[p1][ix];
            }
            // for (p1 = 0; p1 < M; p1++) {
            //    velocity[p1][ix] -= real_vtot[ix].dv(denom);
            // }
        }
        // end update velocity field <main>


        return;
    }

    void updateIsPolymer(vector<Component> & comps) {
        for (int p = 0; p < M; p++) {
            (ne < comps[p].N) ? is_polymer[p] = true : is_polymer[p] = false;
            // condition ? true : false
        }
        return;
    }


public:
    ViscoElastic() {}
    ~ViscoElastic() {}
    ViscoElastic(Mesh in_mesh,
                 double in_dt,
                 int in_M,
                 double in_ne,
                 double in_phic,
                 double in_converge,
                 vector<Component> & comps,
                 string in_model,
                 double in_omega) {
        baseSetup("Visco Elastic", in_mesh, in_dt, in_M);
        ne     = in_ne;
        phic   = in_phic;
        ne3inv = 1.0/(ne*ne*ne);

        if (nx > 1 && ny > 1 && nz > 1) {
            dim_inv = 1.0 / 3.0;
        } else {
            // dim_inv = 0.5;
            warning("invalid simulation box size\nonly support 3D");
            // quasi-2d simulation causes invalid diffusion constant
            // 20160318 youiida
        }

        shear.resize(in_M);
        bulk.resize(in_M);
        for (int p = 0; p < M; p++) {
            shear[p].resize(in_mesh);
            bulk[p].resize(in_mesh);
        }

        temp_shear.resize(in_mesh);
        temp_bulk.resize(in_mesh);

        real_vtot.resize(in_mesh);


        coeff_mat.resize(M*M);
        coeff_inv.resize(M*M);
        vec_vector.resize(M);

        is_polymer.resize(M);
        x.resize(M);
        indx.resize(M);

        // solve pressure paramters
        pressure.resize(in_mesh);
        divs.resize(in_mesh);
        sorce.resize(in_mesh);
        div_vtot.resize(in_mesh);

        dx2inv = 1.0 / (dx*dx);
        dy2inv = 1.0 / (dy*dy);
        dz2inv = 1.0 / (dz*dz);
        d2inv  = 1.0 / (-2.0*(dx2inv + dy2inv + dz2inv));
        dtinv  = 1.0 / dt;

        // set solver method
        if (in_model == "MAC-Gauss_Seidel") {
            solver = new GaussSeidel(size, d2inv, in_mesh, in_converge);
        } else if (in_model == "MAC-Jacobi") {
            caution("no-check");
            solver = new Jacobi(size, d2inv, in_mesh, in_converge);
        } else if (in_model == "MAC-SOR") {
            solver = new SOR(size, d2inv, in_mesh, in_converge, in_omega);
            // else if (in_model == "CGNR") {
            //    caution("no-check");
            //    solver = new CGNR(size, d2inv, in_mesh, in_converge);
            // }
            // else if (in_model == "OSEEN") {
            //    warning("under constraction");
            //    solver = new OSEEN();
            // }
        } else {
            warning("invalid solver method, ViscoElastic");
        }


        // field initialize
        Stress      zero_stress(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        Vector3d    zero_vec(0.0, 0.0, 0.0);
        for (int p = 0; p < M; p++) {
            for (int ix = 0; ix < size; ix++) {
                shear[p][ix]    = zero_stress;
                bulk[p][ix]     = 0.0;
                velocity[p][ix] = zero_vec;
                pressure[ix]    = 0.0;
                real_vtot[ix]   = zero_vec;
            }
        }

        // iteration logs
        max_counter          = 0;
        min_counter          = INT_MAX;
        average_counter      = 0;
        average_counter_flag = false;
        steps                = 0;

        updateIsPolymer(comps);
    }

    double    evolve(vector<ScalarField> & phi,
                     vector<ScalarField> & mu,
                     vector<Component> & comps) {
        // same as Dynamics class
        updatePhiVec(phi);

        // ViscoElastic original
        updateVelocity(phi, mu, comps);

        // same as Dynamics class
        // updateFlow();//merged updateDensity 20160526 youiida
        updateDensity(phi);

        steps++;
        return dt;
    }

    double    evolveWithChangeDt(vector<ScalarField> & phi,
                                 vector<ScalarField> & mu,
                                 vector<Component> & comps) {
        // same as Dynamics class
        updatePhiVec(phi);

        // ViscoElastic original
        updateVelocity(phi, mu, comps);

        // same as Dynamics class
        updateDensityWithChangeDt(phi);

        steps++;
        return dt;
    }

    void logSetting(const string name) {
        string rep = replace2(name, ".bdf", ".udf", "_dynamics_log.csv");
        lobj.open(rep.c_str(), ios::out);
        lobj << "record,averaged iteration,maximum iteration,minimum iteration,steps,dt,max_phi,min_phi\n";
        return;
    }
    bool logger(int record, vector<ScalarField> & phi) {
        double    max = -1.0;
        double    min = 1.0;
        bool      flag_isnum = true;
        for (int p = 0; p < M; p++) {
            for (int ix = 0; ix < size; ix++) {
                if (phi[p][ix] > max)     max = phi[p][ix];
                if (phi[p][ix] < min)     min = phi[p][ix];
                if (isnan(phi[p][ix]))    flag_isnum = false;
            }
        }
        lobj << record << "," << double(average_counter)/double(steps) << ","
             << max_counter << "," << min_counter << "," << steps << "," << dt << ","
             << max << "," << min << endl;

        if (!flag_isnum) {
            error_message = "density field has Inf/NaN";
        }
        max_counter     = 0;
        average_counter = 0;
        min_counter     = INT_MAX;
        steps           = 0;

        return !flag_isnum;
    }
};

#endif
