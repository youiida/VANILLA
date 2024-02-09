#ifndef _INC_FREE_ENERGY_H_
#define _INC_FREE_ENERGY_H_
#include "Field.h"
#include "Component.h"
#ifdef _OPENMP
#include <omp.h>
#endif

class FreeEnergy {
protected:
    int               M;//number of total components
    int               nx, ny, nz;
    int               size;
    double            vinv;
    vector<double>    chi;
    string            name;
    void    setup(string _name,
                  Mesh & _mesh,
                  int _M,
                  vector<double> & _chi) {
        name = _name;
        vinv = 1.0 / (_mesh.nx * _mesh.ny * _mesh.nz * _mesh.dx * _mesh.dy * _mesh.dz);
        nx   = _mesh.nx;
        ny   = _mesh.ny;
        nz   = _mesh.nz;

        M    = _M;
        size = _mesh.nx * _mesh.ny * _mesh.nz;

        chi  = _chi;
        return;
    }
public:
    FreeEnergy(string _name = "NULL") { name = _name; };
    FreeEnergy(string _name, 
               Mesh &_mesh, 
               int _M,
               vector<double> & _chi) {
        setup(_name, _mesh, _M, _chi);
    };
    ~FreeEnergy() {};
    virtual double  get(vector<Component> & comps,
                        vector<ScalarField> & phi) { return PI2; };
    virtual void update(vector<Component> & comps,
                        vector<ScalarField> & phi,
                        vector<ScalarField> & mu) { return; };
    void    setChi(vector<double> in_chi) { chi = in_chi; };
    //definition of chi-array
    //chi[0][0] = chi[0], chi[0][1] = chi[1], ..., chi[0][M-1] = chi[M-1]
    //chi[i][j] = chi[i*M+j]
    virtual void updateVertex(vector<Component> & comps) { return; };
    string    get_name() { return name; };
};

class FHdG : public FreeEnergy {
public:
    FHdG(Mesh &_mesh, 
         int _M,
         vector<double> & _chi) {
        setup("Flory-Huggins-deGennes", _mesh, _M, _chi);
    };
    ~FHdG() {};
    double    get(vector<Component> & comps,
                  vector<ScalarField> & phi) {
        double  fe = 0.0;
        double  uni_entropy  = 0.0;
        double  uni_enthalpy = 0.0;
        double  ninv;

        for (int p = 0; p < M; p++) {
            ninv         = 1.0 / comps[p].N;
            uni_entropy += ninv * comps[p].vf * log(comps[p].vf);
            for (int p2 = 0; p2 < M; p2++) {
                uni_enthalpy += chi[M*p + p2] * comps[p].vf * comps[p2].vf;
            }

            for (int ix = 0; ix < size; ix++) {
                fe += ninv*(phi[p][ix] * log(phi[p][ix]))
                    + comps[p].get_tension() * (phi[p].grad(ix)).length2();
                for (int p2 = 0; p2 < M; p2++) {
                    fe += chi[M*p + p2] * phi[p][ix] * phi[p2][ix];
                }
            }
        }

        return fe*vinv - (uni_entropy + uni_enthalpy);
    };
    void    update(vector<Component> & comps, 
                   vector<ScalarField> & phi,
                   vector<ScalarField> & mu) {
        double ninv;
        for (int p = 0; p < M; p++) {
            ninv = 1.0 / comps[p].N;
            for (int ix = 0; ix < size; ix++) {
                mu[p][ix] = ninv * (log(phi[p][ix]) + 1.0)
                            - 2.0 * comps[p].get_tension() * phi[p].lap(ix);
                for (int p2 = 0; p2 < M; p2++) {
                    mu[p][ix] += chi[M*p + p2] * phi[p2][ix];
                }
            }
        }
        return;
    };

};

class RPA : public FreeEnergy {
private:
    double    alpha;
    double    dqx, dqy, dqz;
    
    //calc inverse matrix
    double            det;
    vector<double>    x;
    vector<int>       indx;

    vector<CScalarField>    phiq;
    vector<ScalarField>     vertex;//GammaKK'(q)   vertex[K*M+K'][iq]
    CScalarField            imag2nd;//sum_K' GammaKK'(q)*phiqK'(q)   imag2nd[iq] <- working array
    ScalarField             real2nd;//real2nd = fft(imag2nd)   real2nd[ix] <- working array
    //vector<ScalarField>        sq;// real correlation, GammaKK'(q) = (sq)^-1   sq[K*M+K'][iq]
    vector<double>          matrix;//work matrix for sq,   matrix[K*M+K']
    vector<double>          invmat;//work matrix for sq,   invmat = (matrix)^-1

    double    pc(const int q,
                 const int size) {//periodic correction
        if (q != 0) {
            return double(size) * sin(double(q)*PI / double(size)) / (double(q)*PI);
        } else {
            return 1.0;
        }
    };
    double    idealCorrelation(vector<Component> & comps,
                               int qx,
                               int qy,
                               int qz,
                               int p1,
                               int p2) {//ideal correlation for linear polymer
        //different polymer
        if (p1 != p2)return 0.0;

        //same polymer(linear)
        double    q2 = qx*dqx*qx*dqx + qy*dqy*qy*dqy + qz*dqz*qz*dqz;
        double    N  = comps[p1].N;
        double    b  = comps[p1].get_b();
        double    x2 = N*b*b*q2 / 6.0;
        if (x2 < 1e-10) {
            return N;
        } else {
            return 2.0*N*(exp(-x2) + x2 - 1.0) / (x2*x2);
        }
    };

    //debug functions and parameters
    fstream lobj;
    void    output(vector<double> matrix, int qx, int qy, int qz) {
        lobj << qx*dqx << "\t" << qy*dqy << "\t" << qz*dqz;
        for (int p1 = 0; p1 < (int)matrix.size(); p1++) {
            lobj << "\t" << matrix[p1];
        }
        lobj << endl;
        return;
    };
    void    lobjReturn(void){lobj<<endl;};
    void    output_ic(vector<Component> & comps, int qx, int qy, int qz) {
        lobj << qx*dqx*qx*dqx + qy*dqy*qy*dqy + qz*dqz*qz*dqz;
        for (int p1 = 0; p1 < (int)comps.size(); p1++) {
            for (int p2 = 0; p2 < (int)comps.size(); p2++) {
                lobj << "\t" << idealCorrelation(comps, qx, qy, qz, p1, p2);
            }
        }
        lobj << endl;
    };

public:
    RPA() {};
    RPA(Mesh & _mesh, 
        int _M, 
        //double _alpha,
        vector<Component> & comps,
        vector<double> & _chi) {
        setup("Random_Phase_Approximation", _mesh, _M, _chi);

        lobj.open("test.dat",ios::out);

        alpha = 1e12;
        dqx = PI2 / (_mesh.nx*_mesh.dx);
        dqy = PI2 / (_mesh.ny*_mesh.dy);
        dqz = PI2 / (_mesh.nz*_mesh.dz);

        phiq.resize(M);
        for (int p = 0; p < M; p++) {
            phiq[p].resize(_mesh);
            phiq[p].fftSetup();
        }
        x.resize(M);
        indx.resize(M);

        imag2nd.resize(_mesh);
        real2nd.resize(_mesh);

        imag2nd.fftSetup();

        //sq.resize(M*M);
        matrix.resize(M*M);
        invmat.resize(M*M);
        vertex.resize(M*M);
        for (int p = 0; p < M*M; p++) {
        //    sq[p].resize(_mesh);
            vertex[p].resize(_mesh, 0.0);
        }
        updateVertex(comps);//bug?? 20160121 youiida
    }
    ~RPA() {};
    double    get(vector<Component> & comps,
                  vector<ScalarField> & phi) {
        double    fe           = 0.0;
        double    entropy      = 0.0;
        double    uni_entropy  = 0.0;
        double    uni_enthalpy = 0.0;
        double    diff, ninv;
        for (int p1 = 0; p1 < M; p1++) {
            ninv = 1.0 / comps[p1].N;
            for (int p2 = 0; p2 < M; p2++) {
                for (int ix = 0; ix < size; ix++) {
                    fe +=  0.5 * vertex[p1*M + p2][ix] * abs(phiq[p1][ix] * phiq[p2][ix]);
                }
                uni_enthalpy += chi[M*p1 + p2] * comps[p1].vf * comps[p2].vf;
            }
            for (int ix = 0; ix < size; ix++) {
                diff     = phi[p1][ix] - comps[p1].vf;
                entropy += (phi[p1][ix] * log(phi[p1][ix]) - 0.5 * diff * diff / comps[p1].vf) * ninv;
            }
            uni_entropy += ninv * comps[p1].vf * log(comps[p1].vf);
        }
        return vinv*(fe + entropy) - uni_entropy - uni_enthalpy;
    };

    void    update(vector<Component> & comps,
                   vector<ScalarField> & phi,
                   vector<ScalarField> & mu) {
        Vector3d    zero(0.0, 0.0, 0.0);
        double      ninv, diff;
        //update phiq
        for (int p = 0; p < M; p++) {
            phi[p].fft(phiq[p]);
        }

        //calculate order phi^2 term in free energy
        for (int p1 = 0; p1 < M; p1++) {
            
            //update Gamma * phi term
            for (int ix = 0; ix < size; ix++) {
                imag2nd[ix] += zero;
            }
            for (int p2 = 0; p2 < M;p2++){
                for (int ix = 0; ix < size; ix++) {
                    imag2nd[ix] += vertex[M*p1 + p2][ix] * phiq[p2][ix];
                }
            }
            imag2nd.fft(real2nd);
            //end update Gamma * phi term

            //update chemical potential
            ninv = 1.0 / comps[p1].N;
            for (int ix = 0; ix < size; ix++) {
                diff = phi[p1][ix] - comps[p1].vf;
                mu[p1][ix] = real2nd[ix] + ninv * (log(phi[p1][ix]) + 1.0 - diff/comps[p1].vf);
            }
        }
    };

    void updateVertex(vector<Component> & comps) {
        double  det;
        int     qx2    = nx / 2;
        int     qy2    = ny / 2;
        int     qz2    = nz / 2;
        int     center = (nx*ny*nz) / 2;

        // T.Kawakatsu, Statistical Physics of Polymers An Introduction(Springer, 2009)
        // Sec. 4.1.3 Evaluation of Expansion Coefficients Using the Random Phase Approximation

        for (int qz = 0; qz < qz2;qz++){
        for (int qy = 0; qy < qy2;qy++){
        for (int qx = 0; qx < qx2;qx++){

            //clear matrix elements
            for (int p1 = 0; p1 < (M*M); p1++) {
                matrix[p1] = 0.0;
            }

            //calculate matrix for real correlation, A_{KK'} = delta_{KK'} + S^(0)_{KK''}(chi_{K''K'} + alpha)
            for (int p1 = 0; p1 < M; p1++) {
                for (int p2 = 0; p2 < M; p2++) {
                    for (int p3 = 0; p3 < M;p3++){
                        matrix[p1*M + p2] += idealCorrelation(comps, qx, qy, qz, p1, p3) * ( chi[p3*M+p2] + alpha );
                    }
                }
                matrix[p1*(M + 1)] += 1.0;
                //delta_{KK'}, matrix[p1*M+p1], p1==p2
            }

            //get A_{KK'}^{-1}
            setInvMatrix(matrix, invmat, x, indx, det, M);
            
            //clear matrix elements
            for (int p1 = 0; p1 < ( M*M ); p1++) {
                matrix[p1] = 0.0;
            }

            //calculate real correlation, S^(real)_{KK'} = A^{-1}_{KK''}S^(0)_{K''K'}
            for (int p1 = 0; p1 < M; p1++) {
                for (int p2 = 0; p2 < M; p2++) {
                    for (int p3 = 0; p3 < M; p3++) {
                        matrix[p1*M + p2] += invmat[p1*M + p3] * idealCorrelation(comps, qx, qy, qz, p3, p2);
                    }
                }
            }

            //output(matrix, qx, qy, qz);

            //set vertex, GammaKK' = (S^(real))^{-1}_{KK'}
            setInvMatrix(matrix, invmat, x, indx, det, M);
            //output(invmat, qx, qy, qz);
            //lobj << qx << "\t" << qy << "\t" << qz << "\t" << invmat[0]/**invmat[3]*/-invmat[1]/**invmat[2]*/ << endl;
            //output_ic(comps, qx, qy, qz);

            for (int p = 0; p < M*M; p++) {
                vertex[p][center - qx - qy*nx - qz*nx*ny] = invmat[p];
                
                //copy (---) -> (--+),(-+-),(+--),(-++),(+-+),(++-),(+++)
                vertex[p][center + qx - qy*nx - qz*nx*ny] = invmat[p];
                vertex[p][center - qx + qy*nx - qz*nx*ny] = invmat[p];
                vertex[p][center - qx - qy*nx + qz*nx*ny] = invmat[p];
                vertex[p][center + qx + qy*nx - qz*nx*ny] = invmat[p];
                vertex[p][center + qx - qy*nx + qz*nx*ny] = invmat[p];
                vertex[p][center - qx + qy*nx + qz*nx*ny] = invmat[p];
                vertex[p][center + qx + qy*nx + qz*nx*ny] = invmat[p];
                //end copy
            }
            

        }//qx
        //lobjReturn();
        }//qy
        //lobjReturn();
        }//qz


        for (int p = 0; p < M*M; p++) {
            vertex[p][center] = 0.0;//origin makes inf/nan
        }
        return;
    };
};

#endif

