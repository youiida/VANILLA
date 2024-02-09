#ifndef _INC_FIELD_H_
#define _INC_FIELD_H_

using namespace std;
//#pragma warning (disable : 4290)//nothrow

#ifndef _DCOMPLEX_
#define _DCOMPLEX_
#include <complex>
typedef complex<double> dcomplex;
#endif
#ifndef PI
#define PI  ( 3.14159265358979 )
#define PI2 ( 6.28318530717959 )
#endif

//#define USE_FFTW


#ifdef USE_FFTW
#include <fftw3.h>
    #ifdef _CONSOLE //windows(defined in Project property(VisualStudio))
        #pragma comment(lib, "libfftw3-3.lib")
        #pragma comment(lib, "libfftw3f-3.lib")
        #pragma comment(lib, "libfftw3l-3.lib")
    #endif //_CONSOLE
#else
//include fft header...

#endif //USE_FFTW

#include <vector>

#ifndef _MESH_
#define _MESH_
struct Mesh {
    int       nx, ny, nz;
    double    dx, dy, dz;
};
#endif

#define DEV2 (0.5)
#define DEV4 (0.25)
#define DEV8 (0.125)


template <class T>
class TemplateField {
protected:
    vector<T>    array;
    int          nx, ny, nz, nxy;
    double       dx, dy, dz;
    double       dxinv, dyinv, dzinv;

    inline int px(const int i) const {
        if ((i + 1) % nx != 0) return i + 1;
        else return i - (nx - 1);
    };
    inline int py(const int i) const {
        if ((i % nxy) < nx * (ny - 1))return i + nx;
        else return i - nx * (ny - 1);
    };
    inline int pz(const int i) const {
        if (i < nxy * (nz - 1)) return i + nxy;
        else return i - nxy*(nz - 1);
    };
    inline int mx(const int i) const {
        if ((i) % nx != 0)    return i - 1;
        else return i + (nx - 1);
    };
    inline int my(const int i) const {
        if (i % nxy >= nx) return i - nx;
        else return i + nx*(ny - 1);
    };
    inline int mz(const int i) const {
        if (i >= nxy) return i - nxy;
        else return i + nxy*(nz - 1);
    };

    //for fft
    bool            fftflag;
#ifdef USE_FFTW
    fftw_complex    *in1, *out1;
    fftw_plan        plan1;
#endif
    inline void init_by_mesh_real(const Mesh & _mesh) {
        nx = _mesh.nx;    ny = _mesh.ny;    nz = _mesh.nz;
        dx = _mesh.dx;    dy = _mesh.dy;    dz = _mesh.dz;
        nxy = nx*ny;
        dxinv = 1.0 / dx;    dyinv = 1.0 / dy;    dzinv = 1.0 / dz;
        array.resize(nx*ny*nz);
        fftflag = false;
        return;
    };
    inline void init_by_mesh_real(const Mesh & _mesh, T & val) {
        nx = _mesh.nx;    ny = _mesh.ny;    nz = _mesh.nz;
        dx = _mesh.dx;    dy = _mesh.dy;    dz = _mesh.dz;
        nxy = nx*ny;
        dxinv = 1.0 / dx;    dyinv = 1.0 / dy;    dzinv = 1.0 / dz;
        array.resize(nx*ny*nz, val);
        fftflag = false;
        return;
    };
    inline void init_by_mesh_imag(const Mesh & _mesh) {
        nx = _mesh.nx;    ny = _mesh.ny;    nz = _mesh.nz;
        dx = PI2 / (_mesh.nx*_mesh.dx);
        dy = PI2 / (_mesh.ny*_mesh.dy);
        dz = PI2 / (_mesh.nz*_mesh.dz);
        nxy = nx*ny;
        //dxinv = 1.0 / dx;    dyinv = 1.0 / dy;    dzinv = 1.0 / dz;
        array.resize(nx*ny*nz);
        fftflag = false;
        return;
    }
    inline void init_by_mesh_imag(const Mesh & _mesh, T & val) {
        nx = _mesh.nx;    ny = _mesh.ny;    nz = _mesh.nz;
        dx = PI2 / (_mesh.nx*_mesh.dx);
        dy = PI2 / (_mesh.ny*_mesh.dy);
        dz = PI2 / (_mesh.nz*_mesh.dz);
        nxy = nx*ny;
        //dxinv = 1.0 / dx;    dyinv = 1.0 / dy;    dzinv = 1.0 / dz;
        array.resize(nx*ny*nz, val);
        fftflag = false;
        return;
    }

public:
    TemplateField() { fftflag = false; };
    TemplateField(const Mesh & _mesh) {
        fftflag = false; 
        init_by_mesh_real(_mesh);
    };
    ~TemplateField() {};

    inline int size() const { return nx*ny*nz; };

    virtual void resize(const Mesh & _mesh) { return; };
    void real_resize(const Mesh & _mesh) {
        init_by_mesh_real(_mesh);
    };
    void imag_resize(const Mesh & _mesh) {
        init_by_mesh_imag(_mesh);
    };
    inline T &          operator [] (const int i);
    inline const T &    operator [] (const int i) const;


};

template <class T>
inline T & TemplateField<T>::operator[](const int i) {
    return array[i];
}

template <class T>
inline const T &TemplateField<T>::operator[](const int i) const {
    return array[i];
}


class ScalarField;
class CScalarField;
class StagScalarField;
class VectorField;
class CVectorField;
class StressField;
//class CStressField;

#include "Vector3d.h"
#include "Vector3c.h"
#include "Stress.h"

class ScalarField : public TemplateField<double> {
public:
    ScalarField() { fftflag = false; };
    ScalarField(const Mesh & _mesh) {
        fftflag = false; 
        resize(_mesh);
    };
    ~ScalarField() {
        if (fftflag) {
#ifdef USE_FFTW
            fftw_destroy_plan(plan1);
            fftw_free(in1);
            fftw_free(out1);
#endif
        }
    };

    inline Vector3d grad(int i) {
        return Vector3d( (array[px(i)] - array[i]) * dxinv, 
                         (array[py(i)] - array[i]) * dyinv,
                         (array[pz(i)] - array[i]) * dzinv  );
    };

    //inline double div(int i) {
    //    return      (array[px(i)] - array[i]) * dxinv
    //            + (array[py(i)] - array[i]) * dyinv
    //            + (array[pz(i)] - array[i]) * dzinv;
    //};

    inline double lap(int i) {
        return    (array[px(i)] - 2.0 * array[i] + array[mx(i)]) * dx2inv
                + (array[py(i)] - 2.0 * array[i] + array[my(i)]) * dy2inv
                + (array[pz(i)] - 2.0 * array[i] + array[mz(i)]) * dz2inv;
    };

    inline Vector3d stag(int i) {
        return Vector3d((array[px(i)] + array[i]) * DEV2,
                        (array[py(i)] + array[i]) * DEV2,
                        (array[pz(i)] + array[i]) * DEV2);
    };
    Vector3d gradDiag(int i) {
        return Vector3d(DEV2 * (array[px(i)] - array[mx(i)]) * dxinv,
                        DEV2 * (array[py(i)] - array[my(i)]) * dyinv,
                        DEV2 * (array[pz(i)] - array[mz(i)]) * dzinv );
    };

    double    macValues(int i) {
        return     (array[px(i)] + array[mx(i)]) * dx2inv
                 + (array[py(i)] + array[my(i)]) * dy2inv
                 + (array[pz(i)] + array[mz(i)]) * dz2inv;
    };

    Stress toTensor(int i) {
        return Stress(array[i],
                      DEV4*(array[i] + array[px(i)] + array[py(i)] + array[px(py(i))]),
                      DEV4*(array[i] + array[px(i)] + array[pz(i)] + array[px(pz(i))]),
                      DEV4*(array[i] + array[px(i)] + array[py(i)] + array[px(py(i))]),
                      array[i],
                      DEV4*(array[i] + array[py(i)] + array[pz(i)] + array[py(pz(i))]),
                      DEV4*(array[i] + array[px(i)] + array[pz(i)] + array[px(pz(i))]),
                      DEV4*(array[i] + array[py(i)] + array[pz(i)] + array[py(pz(i))]),
                      array[i]);
    }

    void resize(const Mesh & _mesh) {
        init_by_mesh_real(_mesh);
        dx2inv = dxinv*dxinv;
        dy2inv = dyinv*dyinv;
        dz2inv = dzinv*dzinv;
    };
    void resize(const Mesh & _mesh, double val) {
        init_by_mesh_real(_mesh, val);
        dx2inv = dxinv*dxinv;
        dy2inv = dyinv*dyinv;
        dz2inv = dzinv*dzinv;
    };

//fft
#ifdef USE_FFTW
    void fftSetup(void) {
        in1   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        out1  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        plan1 = fftw_plan_dft_3d(nx, ny, nz, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
        fftflag = true;
        return;
    };
    void fft(CScalarField & cfield);
    void fftScale(CScalarField & cfield);
#else
    void fftSetup(void) { return; };
    void fft(CScalarField & cfield) { return; };
    void fftScale(CScalarField & cfield) { return; };
#endif

private:
    friend class CScalarField;
    friend class VectorField;
    double dx2inv, dy2inv, dz2inv;
};



class CScalarField : public TemplateField<dcomplex> {
public:
    CScalarField() { fftflag = false; };
    CScalarField(const Mesh & _mesh) {
        fftflag = false; 
        resize(_mesh);
    };
    ~CScalarField() {
        if (fftflag) {
#ifdef USE_FFTW
            fftw_destroy_plan(plan1);
            fftw_free(in1);
            fftw_free(out1);
#endif
        }
    };
#ifdef USE_FFTW
    void fftSetup(void) {
        in1   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        out1  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        plan1 = fftw_plan_dft_3d(nx, ny, nz, in1, out1, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftflag = true;
        return;
    };
    void fft(ScalarField & field);
    void fftScale(ScalarField & field);
#else
    void fftSetup(void) { return; };
    void fft(ScalarField & field) { return; };
    void fftScale(ScalarField & field) { return; };
#endif

    void resize(const Mesh & _mesh) {
        init_by_mesh_imag(_mesh);
    };
    void resize(const Mesh & _mesh, dcomplex & val) {
        init_by_mesh_imag(_mesh, val);
    };
private:
    friend class ScalarField;
};


class VectorField : public TemplateField<Vector3d> {
public:
    VectorField() { fftflag = false; };
    VectorField(const Mesh & _mesh) {
        fftflag = false; 
        resize(_mesh);
    };
    ~VectorField() {
        if (fftflag) {
#ifdef USE_FFTW
            fftw_destroy_plan(plan1);
            fftw_free(in1);
            fftw_free(out1);
            fftw_destroy_plan(plan2);
            fftw_free(in2);
            fftw_free(out2);
            fftw_destroy_plan(plan3);
            fftw_free(in3);
            fftw_free(out3);
#endif
        }
    };

    void resize(const Mesh & _mesh) {
        init_by_mesh_real(_mesh);
        dx2inv = dxinv*dxinv;
        dy2inv = dyinv*dyinv;
        dz2inv = dzinv*dzinv;
    };
    void resize(const Mesh & _mesh, Vector3d & val) {
        init_by_mesh_real(_mesh, val);
        dx2inv = dxinv*dxinv;
        dy2inv = dyinv*dyinv;
        dz2inv = dzinv*dzinv;
    };

#ifdef USE_FFTW
    void fftSetup(void) {
        in1   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        out1  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        plan1 = fftw_plan_dft_3d(nx, ny, nz, in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);

        in2   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        out2  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        plan2 = fftw_plan_dft_3d(nx, ny, nz, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);

        in3   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        out3  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        plan3 = fftw_plan_dft_3d(nx, ny, nz, in3, out3, FFTW_FORWARD, FFTW_ESTIMATE);
        fftflag = true;
        return;
    };
    void fft(CVectorField & cvfield);
    void fftScale(CVectorField & cvfield);
#else
    void fftSetup(void) { return; };
    void fft(CVectorField & cvfield) { return; };
    void fftScale(CVectorField & cvfield) { return; };
#endif

    inline double div(int i) {
        return    (array[i].x - array[mx(i)].x) * dxinv
                + (array[i].y - array[my(i)].y) * dyinv
                + (array[i].z - array[mz(i)].z) * dzinv;
    };
    inline double lap(int i) {
        return    (array[px(i)].x - 2.0 * array[i].x + array[mx(i)].x) * dx2inv
                + (array[py(i)].y - 2.0 * array[i].y + array[my(i)].y) * dy2inv
                + (array[pz(i)].z - 2.0 * array[i].z + array[mz(i)].z) * dz2inv;
    };

    inline Vector3d    to_diag(int i) { //vx,vy,vz -> diag point
        return Vector3d(
            DEV2 * (array[i].x + array[mx(i)].x),
            DEV2 * (array[i].y + array[my(i)].y),
            DEV2 * (array[i].z + array[mz(i)].z)
            );
    };
    inline Vector3d    to_xy(int i) {//vx,vy,vz -> xy point
        return Vector3d(
            DEV2 * (array[i].x + array[py(i)].x),
            DEV2 * (array[i].y + array[px(i)].y),
            DEV8 * (array[i].z + array[mz(i)].z + array[px(i)].z + array[px(mz(i))].z + array[py(i)].z + array[py(mz(i))].z + array[px(py(i))].z + array[px(py(mz(i)))].z)
            );
    };
    inline Vector3d    to_yz(int i) {//vx,vy,vz -> yz point
        return Vector3d(
            DEV8 * (array[i].x + array[mx(i)].x + array[py(i)].x + array[py(mx(i))].x + array[pz(i)].x + array[pz(mx(i))].x + array[py(pz(i))].x + array[py(pz(mx(i)))].x),
            DEV2 * (array[i].y + array[pz(i)].y),
            DEV2 * (array[i].z + array[py(i)].z)
            );
    };
    inline Vector3d    to_zx(int i) {//vx,vy,vz -> zx point
        return Vector3d(
            DEV2 * (array[i].x + array[pz(i)].x),
            DEV8 * (array[i].y + array[my(i)].y + array[px(i)].y + array[px(my(i))].y + array[pz(i)].y + array[pz(my(i))].y + array[pz(px(i))].y + array[pz(px(my(i)))].y),
            DEV2 * (array[i].z + array[px(i)].z)
            );
    };
    inline Stress nabla(int i) {// nabla v == del_i v_j
        return Stress(
            dxinv * (array[i].x - array[mx(i)].x),//xx
            dxinv * (array[px(i)].y - array[i].y),//xy
            dxinv * (array[px(i)].z - array[i].z),//xz
            dyinv * (array[py(i)].x - array[i].x),//yx
            dyinv * (array[i].y - array[my(i)].y),//yy
            dyinv * (array[py(i)].z - array[i].z),//yz
            dzinv * (array[pz(i)].x - array[i].x),//zx
            dzinv * (array[pz(i)].y - array[i].y),//zy
            dzinv * (array[i].z - array[mz(i)].z)//zz
            );
    };

    inline Vector3d nablavx_diag(int i) {//(del/del_x, del/del_y, del/del_z)vx at diag-point
        return Vector3d(
            dxinv * (array[i].x - array[mx(i)].x),
            dyinv * (array[py(i)].x - array[my(i)].x + array[mx(py(i))].x - array[mx(my(i))].x) * DEV4,
            dzinv * (array[pz(i)].x - array[mz(i)].x + array[mx(pz(i))].x - array[mx(mz(i))].x) * DEV4
            );
    };
    inline Vector3d nablavx_xy(int i) {
        return Vector3d(
            dxinv * (array[px(i)].x - array[mx(i)].x + array[py(px(i))].x - array[py(mx(i))].x) * DEV4,
            dyinv * (array[py(i)].x - array[i].x),
            dzinv * (array[pz(i)].x - array[mz(i)].x + array[py(pz(i))].x - array[py(mz(i))].x) * DEV4
            );
    };
    inline Vector3d nablavx_zx(int i) {
        return Vector3d(
            dxinv * (array[px(i)].x - array[mx(i)].x + array[pz(px(i))].x - array[pz(mx(i))].x) * DEV4,
            dyinv * (array[py(i)].x - array[my(i)].x + array[pz(py(i))].x - array[pz(my(i))].x) * DEV4,
            dzinv * (array[pz(i)].x - array[i].x)
            );
    };
    inline Vector3d nablavy_diag(int i) {
        return Vector3d(
            dxinv * (array[px(i)].y - array[mx(i)].y + array[my(px(i))].y - array[my(mx(i))].y) * DEV4,
            dyinv * (array[i].y - array[my(i)].y),
            dzinv * (array[pz(i)].y - array[mz(i)].y + array[my(pz(i))].y - array[my(mz(i))].y) * DEV4
            );
    };
    inline Vector3d nablavy_xy(int i) {
        return Vector3d(
            dxinv * (array[px(i)].y - array[i].y),
            dyinv * (array[py(i)].y - array[my(i)].y + array[px(py(i))].y - array[px(my(i))].y) * DEV4,
            dzinv * (array[pz(i)].y - array[mz(i)].y + array[px(pz(i))].y - array[px(mz(i))].y) * DEV4
            );
    };
    inline Vector3d nablavy_yz(int i) {
        return Vector3d(
            dxinv * (array[px(i)].y - array[mx(i)].y + array[pz(px(i))].y - array[pz(mx(i))].y) * DEV4,
            dyinv * (array[py(i)].y - array[my(i)].y + array[pz(py(i))].y - array[pz(my(i))].y) * DEV4,
            dzinv * (array[pz(i)].y - array[i].y)
            );
    };
    inline Vector3d nablavz_diag(int i) {
        return Vector3d(
            dxinv * (array[px(i)].z - array[mx(i)].x + array[mz(px(i))].z - array[mz(mx(i))].x) * DEV4,
            dyinv * (array[py(i)].z - array[my(i)].x + array[mz(py(i))].z - array[mz(my(i))].x) * DEV4,
            dzinv * (array[i].z - array[mz(i)].z)
            );
    };
    inline Vector3d nablavz_zx(int i) {
        return Vector3d(
            dxinv * (array[px(i)].z - array[i].z),
            dyinv * (array[py(i)].z - array[my(i)].z + array[px(py(i))].z - array[px(my(i))].z) * DEV4,
            dzinv * (array[pz(i)].z - array[mz(i)].z + array[px(pz(i))].z - array[px(mz(i))].z) * DEV4
            );
    };
    inline Vector3d nablavz_yz(int i) {
        return Vector3d(
            dxinv * (array[px(i)].z - array[mx(i)].z + array[py(px(i))].z - array[py(mx(i))].z) * DEV4,
            dyinv * (array[py(i)].z - array[i].z),
            dzinv * (array[pz(i)].z - array[mz(i)].z + array[py(pz(i))].z - array[py(mz(i))].z) * DEV4
            );
    };
private:
    friend class CVectorField;
#ifdef USE_FFTW
    fftw_complex    *in2, *out2, *in3, *out3;
    fftw_plan        plan2, plan3;
#endif
    double            dx2inv, dy2inv, dz2inv;
};


class CVectorField : public TemplateField< Vector3c > {
public:
    CVectorField() { fftflag = false; };
    CVectorField(const Mesh & _mesh) {
        fftflag = false; 
        resize(_mesh);
    };
    ~CVectorField() {
        if (fftflag) {
#ifdef USE_FFTW
            fftw_destroy_plan(plan1);
            fftw_free(in1);
            fftw_free(out1);
            fftw_destroy_plan(plan2);
            fftw_free(in2);
            fftw_free(out2);
            fftw_destroy_plan(plan3);
            fftw_free(in3);
            fftw_free(out3);
#endif
        }
    };

    void resize(const Mesh & _mesh) {
        init_by_mesh_imag(_mesh);
    };
    void resize(const Mesh & _mesh, Vector3c & val) {
        init_by_mesh_imag(_mesh, val);
    };

#ifdef USE_FFTW
    void fftSetup(void) {
        in1   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        out1  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        plan1 = fftw_plan_dft_3d(nx, ny, nz, in1, out1, FFTW_BACKWARD, FFTW_ESTIMATE);

        in2   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        out2  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        plan2 = fftw_plan_dft_3d(nx, ny, nz, in2, out2, FFTW_BACKWARD, FFTW_ESTIMATE);

        in3   = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        out3  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
        plan3 = fftw_plan_dft_3d(nx, ny, nz, in3, out3, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftflag = true;
        return;
    };
    void fft(VectorField & vfield);
    void fftScale(VectorField & vfield);
#else
    void fftSetup(void) { return; };
    void fft(VectorField & vfield) { return; };
    void fftScale(VectorField & vfield) { return; };
#endif

private:
    friend class    VectorField;
#ifdef USE_FFTW
    fftw_complex    *in2, *out2, *in3, *out3;
    fftw_plan        plan2, plan3;
#endif
};


class StressField : public TemplateField<Stress> {
private:
    dcomplex get_x(complex<Vector3d> & cv) {
        return dcomplex(cv.real().x, cv.imag().x);
    };
    dcomplex get_y(complex<Vector3d> & cv) {
        return dcomplex(cv.real().y, cv.imag().y);
    };
    dcomplex get_z(complex<Vector3d> & cv) {
        return dcomplex(cv.real().z, cv.imag().z);
    };

public:
    StressField() { fftflag = false; };
    StressField(const Mesh & _mesh) {
        fftflag = false; 
        resize(_mesh);
    };
    ~StressField() {
        //if (fftflag) {
        //    fftw_destroy_plan(plan1);
        //    fftw_destroy_plan(plan2);
        //    fftw_destroy_plan(plan3);
        //    fftw_destroy_plan(plan4);
        //    fftw_destroy_plan(plan5);
        //    fftw_destroy_plan(plan6);
        //    fftw_free(in1);
        //    fftw_free(in2);
        //    fftw_free(in3);
        //    fftw_free(in4);
        //    fftw_free(in5);
        //    fftw_free(in6);
        //    fftw_free(out1);
        //    fftw_free(out2);
        //    fftw_free(out3);
        //    fftw_free(out4);
        //    fftw_free(out5);
        //    fftw_free(out6);
        //}
    };

    void resize(const Mesh & _mesh) {
        init_by_mesh_real(_mesh);
    };
    void resize(const Mesh & _mesh, Stress & val) {
        init_by_mesh_real(_mesh, val);
    };

    Vector3d rowx_diag(int i) {//(s_xx s_xy s_xz) at diag
        return Vector3d(
            array[i].xx,
            DEV4 * (array[i].xy + array[my(i)].xy + array[mx(i)].xy + array[mx(my(i))].xy),
            DEV4 * (array[i].xz + array[mz(i)].xz + array[mx(i)].xz + array[mx(mz(i))].xz)
            );
    };
    Vector3d colx_diag(int i) {//(s_xx s_yx s_zx) at diag
        return Vector3d(
            array[i].xx,
            DEV4 * (array[i].yx + array[my(i)].yx + array[mx(i)].yx + array[mx(my(i))].yx),
            DEV4 * (array[i].zx + array[mz(i)].zx + array[mx(i)].zx + array[mx(mz(i))].zx)
            );
    };
    Vector3d rowx_xy(int i) {//(s_xx s_xy s_xz) at xy
        return Vector3d(
            DEV4 * (array[i].xx + array[px(i)].xx + array[py(i)].xx + array[py(px(i))].xx),
            array[i].xy,
            DEV4 * (array[i].xz + array[mz(i)].xz + array[py(i)].xz + array[py(mz(i))].xz)
            );
    };
    Vector3d colx_xy(int i) {//(s_xx s_yx s_zx) at xy
        return Vector3d(
            DEV4 * (array[i].xx + array[px(i)].xx + array[py(i)].xx + array[py(px(i))].xx),
            array[i].yx,
            DEV4 * (array[i].zx + array[mz(i)].zx + array[py(i)].zx + array[py(mz(i))].zx)
            );
    };
    Vector3d rowx_zx(int i) {//(s_xx s_xy s_xz) at zx
        return Vector3d(
            DEV4 * (array[i].xx + array[px(i)].xx + array[pz(i)].xx + array[pz(px(i))].xx),
            DEV4 * (array[i].xy + array[my(i)].xy + array[pz(i)].xy + array[pz(my(i))].xy),
            array[i].xz
            );
    };
    Vector3d colx_zx(int i) {//(s_xx s_yx s_zx) at zx
        return Vector3d(
            DEV4 * (array[i].xx + array[px(i)].xx + array[pz(i)].xx + array[pz(px(i))].xx),
            DEV4 * (array[i].yx + array[my(i)].yx + array[pz(i)].yx + array[pz(my(i))].yx),
            array[i].zx
            );
    };
    Vector3d rowy_diag(int i) {//(s_yx s_yy s_yz) at diag
        return Vector3d(
            DEV4 * (array[i].yx + array[mx(i)].yx + array[my(i)].yx + array[my(mx(i))].yx),
            array[i].yy,
            DEV4 * (array[i].yz + array[mz(i)].yz + array[my(i)].yz + array[my(mz(i))].yz)
            );
    };
    Vector3d coly_diag(int i) {//(s_xy s_yy s_zy) at diag
        return Vector3d(
            DEV4 * (array[i].xy + array[mx(i)].xy + array[my(i)].xy + array[my(mx(i))].xy),
            array[i].yy,
            DEV4 * (array[i].zy + array[mz(i)].zy + array[my(i)].zy + array[my(mz(i))].zy)
            );
    };
    Vector3d rowy_xy(int i) {//(s_yx s_yy s_yz) at xy
        return Vector3d(
            array[i].yx,
            DEV4 * (array[i].yy + array[py(i)].yy + array[px(i)].yy + array[px(py(i))].yy),
            DEV4 * (array[i].yz + array[mz(i)].yz + array[px(i)].yz + array[px(mz(i))].yz)
            );
    };
    Vector3d coly_xy(int i) {//(s_xy s_yy s_zy) at xy
        return Vector3d(
            array[i].xy,
            DEV4 * (array[i].yy + array[py(i)].yy + array[px(i)].yy + array[px(py(i))].yy),
            DEV4 * (array[i].zy + array[mz(i)].zy + array[px(i)].zy + array[px(mz(i))].zy)
            );
    };
    Vector3d rowy_yz(int i) {//(s_yx s_yy s_yz) at yz
        return Vector3d(
            DEV4 * (array[i].yx + array[mx(i)].yx + array[pz(i)].yx + array[pz(mx(i))].yx),
            DEV4 * (array[i].yy + array[py(i)].yy + array[pz(i)].yy + array[pz(py(i))].yy),
            array[i].yz
            );
    };
    Vector3d coly_yz(int i) {//(s_xy s_yy s_zy) at zy
        return Vector3d(
            DEV4 * (array[i].xy + array[mx(i)].xy + array[pz(i)].xy + array[pz(mx(i))].xy),
            DEV4 * (array[i].yy + array[py(i)].yy + array[pz(i)].yy + array[pz(py(i))].yy),
            array[i].zy
            );
    };
    Vector3d rowz_diag(int i) {//(s_zx s_zy s_zz) at diag
        return Vector3d(
            DEV4 * (array[i].zx + array[mx(i)].zx + array[mz(i)].zx + array[mz(mx(i))].zx),
            DEV4 * (array[i].zy + array[my(i)].zy + array[mz(i)].zy + array[mz(my(i))].zy),
            array[i].zz
            );
    };
    Vector3d colz_diag(int i) {//(s_xz s_yz s_zz) at diag
        return Vector3d(
            DEV4 * (array[i].xz + array[mx(i)].xz + array[mz(i)].xz + array[mz(mx(i))].xz),
            DEV4 * (array[i].yz + array[my(i)].yz + array[mz(i)].yz + array[mz(my(i))].yz),
            array[i].zz
            );
    };
    Vector3d rowz_zx(int i) {//(s_zx s_zy s_zz) at zx
        return Vector3d(
            array[i].zx,
            DEV4 * (array[i].zy + array[my(i)].zy + array[px(i)].zy + array[px(my(i))].zy),
            DEV4 * (array[i].zz + array[pz(i)].zz + array[px(i)].zz + array[px(pz(i))].zz)
            );
    };
    Vector3d colz_zx(int i) {//(s_xz s_yz s_zz) at xz
        return Vector3d(
            array[i].xz,
            DEV4 * (array[i].yz + array[my(i)].yz + array[px(i)].yz + array[px(my(i))].yz),
            DEV4 * (array[i].zz + array[pz(i)].zz + array[px(i)].zz + array[px(pz(i))].zz)
            );
    };
    Vector3d rowz_yz(int i) {//(s_zx s_zy s_zz) at yz
        return Vector3d(
            DEV4 * (array[i].zx + array[mx(i)].zx + array[py(i)].zx + array[py(mx(i))].zx),
            array[i].zy,
            DEV4 * (array[i].zz + array[pz(i)].zz + array[py(i)].zz + array[py(pz(i))].zz)
            );
    };
    Vector3d colz_yz(int i) {//(s_xz s_yz s_zz) at yz
        return Vector3d(
            DEV4 * (array[i].xz + array[mx(i)].xz + array[py(i)].xz + array[py(mx(i))].xz),
            array[i].yz,
            DEV4 * (array[i].zz + array[pz(i)].zz + array[py(i)].zz + array[py(pz(i))].zz)
            );
    };

    //check 20160810 youiida
    Vector3d nabladot(int i) {//operate (del/del_x, del/del_y, del/del_z) to (s_xx s_yx s_zx) at vx,vy,vz
        return Vector3d(
            dxinv * (array[px(i)].xx - array[i].xx) + dyinv * (array[i].yx - array[my(i)].yx) + dzinv * (array[i].zx - array[mz(i)].zx),
            dxinv * (array[i].xy - array[mx(i)].xy) + dyinv * (array[py(i)].yy - array[i].yy) + dzinv * (array[i].zy - array[mz(i)].zy),
            dxinv * (array[i].xz - array[mx(i)].xz) + dyinv * (array[i].yz - array[my(i)].yz) + dzinv * (array[pz(i)].zz - array[i].zz)
            );
    };

    Vector3d grad_pm_xx(int i) {// (s_xx[px]-s_xx[mx])/2*dx , (s_xx[py]-s_xx[my])/2*dy, (s_xx[pz]-s_xx[mz])/2*dz
        return Vector3d(
            DEV2 * dxinv * (array[px(i)].xx - array[mx(i)].xx),
            DEV2 * dyinv * (array[py(i)].xx - array[my(i)].xx),
            DEV2 * dzinv * (array[pz(i)].xx - array[mz(i)].xx)
            );
    };
    Vector3d grad_pm_xy(int i) {
        return Vector3d(
            DEV2 * dxinv * (array[px(i)].xy - array[mx(i)].xy),
            DEV2 * dyinv * (array[py(i)].xy - array[my(i)].xy),
            DEV2 * dzinv * (array[pz(i)].xy - array[mz(i)].xy)
            );
    };
    Vector3d grad_pm_xz(int i) {
        return Vector3d(
            DEV2 * dxinv*(array[px(i)].xz - array[mx(i)].xz),
            DEV2 * dyinv*(array[py(i)].xz - array[my(i)].xz),
            DEV2 * dzinv*(array[pz(i)].xz - array[mz(i)].xz)
            );
    };
    Vector3d grad_pm_yx(int i) {
        return Vector3d(
            DEV2 * dxinv * (array[px(i)].yx - array[mx(i)].yx),
            DEV2 * dyinv * (array[py(i)].yx - array[my(i)].yx),
            DEV2 * dzinv * (array[pz(i)].yx - array[mz(i)].yx)
            );
    };
    Vector3d grad_pm_yy(int i) {
        return Vector3d(
            DEV2 * dxinv * (array[px(i)].yy - array[mx(i)].yy),
            DEV2 * dyinv * (array[py(i)].yy - array[my(i)].yy),
            DEV2 * dzinv * (array[pz(i)].yy - array[mz(i)].yy)
            );
    };
    Vector3d grad_pm_yz(int i) {
        return Vector3d(
            DEV2 * dxinv * (array[px(i)].yz - array[mx(i)].yz),
            DEV2 * dyinv * (array[py(i)].yz - array[my(i)].yz),
            DEV2 * dzinv * (array[pz(i)].yz - array[mz(i)].yz)
            );
    };
    Vector3d grad_pm_zx(int i) {
        return Vector3d(
            DEV2 * dxinv * (array[px(i)].zx - array[mx(i)].zx),
            DEV2 * dyinv * (array[py(i)].zx - array[my(i)].zx),
            DEV2 * dzinv * (array[pz(i)].zx - array[mz(i)].zx)
            );
    };
    Vector3d grad_pm_zy(int i) {
        return Vector3d(
            DEV2 * dxinv * (array[px(i)].zy - array[mx(i)].zy),
            DEV2 * dyinv * (array[py(i)].zy - array[my(i)].zy),
            DEV2 * dzinv * (array[pz(i)].zy - array[mz(i)].zy)
            );
    };
    Vector3d grad_pm_zz(int i) {
        return Vector3d(
            DEV2 * dxinv * (array[px(i)].zz - array[mx(i)].zz),
            DEV2 * dyinv * (array[py(i)].zz - array[my(i)].zz),
            DEV2 * dzinv * (array[pz(i)].zz - array[mz(i)].zz)
            );
    };


    //void fftSetup(void) {
    //    in1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    in3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    in4 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    in5 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    in6 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    out1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    out3 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    out4 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    out5 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    out6 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * size());
    //    plan1 = fftw_plan_dft_3d(nx, ny, nz, in,  out,  FFTW_FORWARD, FFTW_ESTIMATE);
    //    plan2 = fftw_plan_dft_3d(nx, ny, nz, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
    //    plan3 = fftw_plan_dft_3d(nx, ny, nz, in3, out3, FFTW_FORWARD, FFTW_ESTIMATE);
    //    plan4 = fftw_plan_dft_3d(nx, ny, nz, in4, out4, FFTW_FORWARD, FFTW_ESTIMATE);
    //    plan5 = fftw_plan_dft_3d(nx, ny, nz, in5, out5, FFTW_FORWARD, FFTW_ESTIMATE);
    //    plan6 = fftw_plan_dft_3d(nx, ny, nz, in6, out6, FFTW_FORWARD, FFTW_ESTIMATE);
    //    fftflag = true;
    //    return;
    //};
    ////void fft(CStressField & cs);
    //void fft_scale(CStressField & cs);

private:
    //fftw_complex    *in2, *in3, *in4, *in5, *in6, *out2, *out3, *out4, *out5, *out6;
    //fftw_plan        plan2, plan3, plan4, plan5, plan6;
};

///////////////////////////////////////////////////////////////////////////////////////////
//
//    field operations
//
///////////////////////////////////////////////////////////////////////////////////////////

//check 20160810 youiida
Vector3d nablaDotStrss(StressField & st, ScalarField & sc, int i);
double   vectorDotNablaScalar(VectorField & vf, ScalarField & sf, int i);
Stress   vectorDotNablaStress(VectorField & vf, StressField & sf, int i);
Stress   stressDotNablaVectorSymmetry(StressField & sf, VectorField & vf, int i);
#endif
