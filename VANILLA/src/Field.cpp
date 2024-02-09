#include "Field.h"

#ifdef USE_FFTW
void ScalarField::fft(CScalarField &cs) {
    for (int i = 0; i < size(); i++) {
        in1[i][0] = array[i];
        in1[i][1] = 0.0e0;
    }

    fftw_execute(plan1);
    for (int i = 0; i<size(); i++) {
        cs[i] = dcomplex(out1[i][0], out1[i][1]);
    }
    return;
}

void ScalarField::fftScale(CScalarField &cs) {
    for (int i = 0; i < size(); i++) {
        in1[i][0] = array[i];
        in1[i][1] = 0.0e0;
    }

    fftw_execute(plan1);
    double scale = 1.0 / double(size());
    for (int i = 0; i<size(); i++) {
        cs[i] = dcomplex(out1[i][0], out1[i][1]) * scale;
    }
    return;
}

void CScalarField::fft(ScalarField &s) {
    for (int i = 0; i < size(); i++) {
        in1[i][0] = array[i].real();
        in1[i][1] = array[i].imag();
    }

    fftw_execute(plan1);
    for (int i = 0; i < size(); i++) {
        s[i] = out1[i][0];
    }
    return;
}
void CScalarField::fftScale(ScalarField &s) {
    for (int i = 0; i < size(); i++) {
        in1[i][0] = array[i].real();
        in1[i][1] = array[i].imag();
    }

    fftw_execute(plan1);
    double scale = 1.0 / double(size());
    for (int i = 0; i < size(); i++) {
        s[i] = out1[i][0]*scale;
    }
    return;
}

void VectorField::fft(CVectorField &cv) {
    for (int i = 0; i < size(); i++) {
        in1[i][0]  = array[i].x;
        in1[i][1]  = 0.0e0;
        in2[i][0] = array[i].y;
        in2[i][1] = 0.0e0;
        in3[i][0] = array[i].z;
        in3[i][1] = 0.0e0;
    }

    fftw_execute(plan1);
    fftw_execute(plan2);
    fftw_execute(plan3);
    
    for (int i = 0; i<size(); i++) {
        cv[i].x = dcomplex(out1[i][0], out1[i][1]);
        cv[i].y = dcomplex(out2[i][0], out2[i][1]);
        cv[i].z = dcomplex(out3[i][0], out3[i][1]);
    }
    return;
}

void VectorField::fftScale(CVectorField &cv) {
    for (int i = 0; i < size(); i++) {
        in1[i][0] = array[i].x;
        in1[i][1] = 0.0e0;
        in2[i][0] = array[i].y;
        in2[i][1] = 0.0e0;
        in3[i][0] = array[i].z;
        in3[i][1] = 0.0e0;
    }

    fftw_execute(plan1);
    fftw_execute(plan2);
    fftw_execute(plan3);

    double scale = 1.0 / double(size());

    for (int i = 0; i<size(); i++) {
        cv[i].x = dcomplex(out1[i][0], out1[i][1]) * scale;
        cv[i].y = dcomplex(out2[i][0], out2[i][1]) * scale;
        cv[i].z = dcomplex(out3[i][0], out3[i][1]) * scale;
    }
    return;
}

void CVectorField::fft(VectorField &v) {
    for (int i = 0; i < size(); i++) {
        in1[i][0] = array[i].x.real();
        in1[i][1] = array[i].x.imag();
        in2[i][0] = array[i].y.real();
        in2[i][1] = array[i].y.imag();
        in3[i][0] = array[i].z.real();
        in3[i][1] = array[i].z.imag();
    }

    fftw_execute(plan1);
    fftw_execute(plan2);
    fftw_execute(plan3);
    
    for (int i = 0; i<size(); i++) {
        v[i] = Vector3d(out1[i][0], out2[i][0], out3[i][0]);
    }
    return;
}

void CVectorField::fftScale(VectorField &v) {
    for (int i = 0; i < size(); i++) {
        in1[i][0] = array[i].x.real();
        in1[i][1] = array[i].x.imag();
        in2[i][0] = array[i].y.real();
        in2[i][1] = array[i].y.imag();
        in3[i][0] = array[i].z.real();
        in3[i][1] = array[i].z.imag();
    }

    fftw_execute(plan1);
    fftw_execute(plan2);
    fftw_execute(plan3);
    double scale = 1.0 / double(size());

    for (int i = 0; i<size(); i++) {
        v[i] = Vector3d(out1[i][0] * scale, out2[i][0] * scale, out3[i][0] * scale);
    }
    return;
}

//void StressField::fft_scale(CStressField & cs) {
//    for (int i = 0; i < size(); i++) {
//        in[i][0] = array[i].xx;
//        in[i][1] = 0.0e0;
//        in2[i][0] = array[i].xy;
//        in2[i][1] = 0.0e0;
//        in3[i][0] = array[i].xz;
//        in3[i][1] = 0.0e0;
//        in4[i][0] = array[i].yy;
//        in4[i][1] = 0.0e0;
//        in5[i][0] = array[i].yz;
//        in5[i][1] = 0.0e0;
//        in6[i][0] = array[i].zz;
//        in6[i][1] = 0.0e0;
//    }
//    
//    fftw_execute(plan);
//    fftw_execute(plan2);
//    fftw_execute(plan3);
//    fftw_execute(plan4);
//    fftw_execute(plan5);
//    fftw_execute(plan6);
//
//    double scale = 1.0 / double(size());
//
//    for (int i = 0; i < size(); i++) {
//        cs[i].xx = scale * dcomplex(out[i][0], out[i][1]);
//        cs[i].xy = scale * dcomplex(out2[i][0], out2[i][1]);
//        cs[i].xz = scale * dcomplex(out3[i][0], out3[i][1]);
//        cs[i].yy = scale * dcomplex(out4[i][0], out4[i][1]);
//        cs[i].yz = scale * dcomplex(out5[i][0], out5[i][1]);
//        cs[i].zz = scale * dcomplex(out6[i][0], out6[i][1]);
//
//        cs[i].yx = cs[i].xy;
//        cs[i].zx = cs[i].zx;
//        cs[i].zy = cs[i].yz;
//    }
//    return;
//}
#endif

///////////////////////////////////////////////////////////////////////////////////////////
//
//    field operations
//
///////////////////////////////////////////////////////////////////////////////////////////

//check 20160810 youiida
Vector3d    nablaDotStrss(StressField & st,
                          ScalarField & sc,
                          int i) {
    return st.nabladot(i) + sc.grad(i);
};

double    vectorDotNablaScalar(VectorField & vf,
                               ScalarField & sf,
                               int i) {
    return vf.to_diag(i) * sf.gradDiag(i);
};

Stress    vectorDotNablaStress(VectorField & vf,
                               StressField & sf,
                               int i) {
    Vector3d    v_diag = vf.to_diag(i);
    Vector3d    v_xy = vf.to_xy(i);
    Vector3d    v_yz = vf.to_yz(i);
    Vector3d    v_zx = vf.to_zx(i);
    return Stress(
        v_diag * sf.grad_pm_xx(i),
        v_xy   * sf.grad_pm_xy(i),
        v_zx   * sf.grad_pm_xz(i),
        v_xy   * sf.grad_pm_yx(i),
        v_diag * sf.grad_pm_yy(i),
        v_yz   * sf.grad_pm_yz(i),
        v_zx   * sf.grad_pm_zx(i),
        v_yz   * sf.grad_pm_zy(i),
        v_diag * sf.grad_pm_zz(i)
    );
};

//stress dot (nabla v) + (nabla v)^t dot stress
Stress    stressDotNablaVectorSymmetry(StressField & sf,
                                       VectorField & vf,
                                       int i) {
    return Stress(
        (sf.rowx_diag(i) + sf.colx_diag(i)) * vf.nablavx_diag(i),
        (sf.rowx_xy(i) + sf.coly_xy(i)) * vf.nablavx_xy(i),
        (sf.rowx_zx(i) + sf.colz_zx(i)) * vf.nablavx_zx(i),
        (sf.rowy_xy(i) + sf.colx_xy(i)) * vf.nablavy_xy(i),
        (sf.rowy_diag(i) + sf.coly_diag(i)) * vf.nablavy_diag(i),
        (sf.rowy_yz(i) + sf.colz_yz(i)) * vf.nablavy_yz(i),
        (sf.rowz_zx(i) + sf.colx_zx(i)) * vf.nablavz_zx(i),
        (sf.rowz_yz(i) + sf.coly_yz(i)) * vf.nablavz_yz(i),
        (sf.rowz_diag(i) + sf.colz_diag(i)) * vf.nablavz_diag(i)
    );//miss type * and +, bug-fix 20160613 youiida
};

