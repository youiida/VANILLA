#ifndef _INC_STRESS_H_
#define _INC_STRESS_H_


//#include <math.h>    // for sqrt()
#include <iostream>    // for << and >>#include <iostream>
#include "Vector3d.h"
#pragma warning (disable : 4290)//nothrow
#include "udfstream.h"

#if !defined(_WITHOUT_USING_NAMESPACE_STD)
using namespace std;
#endif

class Stress{
private:
public:
    //Stress(){};
    //~Stress(){};
    double    xx, xy, xz, yx, yy, yz, zx, zy, zz;
    //----- Constructors -----
    Stress(){};
    Stress(double _xx, double _xy, double _xz, 
           double _yx, double _yy, double _yz,
           double _zx, double _zy, double _zz)
                    : xx(_xx), xy(_xy), xz(_xz), yx(_yx), yy(_yy), yz(_yz), zx(_zx), zy(_zy), zz(_zz){}
    //Stress(const Stress& _in): xx(_in.xx), xy(_in.xy), xz(_in.xz), yx(_in.yx), yy(_in.yy), yz(_in.yz), zx(_in.zx), zy(_in.zy), zz(_in.zz) {}
    //~Stress(){cout << "delete: stress" << endl;};
    //----- get/set values -----

    void set(double _xx, double _xy, double _xz, 
             double _yx, double _yy, double _yz, 
             double _zx, double _zy, double _zz){
        xx = _xx;    xy = _xy;    xz = _xz;
        yx = _yx;    yy = _yy;    yz = _yz;
        zx = _zx;    zy = _zy;    zz = _zz;
        return;
    };


    //----- operators -----

    Stress& operator = (const Stress& t);

    friend Stress operator + (const Stress& t1, const Stress& t2);
    friend Stress operator - (const Stress& t1, const Stress& t2);
    //friend Stress operator * (const Stress& t1, const Stress& t2);
    Stress& operator += (const Stress& t);
    Stress& operator -= (const Stress& t);

    friend Stress operator * (const Stress& t, const double& s);
    friend Stress operator * (const double& s, const Stress &t);
    //friend Stress operator / (const Stress& t, const double& s);
    //friend Stress operator * (const Stress& t, const int& s);
    //friend Stress operator * (const int& s, const Stress &t);
    //friend Stress operator / (const Stress& t, const int& s);
    Stress& operator *= (const double& s);
    Stress& operator /= (const double& s);

    //friend Stress operator * (const Stress& t1, const Stress& t2);

    //----- functions -----
    Stress trans(void) {
        return Stress(xx, yx, zx, xy, yy, zy, xz, yz, zz);
    };
    double trace(void) { 
        return xx + yy + zz; 
    };
    Stress square(void) {
        return Stress(xx*xx, xy*xy, xz*xz, yx*yx, yy*yy, yz*yz, zx*zx, zy*zy, zz*zz);
    };
    Stress cube(void) {
        return Stress(xx*xx*xx, xy*xy*xy, xz*xz*xz, yx*yx*yx, yy*yy*yy, yz*yz*yz, zx*zx*zx, zy*zy*zy, zz*zz*zz);
    };


//#ifdef _INC_VECTOR3d_H_
//    friend Vector3d operator * (const Stress & t, const Vector3d & v);
//    friend Vector3d operator * (const Vector3d & v, const Stress & t);
//#endif

    Stress operator - () const;

    friend bool operator == (const Stress& t1, const Stress& t2);

    //-----  set/reset/get outputParen  -----

    static void parenOn() { outputParen = true; }
    static void parenOff() { outputParen = false; }
    static bool parenIsOn() { return outputParen; }
    static void xmlOn() { outputXML = true; }
    static void xmlOff() { outputXML = false; }
    static bool xmlIsOn() { return outputXML; }

    //----- set/get output precision -----

    static std::streamsize setPrecision(std::streamsize prec) {
        std::streamsize p = precision;
        precision = prec;
        return p;
    }
    static std::streamsize getPrecision(void) { return precision; }

    static void setReducedPrecision(std::streamsize prec) {
        reducedPrecision = prec;
    }
    static std::streamsize getReducedPrecision(void) { return reducedPrecision; }

private:

    static bool outputParen;    // default value is true (Vector3d.cpp)
    static bool outputXML;    // default value is false (Vector3d.cpp)
    static std::streamsize precision;    // default UDF_DOUBLE_PRECISION (14)
    static std::streamsize reducedPrecision;    // reduced precision used for Position, Velocity and Force in Struct
};

ostream& operator << (ostream& os, const Stress & t);
istream& operator >> (istream& is, Stress& t);
UDFistream& operator >> (UDFistream& is, Stress& t);


inline Stress& Stress::operator = (const Stress & t){
    xx = t.xx;    xy = t.xy;    xz = t.xz;
    yx = t.yx;    yy = t.yy;    yz = t.yz;
    zx = t.zx;    zy = t.zy;    zz = t.zz;
    return *this;
}

inline Stress operator + (const Stress & t1, const Stress & t2){
    return Stress(t1.xx+t2.xx, t1.xy+t2.xy, t1.xz+t2.xz,
                  t1.yx+t2.yx, t1.yy+t2.yy, t1.yz+t2.yz, 
                  t1.zx+t2.zx, t1.zy+t2.zy, t1.zz+t2.zz);
}

inline Stress operator - (const Stress & t1, const Stress & t2){
    return Stress(t1.xx - t2.xx, t1.xy - t2.xy, t1.xz - t2.xz,
                  t1.yx - t2.yx, t1.yy - t2.yy, t1.yz - t2.yz,
                  t1.zx - t2.zx, t1.zy - t2.zy, t1.zz - t2.zz);
}

inline Stress & Stress::operator += (const Stress& t){
    xx += t.xx;    xy += t.xy;    xz += t.xz;
    yx += t.yx;    yy += t.yy;    yz += t.yz;
    zx += t.zx; zy += t.zy;    zz += t.zz;
    return *this;
}

inline Stress & Stress::operator -= (const Stress & t){
    xx -= t.xx;    xy -= t.xy;    xz -= t.xz;
    yx -= t.yx;    yy -= t.yy;    yz -= t.yz;
    zx -= t.zx; zy -= t.zy;    zz -= t.zz;
    return *this;
}

inline Stress operator * (const Stress & t, const double& s){
    return Stress(s*t.xx, s*t.xy, s*t.xz,
                  s*t.yx, s*t.yy, s*t.yz,
                  s*t.zx, s*t.zy, s*t.zz);
}

inline Stress operator * (const double& s, const Stress &t){
    return Stress(s*t.xx, s*t.xy, s*t.xz,
                  s*t.yx, s*t.yy, s*t.yz,
                  s*t.zx, s*t.zy, s*t.zz);
}

inline Stress operator / (const Stress & t, const double& s){
    return Stress(t.xx/s, t.xy/s, t.xz/s,
                  t.yx/s, t.yy/s, t.yz/s,
                  t.zx/s, t.zy/s, t.zz/s);
}

inline Stress & Stress::operator *= (const double& s){
    xx *= s;    xy *= s;    xz *= s;
    yx *= s;    yy *= s;    yz *= s;
    zx *= s;    zy *= s;    zz *= s;
    return *this;
}

inline Stress& Stress::operator /= (const double& s){
    xx /= s;    xy /= s;    xz /= s;
    yx /= s;    yy /= s;    yz /= s;
    zx /= s;    zy /= s;    zz /= s;
    return *this;
}

inline Stress Stress::operator - () const{
    return Stress(-xx, -xy, -xz, -yx, -yy, -yz, -zx, -zy, -zz);
}

inline bool operator == (const Stress & t1, const Stress & t2){
    return ( t1.xx == t2.xx && t1.xy == t2.xy && t1.xz == t2.xz &&
             t1.yx == t2.yx && t1.yy == t2.yy && t1.yz == t2.yz &&
             t1.zx == t2.zx && t1.zy == t2.zy && t1.zz == t2.zz );
}

#endif