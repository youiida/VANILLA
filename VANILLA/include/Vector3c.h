//===== Vector3c.h ===========================================================
//
//        3D vector class
//
//
//===========================================================================

#ifndef _INC_VECTOR3C_H_
#define _INC_VECTOR3C_H_

#include <math.h>    // for sqrt()
#include <iostream>    // for << and >>
//ignore the setting in configenv.h and use Vector3d::precision instead
//#include "configenv.h"  // for output precision of double value
#include "udfstream.h"
#include "Vector3d.h"
#include "Stress.h"

#if !defined(_WITHOUT_USING_NAMESPACE_STD)
using namespace std;
#endif

//============================================================================
//
//  class Vector3c
//
//    data members:
//        double x, y, z;            // public
//    operators:
//        vec = vec,    - vec
//        vec + vec,    vec - vec,    vec += vec,    vec -= vec
//        double * vec,    vec * double,    vec / double
//                vec *= double,    vec /= double
//        vec * vec    (inner product)
//        (vec^vec)    (outer product)
//
//        NOTE: The operator^ has very low precedence, so
//              parentheses must be used almost always.
//
//    member functions:
//        double length()        = sqrt( vec * vec )
//        double length2()    = vec * vec
//
//
// 'Vector3d.cpp' is necessary only if the followings are used:
//
//    output to / input from stream
//        ostream& operator << (ostream& os, const Vector3d& v);
//        istream& operator >> (istream& is, Vector3d& v);
//
//    static member functions which control output format:
//        static void parenOn()    { outputParen = true; }
//        static void parenOff()    { outputParen = false; }
//        static bool parenIsOn()    { return outputParen; }
//
//    'cout << vec' writes
//        '( x, y, z )'    if parenIsOn()    ( default )
//        'x y z'        otherwise
//
//    'cin >> vec' can read either '(x, y, z)' or 'x y z'
//
//============================================================================

#include <complex>
#ifndef _DCOMPLEX_
#define _DCOMPLEX_
typedef complex<double> dcomplex;
#endif
class Vector3c {

public:
    dcomplex x, y, z;            // data members are public

                            //----- Constructors -----

    Vector3c(dcomplex vx = (0.0,0.0), dcomplex vy = (0.0, 0.0), dcomplex vz = (0.0, 0.0))
        : x(vx), y(vy), z(vz) {}

    Vector3c(const Vector3c& v) : x(v.x), y(v.y), z(v.z) {}

    //----- get/set values -----

    void set(dcomplex vx, dcomplex vy, dcomplex vz) { x = vx; y = vy; z = vz; }

    //----- operators -----

    Vector3c& operator = (const Vector3c& v);
    //    Vector3d& operator = (const double& v);

    friend Vector3c operator + (const Vector3c& v1, const Vector3c& v2);
    friend Vector3c operator - (const Vector3c& v1, const Vector3c& v2);
    Vector3c& operator += (const Vector3c& v);
    Vector3c& operator -= (const Vector3c& v);

    friend Vector3c operator * (const Vector3c& v, const double& s);
    friend Vector3c operator * (const double& s, const Vector3c &v);
    friend Vector3c operator / (const Vector3c& v, const double& s);
    friend Vector3c operator * (const Vector3c& v, const int& s);
    friend Vector3c operator * (const int& s, const Vector3c &v);
    friend Vector3c operator / (const Vector3c& v, const int& s);
    Vector3c& operator *= (const double& s);
    Vector3c& operator /= (const double& s);

    friend dcomplex operator * (const Vector3c& v1, const Vector3c& v2);
    friend dcomplex operator * (const Vector3c& v1, const Vector3d& v2);
    friend dcomplex operator * (const Vector3d& v1, const Vector3c& v2);
    friend Vector3c operator ^ (const Vector3c& v1, const Vector3c& v2);

    Vector3c operator - () const;

    friend Vector3c operator * (const Stress& s, const Vector3c& v);// vector_i = stress_ij * vector_j
    friend Vector3c operator * (const Vector3c& v, const Stress& s);// vector_j = vector_i * stress_ij

    // add for platform 2001/1/26 T.Aoyagi
    //friend bool operator == (const Vector3c& v1, const Vector3c& v2);
    //friend bool operator < (const Vector3c& v1, const Vector3c& v2);

    //-----  length of vector  -----

    double length() const { return sqrt((x*conj(x) + y*conj(y) + z*conj(z)).real()); }
    double length2() const { return (x*conj(x) + y*conj(y) + z*conj(z)).real(); }

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

    // finally modified for datatable. 2000/7/30 T.Aoyagi
    //Vector3c square(void) const { return Vector3c(x*x, y*y, z*z); }
    //Vector3c root(void) const { return Vector3c(sqrt(x), sqrt(y), sqrt(z)); }
    void clear(void) { x = y = z = 0.0; }
    // add 2003/10/28 T.Aoyagi
    //Vector3c cube(void) const { return Vector3c(x*x*x, y*y*y, z*z*z); }
    Vector3c mp(const Vector3c& v) { return Vector3c(x*v.x, y*v.y, z*v.z); }

    //youiida 20151211
    Vector3c inv() { return Vector3c(1.0 / x, 1.0 / y, 1.0 / z); }
    Vector3c mp(const Vector3c& v1, const Vector3c & v2) {
        return Vector3c(x * v1.x * v2.x,
            y * v1.y * v2.y,
            z * v1.z * v2.z);
    }
    //Vector3c mp(const Vector3c& v1, const Vector3c & v2, const Vector3c & v3) {
    //    return Vector3c(x * v1.x * v2.x * v3.x,
    //        y * v1.y * v2.y * v3.y,
    //        z * v1.z * v2.z * v3.z);
    //}
    //Vector3c mp(const Vector3c& v1, const Vector3c & v2, const Vector3c & v3, const Vector3c & v4) {
    //    return Vector3c(x * v1.x * v2.x * v3.x * v4.x,
    //        y * v1.y * v2.y * v3.y * v4.y,
    //        z * v1.z * v2.z * v3.z * v4.z);
    //}
    Vector3c dv(const Vector3c & v) { return Vector3c(x / v.x, y / v.y, z / v.z); };


    //operator double() const { return length(); }
private:

    static bool outputParen;    // default value is true (Vector3c.cpp)
    static bool outputXML;    // default value is false (Vector3c.cpp)
    static std::streamsize precision;    // default UDF_DOUBLE_PRECISION (14)
    static std::streamsize reducedPrecision;    // reduced precision used for Position, Velocity and Force in Struct
};


//============================================================================
//
//    Input / Output Vector3c objects
//
//        these two functions are defined in Vector3c.cpp
//
//============================================================================


//ostream& operator << (ostream& os, const Vector3c& v);
//istream& operator >> (istream& is, Vector3c& v);
//UDFistream& operator >> (UDFistream& is, Vector3c& v);

//============================================================================
//
//    Definition of inline operators
//
//============================================================================


//=====  Assignment  =====

inline Vector3c& Vector3c::operator = (const Vector3c& v)
{
    x = v.x; y = v.y; z = v.z;
    return *this;
}

//inline Vector3c& Vector3c::operator = (const double& d)
//{
//    x = y = z = d;
//    return *this;
//}

//=====  Vec + Vec , Vec - Vec  =====

inline Vector3c operator + (const Vector3c& v1, const Vector3c& v2)
{
    return Vector3c(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector3c operator - (const Vector3c& v1, const Vector3c& v2)
{
    return Vector3c(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector3c& Vector3c::operator += (const Vector3c& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

inline Vector3c& Vector3c::operator -= (const Vector3c& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

//=====  scalar*Vec, Vec*scalar, Vec/scalar  =====

inline Vector3c operator * (const Vector3c& v, const double& s)
{
    return Vector3c(s*v.x, s*v.y, s*v.z);
}

inline Vector3c operator * (const double& s, const Vector3c& v)
{
    return Vector3c(s*v.x, s*v.y, s*v.z);
}

inline Vector3c operator / (const Vector3c& v, const double& s)
{
    return Vector3c(v.x / s, v.y / s, v.z / s);
}

inline Vector3c operator * (const Vector3c& v, const int& s)
{
    return Vector3c((double)s*v.x, (double)s*v.y, (double)s*v.z);
}

inline Vector3c operator * (const int& s, const Vector3c& v)
{
    return Vector3c((double)s*v.x, (double)s*v.y, (double)s*v.z);
}

inline Vector3c operator / (const Vector3c& v, const int& s)
{
    return Vector3c(v.x / (double)s, v.y / (double)s, v.z / (double)s);
}

inline Vector3c& Vector3c::operator *= (const double& s)
{
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

inline Vector3c& Vector3c::operator /= (const double& s)
{
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

inline Vector3c operator *(const Stress & s, const Vector3c & v)
{
    return Vector3c(
        s.xx*v.x + s.xy*v.y + s.xz*v.z,
        s.yx*v.x + s.yy*v.y + s.yz*v.z,
        s.zx*v.x + s.zy*v.y + s.zz*v.z);
}

inline Vector3c operator *(const Vector3c & v, const Stress & s)
{
    return Vector3c(
        v.x * ( s.xx + s.xy + s.xz),
        v.y * ( s.yx + s.yy + s.yz),
        v.z * ( s.zx + s.zy + s.zz) );
}



//======  inner product  =====

inline dcomplex operator * (const Vector3c& v1, const Vector3c& v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

inline dcomplex operator * (const Vector3c& v1, const Vector3d& v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}
inline dcomplex operator * (const Vector3d& v1, const Vector3c& v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

//======  outer product  =====

inline Vector3c operator ^ (const Vector3c& v1, const Vector3c& v2)
{
    return Vector3c(v1.y*v2.z - v1.z*v2.y,
        v1.z*v2.x - v1.x*v2.z,
        v1.x*v2.y - v1.y*v2.x);
}

//=====  Unary minus  =====

inline Vector3c Vector3c::operator - () const
{
    return Vector3c(-x, -y, -z);
}

//===== Add for Platform(STL) =====
//inline bool operator==(const Vector3c& v1, const Vector3c& v2)
//{
//    return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
//}
//inline bool operator<(const Vector3c& v1, const Vector3c& v2)
//{
//    // This definition is not universal.
//    return v1.length2() < v2.length2();
//}
//


#endif    // _Vector3c_H_
