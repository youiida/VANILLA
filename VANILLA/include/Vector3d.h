//===== Vector3d.h ===========================================================
//
//        3D vector class
//
//                ver.1.0        Aug.10, 1998    J. Takimoto
//
//===========================================================================

#ifndef _INC_VECTOR3D_H_
#define _INC_VECTOR3D_H_

#include <math.h>    // for sqrt()
#include <iostream>    // for << and >>
//ignore the setting in configenv.h and use Vector3d::precision instead
//#include "configenv.h"  // for output precision of double value
#pragma warning (disable : 4290)//nothrow
#include "udfstream.h"

#if !defined(_WITHOUT_USING_NAMESPACE_STD)
using namespace std;
#endif

//============================================================================
//
//  class Vector3d
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

class Vector3d {

public:
    double x, y, z;            // data members are public

    //----- Constructors -----

    Vector3d(double vx = 0.0, double vy = 0.0, double vz = 0.0)
                    : x(vx), y(vy), z(vz)    {}

    Vector3d(const Vector3d& v) : x(v.x), y(v.y), z(v.z) {}

    //----- get/set values -----

    void set(double vx, double vy, double vz) { x = vx; y = vy; z = vz; }

    //----- operators -----

    Vector3d& operator = (const Vector3d& v);
//    Vector3d& operator = (const double& v);

    friend Vector3d operator + (const Vector3d& v1, const Vector3d& v2);
    friend Vector3d operator - (const Vector3d& v1, const Vector3d& v2);
    Vector3d& operator += (const Vector3d& v);
    Vector3d& operator -= (const Vector3d& v);

    friend Vector3d operator * (const Vector3d& v, const double& s);
    friend Vector3d operator * (const double& s, const Vector3d &v);
    friend Vector3d operator / (const Vector3d& v, const double& s);
    friend Vector3d operator * (const Vector3d& v, const int& s);
    friend Vector3d operator * (const int& s, const Vector3d &v);
    friend Vector3d operator / (const Vector3d& v, const int& s);
    Vector3d& operator *= (const double& s);
    Vector3d& operator /= (const double& s);

    friend double operator * (const Vector3d& v1, const Vector3d& v2);
    friend Vector3d operator ^ (const Vector3d& v1, const Vector3d& v2);

    Vector3d operator - () const;

    // add for platform 2001/1/26 T.Aoyagi
    friend bool operator == (const Vector3d& v1, const Vector3d& v2);
    friend bool operator < (const Vector3d& v1, const Vector3d& v2);

    //-----  length of vector  -----

    double length() const { return sqrt(x*x + y*y + z*z); }
    double length2() const { return x*x + y*y + z*z; }

    //-----  set/reset/get outputParen  -----

    static void parenOn()    { outputParen = true; }
    static void parenOff()    { outputParen = false; }
    static bool parenIsOn()    { return outputParen; }
    static void xmlOn()    { outputXML = true; }
    static void xmlOff()    { outputXML = false; }
    static bool xmlIsOn()    { return outputXML; }

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
    Vector3d square(void) const {return Vector3d(x*x,y*y,z*z);}
    Vector3d root(void) const {return Vector3d(sqrt(x),sqrt(y),sqrt(z));}
    void clear(void){x=y=z=0.0;}
    // add 2003/10/28 T.Aoyagi
    Vector3d cube(void) const {return Vector3d(x*x*x,y*y*y,z*z*z);}
    Vector3d mp(const Vector3d& v) {return Vector3d(x*v.x, y*v.y, z*v.z);}

    //youiida 20151211
    Vector3d inv() { return Vector3d(1.0 / x, 1.0 / y, 1.0 / z); }
    Vector3d mp(const Vector3d& v1, const Vector3d & v2) {
        return Vector3d(x * v1.x * v2.x, 
                        y * v1.y * v2.y, 
                        z * v1.z * v2.z); 
    }
    //Vector3d mp(const Vector3d& v1, const Vector3d & v2, const Vector3d & v3) {
    //    return Vector3d(x * v1.x * v2.x * v3.x, 
    //                    y * v1.y * v2.y * v3.y, 
    //                    z * v1.z * v2.z * v3.z); 
    //}
    //Vector3d mp(const Vector3d& v1, const Vector3d & v2, const Vector3d & v3, const Vector3d & v4) {
    //    return Vector3d(x * v1.x * v2.x * v3.x * v4.x,
    //                    y * v1.y * v2.y * v3.y * v4.y,
    //                    z * v1.z * v2.z * v3.z * v4.z);
    //}
    Vector3d dv(const Vector3d & v) { return Vector3d(x/v.x, y/v.y, z/v.z); };

    //friend Vector3d operator * (const Stress& s, const Vector3d& v);// vector_i = stress_ij * vector_j
    //friend Vector3d operator * (const Vector3d& v, const Stress& s);// vector_j = vector_i * stress_ij

    operator double() const {return length();}
private:

    static bool outputParen;    // default value is true (Vector3d.cpp)
    static bool outputXML;    // default value is false (Vector3d.cpp)
    static std::streamsize precision;    // default UDF_DOUBLE_PRECISION (14)
    static std::streamsize reducedPrecision;    // reduced precision used for Position, Velocity and Force in Struct
};


//============================================================================
//
//    Input / Output Vector3d objects
//
//        these two functions are defined in Vector3d.cpp
//
//============================================================================


ostream& operator << (ostream& os, const Vector3d& v);
istream& operator >> (istream& is, Vector3d& v);
UDFistream& operator >> (UDFistream& is, Vector3d& v);

//============================================================================
//
//    Definition of inline operators
//
//============================================================================


//=====  Assignment  =====

inline Vector3d& Vector3d::operator = (const Vector3d& v)
{
    x = v.x; y = v.y; z = v.z;
    return *this;
}

//inline Vector3d& Vector3d::operator = (const double& d)
//{
//    x = y = z = d;
//    return *this;
//}

//=====  Vec + Vec , Vec - Vec  =====

inline Vector3d operator + (const Vector3d& v1, const Vector3d& v2)
{
    return Vector3d(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

inline Vector3d operator - (const Vector3d& v1, const Vector3d& v2)
{
    return Vector3d(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

inline Vector3d& Vector3d::operator += (const Vector3d& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

inline Vector3d& Vector3d::operator -= (const Vector3d& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

//=====  scalar*Vec, Vec*scalar, Vec/scalar  =====

inline Vector3d operator * (const Vector3d& v, const double& s)
{
    return Vector3d(s*v.x, s*v.y, s*v.z);
}

inline Vector3d operator * (const double& s, const Vector3d& v)
{
    return Vector3d(s*v.x, s*v.y, s*v.z);
}

inline Vector3d operator / (const Vector3d& v, const double& s)
{
    return Vector3d(v.x/s, v.y/s, v.z/s);
}

inline Vector3d operator * (const Vector3d& v, const int& s)
{
    return Vector3d((double)s*v.x, (double)s*v.y, (double)s*v.z);
}

inline Vector3d operator * (const int& s, const Vector3d& v)
{
    return Vector3d((double)s*v.x, (double)s*v.y, (double)s*v.z);
}

inline Vector3d operator / (const Vector3d& v, const int& s)
{
    return Vector3d(v.x/(double)s, v.y/(double)s, v.z/(double)s);
}

inline Vector3d& Vector3d::operator *= (const double& s)
{
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

inline Vector3d& Vector3d::operator /= (const double& s)
{
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

//======  inner product  =====

inline double operator * (const Vector3d& v1, const Vector3d& v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

//======  outer product  =====

inline Vector3d operator ^ (const Vector3d& v1, const Vector3d& v2)
{
    return Vector3d(v1.y*v2.z - v1.z*v2.y ,
            v1.z*v2.x - v1.x*v2.z ,
            v1.x*v2.y - v1.y*v2.x );
}

//=====  Unary minus  =====

inline Vector3d Vector3d::operator - () const
{
    return Vector3d(-x, -y, -z);
}

//===== Add for Platform(STL) =====
inline bool operator==(const Vector3d& v1, const Vector3d& v2)
{
    return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
}
inline bool operator<(const Vector3d& v1, const Vector3d& v2) 
{
    // This definition is not universal.
    return v1.length2() < v2.length2();
}



#endif    // _Vector3d_H_
