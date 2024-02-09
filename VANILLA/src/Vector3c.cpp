#if defined(_CONSOLE)
#pragma warning(disable:4786)
#endif
#include <iomanip>
//===== Vector3c.cpp ========================================================
//
//        Input / Output Vector3c object
//
//
//                ver.1.0        Aug.10, 1998    J. Takimoto
//
//===========================================================================

#include "Vector3c.h"

//=====    static data member of class Vector3c ================================
//
//    select output format:
//        outputParen ==    true:    ( x, y, z )
//                false:    x y z

bool Vector3c::outputParen = true;
bool Vector3c::outputXML = false;

// output precision

std::streamsize Vector3c::precision = UDF_DOUBLE_PRECISION;
std::streamsize Vector3c::reducedPrecision = 6;

//=====  output  ============================================================

ostream& operator << (ostream& os, const Vector3c& v)
{
    int prec = (int)os.precision(Vector3c::getPrecision());
    if (Vector3c::xmlIsOn()) {
        os << "<x>" << v.x << "</x><y>" << v.y << "</y><z>" << v.z << "</z>";
    }
    else if (Vector3c::parenIsOn()) {
        os << "{ " << v.x << ", " << v.y << ", " << v.z << " }" << endl;
    }
    else {
        os.precision(6);
        os << setw(8) << v.x << " " << setw(8) << v.y << " " << setw(8) << v.z;
    }
    os.precision(prec);
    return os;
}

//=====  input ==============================================================
//
//    both    ( x, y, z )    and    x y z    are acceptable
//
//
//        see "The C++ Programing Language (3rd ed.)" by B.Stroustrup,
//            Section 21.3.5
//
//istream& operator >> (istream& is, Vector3c& v)
//{
//    float    x, y, z;
//    char    c = 0;
//
//    is >> c;
//    if (c == '{') {
//        is >> x >> c;
//        if (c == ',') {
//            is >> y >> c;
//            if (c == ',') {
//                is >> z >> c;
//                if (is && c == '}') {
//                    v.set(x, y, z);
//                    return is;
//                }
//            }
//        }
//    }
//    else {
//        is.putback(c);
//        is >> x >> y >> z;
//        if (is) {
//            v.set(x, y, z);
//            return is;
//        }
//    }
//
//    is.clear(ios::badbit);        // set badbit ON
//
//    return is;
//}
//
//UDFistream& operator >> (UDFistream& is, Vector3c& v)
//{
//    //    float    x, y, z;
//    double    x, y, z;
//
//    if (is.mode() == UDFistream::BINARY) {
//        is.read((char*)&x, sizeof(x));
//        is.read((char*)&y, sizeof(y));
//        is.read((char*)&z, sizeof(z));
//    }
//    else {
//        is >> x >> y >> z;
//    }
//    v.set((double)x, (double)y, (double)z);
//
//    return is;
//}
//
