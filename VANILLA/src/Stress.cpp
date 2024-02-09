#if defined(_CONSOLE)
#pragma warning(disable:4786)
#endif
#include <iomanip>
#include "Stress.h"
//=====    static data member of class Vector3d ================================
//
//    select output format:
//        outputParen ==    true:    ( x, y, z )
//                false:    x y z

bool Stress::outputParen = true;
bool Stress::outputXML = false;

// output precision

std::streamsize Stress::precision = UDF_DOUBLE_PRECISION;
std::streamsize Stress::reducedPrecision = 6;

ostream& operator << (ostream & os, const Stress & t) {
    int prec = (int)os.precision(Stress::getPrecision());
    if (Stress::xmlIsOn()) {
        os << "<xx>" << t.xx << "</xx><xy>" << t.xy << "</xy><xz>" << t.xz << "</xz>";
        os << "<yx>" << t.yx << "</yx><yy>" << t.yy << "</yy><yz>" << t.yz << "</yz>";
        os << "<zx>" << t.zx << "</zx><zy>" << t.zy << "</zy><zz>" << t.zz << "</zz>";
    } else if (Stress::parenIsOn()) {
        os << "{ "; 
        os << t.xx << ", " << t.xy << ", " << t.xz << ", ";
        os << t.yx << ", " << t.yy << ", " << t.yz << ", ";
        os << t.zx << ", " << t.zy << ", " << t.zz << ", ";
        os << "}" << endl;
    } else {
        os.precision(6);
        os << setw(8) << t.xx << " " << setw(8) << t.xy << " " << setw(8) << t.xz << " ";
        os << setw(8) << t.yx << " " << setw(8) << t.yy << " " << setw(8) << t.yz << " ";
        os << setw(8) << t.zx << " " << setw(8) << t.zy << " " << setw(8) << t.zz << " ";
    }
    os.precision(prec);
    return os;
}

istream& operator >> (istream& is, Stress& t)
{
    float    xx, xy, xz, yx, yy, yz, zx, zy, zz;
    char    c = 0;

    is >> c;
    if (c == '{') { is >> xx >> c; if (c == ',') { is >> xy >> c; if (c == ',') { is >> xz >> c;
    if (c == ',') { is >> yx >> c; if (c == ',') { is >> yy >> c; if (c == ',') { is >> yz >> c;
    if (c == ',') { is >> zx >> c; if (c == ',') { is >> zy >> c; if (c == ',') { is >> zz >> c;
        if (is && c == '}') {
            t.set(xx, xy, xz, yx, yy, yz, zx, zy, zz);
            return is;
        }
    }}}
    }}}
    }}} else {
        is.putback(c);
        is >> xx >> xy >> xz >> yx >> yy >> yz >> zx >> zy >> zz;
        if (is) {
            t.set(xx, xy, xz, yx, yy, yz, zx, zy, zz);
            return is;
        }
    }

    is.clear(ios::badbit);        // set badbit ON

    return is;
}

UDFistream& operator >> (UDFistream& is, Stress& t)
{
    //    float    x, y, z;
    double    xx, xy, xz, yx, yy, yz, zx, zy, zz;

    if (is.mode() == UDFistream::BINARY) {
        is.read((char*)&xx, sizeof(xx));
        is.read((char*)&xy, sizeof(xy));
        is.read((char*)&xz, sizeof(xz));
        is.read((char*)&yx, sizeof(yx));
        is.read((char*)&yy, sizeof(yy));
        is.read((char*)&yz, sizeof(yz));
        is.read((char*)&zx, sizeof(zx));
        is.read((char*)&zy, sizeof(zy));
        is.read((char*)&zz, sizeof(zz));
    }
    else {
        is >> xx >> xy >> xz >> yx >> yy >> yz >> zx >> zy >> zz;
    }
    t.set((double)xx, (double)xy, (double)xz, 
          (double)yx, (double)yy, (double)yz, 
          (double)zx, (double)zy, (double)zz);

    return is;
}

