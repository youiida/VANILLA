#ifndef _INC_COMPONENT_H_
#define _INC_COMPONENT_H_

#include "MyTools.h"
using namespace std;

//youiida 20160301 add gb,gs,tb,ts
class Segment {
private:
protected:

public:
    Segment() {};
    Segment(const string _name) { name = _name; };
    //youiida 20160301
    //Segment(const string _name,
    //        const double _size,
    //        const double _zeta,
    //        const double _tension)
    //        : name(_name), size(_size), zeta(_zeta), tension(_tension) {};
    void    set(const string _name,
                const double _size,
                const double _gb,
                const double _gs,
                const double _tb,
                const double _ts,
                const double _zeta,
                //const double _eta,
                const double _tension) {
        name    = _name;
        size    = _size;
        gb      = _gb;
        gs      = _gs;
        tb      = _tb;
        ts      = _ts;
        zeta    = _zeta;
        tension = _tension;
        //eta        = gs/ts;
        return;
    };
    string    output(const int id) {
        //youiida 20160301
        string head = "========== Segment Condition (ID:" + itos(id) + ")==========\n";
        string  mid = "Name .................... " + name + "\n";
                mid +="Size .................... " + dtos(size) + "\n";
                mid +="Bulk Modulus ............ " + dtos(gb) + "\n";
                mid +="Shear Modulus ........... " + dtos(gs) + "\n";
                mid +="Bulk Relaxation Time .... " + dtos(tb) + "\n";
                mid +="Shear Relaxation Time ... " + dtos(ts) + "\n";
                mid +="Friction ................ " + dtos(zeta) + "\n";
                mid +="Interfacial Tension ..... " + dtos(tension) + "\n";
        return head + mid;
    };
    string    name;
    double    size;//segment size, b
    double    zeta;//friction
    //double    eta;
    double    tension;//interfacial tension (for FHdG energy)

    //youiida 20160301
    //double    eta;//viscosity calc from modulus and tau
    double    gb, gs, tb, ts;//bulk, shear modulus and tau
};

class Component {
private:
protected:
    string    name;
    Segment   seg;
    double    mzeta;//polymer -> 0, monomer -> zeta0
    double    pzeta;//polymer -> N*zeta0/Ne, monomer -> 0
    double    volume;//unit volume for reaction
    //polymer N > Ne
    //monomer N < Ne

public:
    Component() {};
    Component(const double _N,
              const double _vf,
              const string _name,
              const Segment _seg)
              : N(_N), vf(_vf), name(_name), seg(_seg) {};
    ~Component() {};

    double    N;// degree of polymerization
    double    vf;//volume fraction

    double    mz() { return mzeta; };
    double    pz() { return pzeta; };
    double    get_tension() { return seg.tension; };
    double    get_b() { return seg.size; };
    double    get_gb() { return seg.gb; };
    double    get_gs() { return seg.gs; };
    double    get_tb() { return seg.tb; };
    double    get_ts() { return seg.ts; };
//    double    get_eta() { return seg.eta; };
    double    get_volume() { return volume; };
    string    get_segname() { return seg.name; };
    string    get_name() { return name; };
    //for visco elastic dynamics
    void    set_mpzeta(double Ne) {
        if (N > Ne) {
            pzeta = N*seg.zeta/Ne;
            mzeta = 0.0;
        }
        else {
            pzeta = 0.0;
            mzeta = seg.zeta;
        }
    };

    void    set(const double _N,
                const double _vf,
                const string _name,
                const Segment _seg) {
        N       = _N;
        vf      = _vf;
        name    = _name;
        seg     = _seg;
        volume  = 1.0;//temporary setting
        return;
    };
    string    output(const int id) {
        string    head = "========== Component Condition (ID:" + itos(id) + ")==========\n";
        string    mid  = "Name ....................... " + name + "\n";
                  mid += "Volume Fraction ............ " + dtos(vf) + "\n";
                  mid += "Degree of Polymerization ... " + dtos(N) + "\n";
                  mid += "Segment Name ............... " + seg.name + "\n";
                  mid += "Unit Volume ................ " + dtos(volume) + "\n";
        return head + mid;
    };
};

#endif