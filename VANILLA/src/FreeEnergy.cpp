//#include "FreeEnergy.h"

//double FHdG::get(vector<Component> & comps,
//                 vector<ScalarField> & phi) {
//    double        fe            = 0.0;
//    double        uni_entropy    = 0.0;
//    double        uni_enthalpy= 0.0;
//    double        ninv;
//
//    for (int p = 0; p < M; p++) {
//        ninv            =  1.0 / comps[p].N;
//        uni_entropy        += ninv * comps[p].vf * log(comps[p].vf);
//        for (int p2 = 0; p2 < M; p2++) {
//            uni_enthalpy += chi[M*p + p2] * comps[p].vf * comps[p2].vf;
//        }
//
//        for (int ix = 0; ix < size; ix++) {
//            fe += ninv*(phi[p][ix] * log(phi[p][ix]))
//                  + comps[p].get_tension() * ( phi[p].grad(ix) ).length2();
//            for (int p2 = 0; p2 < M; p2++){
//                fe += chi[M*p + p2] * phi[p][ix] * phi[p2][ix];
//            }
//        }
//    }
//
//    return fe*vinv - (uni_entropy+uni_enthalpy);
//}

//void FHdG::update(vector<Component> & comps,
//                  vector<ScalarField> & phi,
//                  vector<ScalarField> & mu) {
//    double ninv;
//    for (int p = 0; p < M; p++) {
//        ninv = 1.0 / comps[p].N;
//        for (int ix = 0; ix < size; ix++) {
//            mu[p][ix] = ninv * (log(phi[p][ix]) + 1.0) 
//                        - 2.0*comps[p].get_tension()*phi[p].lap(ix);
//            for (int p2 = 0; p2 < M; p2++) {
//                mu[p][ix] += chi[M*p + p2] * phi[p2][ix];
//            }
//        }
//    }
//    return;
//}