#ifndef _INC_REACTION_H_
#define _INC_REACTION_H_

#include "Field.h"
#include "Component.h"

class BaseReaction {
protected:
    int        interval;
    int        size;//nx*ny*nz
    string     type;
    void       baseSet(int _size,
                       int _interval) {
        size     = _size;
        interval = _interval;
        return;
    }
public:
    BaseReaction() { type = "BaseReaction(NULL)"; };
    ~BaseReaction() {};
    virtual void    react(vector<ScalarField> & phi,
                          vector<Component> & comps) {
        return;
    };
    virtual void    react(vector<ScalarField> & phi,
                          vector<Component> & comps,
                          double dt) {
        return;
    };
    virtual string    output(int id) {
        string head = "========== Reaction Condition (ID:" + itos(id) + ")==========\n";
        string mid  = "Type ................ " + type + "\n";
              // mid  +="Reaction Interval ... " + itos(interval) +"\n";
        return head + mid;
    };
};

class FirstOrder : public BaseReaction {
protected:
    double    k1; // k1
    int       react_id;
    int       product_id;
    // A(react) -> B(product)
    // d[B]/dt =  k1[A]    phi_B += k1*phi_A * dt
    // d[A]/dt = -k1[A]   phi_A -= k1*phi_A * dt

public:
    FirstOrder() {};
    FirstOrder(int _react_id, 
                int _product_id,
                double _k1,
                int _interval,
                int _size) {
        baseSet(_size, _interval);
        react_id   = _react_id;
        product_id = _product_id;
        k1         = _k1;
        type       = "First Order Reaction";
    };
    ~FirstOrder() {};
    //void    react(vector<ScalarField> & phi,
    //                vector<Component> & comps) {
    //    double sum = 0.0;
    //    double val;
    //    for (int ix = 0; ix < size; ix++) {
    //        val = phi[react_id][ix] * k1*dt * (double)interval;
    //        phi[react_id][ix]    -= val;
    //        phi[product_id][ix] += val;
    //        sum += val;
    //    }
    //    sum /= volume;
    //    comps[react_id].vf        -= sum;
    //    comps[product_id].vf    += sum;
    //    return;
    //};
    void    react(vector<ScalarField> & phi,
                  vector<Component> & comps,
                  double dt) {
        double sum = 0.0;
        double val;
        for (int ix = 0; ix < size; ix++) {
            val = phi[react_id][ix] * k1 * dt;
            phi[react_id][ix]   -= val;
            phi[product_id][ix] += val;
            sum                 += val;
        }
        sum /= double(size);
        comps[react_id].vf   -= sum;
        comps[product_id].vf += sum;
        return;
    };
    string output(int id) {
        string    head   = BaseReaction::output(id);
        string    mid    =  "Reactant ID ......... " + itos(react_id) + "\n";
                  mid    += "Product ID .......... " + itos(product_id) + "\n";
                  mid    += "k1 .................. " + dtos(k1) + "\n";
        return head + mid;
    };
};

class SecondOrder : public BaseReaction {
protected:
    int        react_id1, react_id2, product_id;
    double    k2;
    // A1(react1) + A1(react2) -> B(product)
    // d[B]/dt =   k2[A1][A2]    phi_B  += k2*phi_A1*phi_A2 * dt
    // d[A1]/dt = -k2[A1][A2]    phi_A1 -= k2*phi_A1*phi_A2*V_A1/(V_A1+V_A2) * dt
    // d[A2]/dt = -k2[A1][A2]    phi_A2 -= k2*phi_A1*phi_A2*V_A2/(V_A1+V_A2) * dt
public:
    SecondOrder() {};
    SecondOrder(int _react_id1,
                int _react_id2,
                int _product_id,
                double _k2,
                int _interval,
                int _size) {
        baseSet(_size, _interval);
        react_id1  = _react_id1;
        react_id2  = _react_id2;
        product_id = _product_id;
        k2         = _k2;
        type       = "Second Order Reaction";
    };
    ~SecondOrder() {};

    //void react(vector<ScalarField> & phi,
    //            vector<Component> & comps){
    //    double sum = 0.0;
    //    double val;
    //    double rate = comps[react_id1].get_volume() / (comps[react_id1].get_volume() + comps[react_id2].get_volume());
    //    for (int ix = 0; ix < size; ix++) {
    //        val = phi[react_id1][ix] * phi[react_id2][ix] * k2dt * (double)interval;
    //        phi[react_id1][ix] -= rate*val;
    //        phi[react_id2][ix] -= (1.0-rate)*val;
    //        phi[product_id][ix] -= val;
    //        sum += val;
    //    }
    //    sum /= volume;
    //    comps[react_id1].vf -= sum * rate;
    //    comps[react_id2].vf -= sum * (1.0 - rate);
    //    comps[product_id].vf += sum;
    //    return;
    //};
    void react(vector<ScalarField> & phi,
               vector<Component> & comps,
               double dt) {
        double sum = 0.0;
        double val;
        double rate = comps[react_id1].get_volume() / (comps[react_id1].get_volume() + comps[react_id2].get_volume());
        for (int ix = 0; ix < size; ix++) {
            val = phi[react_id1][ix] * phi[react_id2][ix] * k2 * dt;
            phi[react_id1][ix]  -= rate*val;
            phi[react_id2][ix]  -= (1.0 - rate)*val;
            phi[product_id][ix] -= val;
            sum += val;
        }
        sum /= double(size);
        comps[react_id1].vf  -= sum * rate;
        comps[react_id2].vf  -= sum * (1.0 - rate);
        comps[product_id].vf += sum;
        return;
    };
    string output(int id) {
        string    head   =  BaseReaction::output(id);
        string    mid    =  "Reactant1 ID ........ " + itos(react_id1) + "\n";
                  mid    += "Reactant2 ID ........ " + itos(react_id2) + "\n";
                  mid    += "Product ID .......... " + itos(product_id) + "\n";
                  mid    += "k2 .................. " + dtos(k2) + "\n";
        return head + mid;
    };
};

class ChainGrowth : public BaseReaction {
protected:
    int       react_id, product_id;
    double    k;
    // A(react) -> B(product)
    // d[B]/dt =  k1[A]   phi_B += k1*phi_A * dt
    // d[A]/dt = -k1[A]   phi_A -= k1*phi_A * dt
    // V_i : segment volume
    // N_i : degree of polymerization
    // M_i : number of polymers in unit volume
    // 
    // N_A*M_A*V_A = 1, N_A*M_V_A = phi
    // dN = N_A*M = N_A*M_A*phi = phi/V_A


public:
    ChainGrowth() {};
    ChainGrowth(int _react_id,
                int _product_id,
                double _k,
                int _interval,
                int _size) {
        baseSet(_size, _interval);
        react_id   = _react_id;
        product_id = _product_id;
        k          = _k;
        type       = "Chain Growth Reaction";
    };
    ~ChainGrowth() {};
    string output(int id) {
        string    head = BaseReaction::output(id);
        string    mid =  "Reactant ID ......... " + itos(react_id) + "\n";
                  mid += "Product ID .......... " + itos(product_id) + "\n";
                  mid += "kp .................. " + dtos(k) + "\n";
        return head + mid;
    };

    void react(vector<ScalarField> & phi,
                vector<Component> & comps,
                double dt) {
        double sum = 0.0;
        double val;
        for (int ix = 0; ix < size; ix++) {
            val = phi[react_id][ix] * k*dt * (double)interval;
            phi[react_id][ix]   -= val;
            phi[product_id][ix] += val;
            sum += val;
        }
        sum /= double(size);

        comps[product_id].N  += sum;
        comps[react_id].vf   -= sum;
        comps[product_id].vf += sum;
        return;
    };
    //void react(vector<ScalarField> & phi,
    //    vector<Component> & comps) {
    //    double sum = 0.0;
    //    double val;
    //    for (int ix = 0; ix < size; ix++) {
    //        val = phi[react_id][ix] * k*dt * (double)interval;
    //        phi[react_id][ix] -= val;
    //        phi[product_id][ix] += val;
    //        sum += val;
    //    }
    //    sum /= volume;
    //    comps[product_id].N += sum;
    //    comps[react_id].vf -= sum;
    //    comps[product_id].vf += sum;
    //    return;
    //};
};

class Reaction {
private:
    int    num_react, interval;
public:
    Reaction() { num_react = 0; };
    ~Reaction() {
        for (int i = 0; i < num_react; i++) {
            delete reactions[i];
        }
    };
    void    set(int _num_react,
                int _interval) {
        num_react = _num_react;
        reactions.resize( num_react );
        //flags.resize( num_react );
        interval = _interval;
    }
    vector<BaseReaction*>    reactions;
    //vector<bool>            flags;

    void    react(vector<ScalarField> & phi,
                  vector<Component> & comps,
                  double dt) {
        for (int i = 0; i < num_react; i++) {
                reactions[i]->react( phi, comps, dt );
        }
        return;
    };
    //void    react(vector<ScalarField> & phi,
    //    vector<Component> & comps) {
    //    for (int i = 0; i < num_react; i++) {
    //        reactions[i]->react(phi, comps);
    //    }
    //    return;
    //};
    //bool    flag(int dyn_counter) {
    //    if (num_react == 0)
    //        return false;
    //    if (dyn_counter%interval == 0)
    //        return true;
    //    else return false;
    //};

    //bool    flag(int dyn_conuter) {
    //    bool react_flag = false;
    //    for (int i = 0; i < num_react; i++) {
    //        if (reactions[i]->flag(dyn_conuter)) {
    //            flags[i] = true;
    //            react_flag = true;
    //        } else {
    //            flags[i] = false;
    //        }
    //    }
    //    return react_flag;
    //};

    string output(void) {
        string out;
        if (num_react == 0) {
            out =  "========== Reaction Condition ==========\n";
            out += "No Reactions\n";
        } else {
            for (int i = 0; i < num_react; i++) {
                out += reactions[i]->output(i);
            }
        }
        return out;
    };
};

#endif