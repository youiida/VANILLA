#ifndef _INC_DATAIO_H_
#define _INC_DATAIO_H_

#pragma warning (disable : 4290)//nothrow
#include "udfmanager.h"
#include "MyTools.h"
#include "Field.h"
#include "Dynamics.h"
#include "FreeEnergy.h"
#include "Reaction.h"

class DataIO {
private:

    //bool    readScalarField(UDFManager & uobj,
    //                        Location cursor,
    //                        ScalarField & sf);
    //void    putScalarField(UDFManager & uobj,
    //                        Location cursor,
    //                        const ScalarField & sf);

    //some bugs in Vector/StressField to make template <class T>...
    //20160121 youiida
    bool    readVectorField(UDFManager & uobj,
                            Location cursor,
                            VectorField & vf);
    void    putVectorField(UDFManager & uobj,
                            Location cursor,
                            const VectorField & vf);
    bool    readStressField(UDFManager & uobj,
                            Location cursor,
                            StressField & sf);
    void    putStressField(UDFManager & uobj,
                            Location cursor,
                            const StressField & sf);
    template <class T>
    bool    readField(UDFManager & uobj,
        Location cursor,
        T & sf);
    template <class T>
    void    putField(UDFManager & uobj,
        Location cursor,
        const T & sf);

    void    readSimulationConditions(UDFManager & uobj);
    void    readDynamics(UDFManager & uobj,
                         Dynamics* & dynamics,
                         vector<Component> & components);
    void    readFreeEnergyType(UDFManager & uobj,
                               FreeEnergy* & energy,
                               vector<Component> & components);
    void    readComponentProperties(UDFManager & uobj,
                                    vector<Component> & components);
    void    readInitialConditions(UDFManager & uobj,
                                  vector<ScalarField> & phi,
                                  vector<Component> & components);
    void    readReactions(UDFManager & uobj,
                          Reaction & reaction);

    Segment    find(string name, vector<Segment> & seg) {
        for (int i = 0; i < (int)seg.size(); i++) {
            if (seg[i].name == name)    return seg[i];
        }
        return Segment("not found");
    };
    vector<int>        find(string _name, vector<Component> & comps) {
        vector<int>    index;
        for (int i = 0; i < (int)comps.size(); i++) {
            if (comps[i].get_segname() == _name)    index.push_back(i);
        }
        return index;
    }

    void    randomInit(vector<Component> & _comps,
                       vector<ScalarField> & _phi,
                       string seed,
                       double sigma);

    Mesh            mesh;
    double          dt;
    int             totalRecords, interval, M;

    int             init_step;
    double          init_time;

    vector<double>  chi;

    fstream         lobj;
    void            open(const string & log_file) { lobj.open(log_file.c_str(), ios::out); };
    string          start_condition;
    bool            vector_flag, stress_flag;
public:
    DataIO(const string & in_udf,
           const string & out_udf,
           const string & log_file) {
        open(log_file);
        init_step   = 0;
        init_time   = 0.0;
        vector_flag = false;
        stress_flag = false;
        string head =  "========== File Names ==========\n";
        string  msg =  "Input UDF .... " + in_udf + "\n";
                msg += "Output UDF ... " + out_udf + "\n";
                msg += "Log File ..... " + log_file + "\n";
        logger(head + msg);
    };
    ~DataIO() {};
    void    setup(UDFManager & uobj,
                  vector<Component> & comps,
                  vector<ScalarField> & phi,
                  Dynamics* & dynamics,
                  FreeEnergy* & energy,
                  Reaction & reaction,
                  vector<ScalarField> & mu,
                  int omp_threads);
    void    write(UDFManager & uobj,
                  const int record,
                  const double free_energy,
                  vector<Component> & comps,
                  vector<ScalarField> & phi,
                  Dynamics* & dynamics,
                  const string & lap_time,
                  const string & elapse_time,
                  const string & remaining_time,
                  const string & comment = "record");
    int        getTotalRecords(void) { return totalRecords; };
    int        getInterval(void) { return interval; };
    double     getIntervalTime(void){return interval*dt;};
    string     getStartCondition(void) { return start_condition; };
    void       logger(const string & str) { cout << str << endl; lobj << str << endl; };
    bool       variable_dt_flag;
};
#endif