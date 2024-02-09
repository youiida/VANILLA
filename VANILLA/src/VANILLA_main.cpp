#include "MyTools.h"
#include "Field.h"
#include "Component.h"
#include "FreeEnergy.h"
#include "Dynamics.h"
#include "DataIO.h"
#include "MyUDFIO.h"
#include "Reaction.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void small_dt_warning(UDFManager & uobj, DataIO & io, FreeEnergy* & energy,
                      Dynamics* & dynamics, vector<Component> & component,
                      vector<ScalarField> & phi, int record, double dt);

int main(int argc, char* argv[]) {
    Timer    timer;
    clock_t  present;
    string   in_udf, out_udf, log_file;
    int      omp_threads = 1;
    double   dt, interval_time, time;
    int      counter;
    string   msg;

    vector<Component>      component;
    vector<ScalarField>    phi;
    vector<ScalarField>    mu;

    Dynamics*      dynamics;
    FreeEnergy*    energy;

    Reaction       reaction;

    parseArg(argc, argv, in_udf, out_udf, log_file, omp_threads);

    UDFManager    in_uobj(in_udf);
    UDFManager    out_uobj(out_udf, in_udf);
    #ifdef _OPENMP
        omp_set_num_threads(omp_threads);
    #endif


    try {
        // read udf and setup
        DataIO    io(in_udf, out_udf, log_file);
        io.setup(in_uobj,
                 component,
                 phi,
                 dynamics,
                 energy,
                 reaction,
                 mu,
                 omp_threads);

        int total_records, interval;
        total_records    = io.getTotalRecords();
        interval         = io.getInterval();
        interval_time    = io.getIntervalTime();
        present          = clock();

        dynamics->logSetting(in_udf);
        energy->update(component, phi, mu);

        io.write(out_uobj,
                 0,
                 energy->get(component, phi),
                 component,
                 phi,
                 dynamics,
                 timer.lapTime(present),
                 timer.elapseTime(present),
                 timer.estTime(present, 0, total_records),
                 io.getStartCondition() );

        // start time evolution
        for (int record = 1; record <= total_records; record++) {
            time    = 0.0;
            counter = 0;
            if (io.variable_dt_flag) {
                while (time <= interval_time) {
                    // dynamics main
                    energy->update(component, phi, mu);
                    dt = dynamics->evolveWithChangeDt(phi, mu, component);

                    // reaction-depend block
                    // 20160715 youiida
                    // adapt react, updateVertex, updateIsPolymer
                    // to updateDensityWithChangeDt()
                    // -> DONE 0719
                    reaction.react(phi, component, dt);
                    energy->updateVertex(component);
                    dynamics->updateIsPolymer(component);

                    time += dt;
                    counter++;
                    if (counter > interval * 10) {
                        small_dt_warning(out_uobj, io, energy,
                                         dynamics, component, phi, record, dt);
                        counter = 0;
                    }
                }
            } else {
                // change from <for-loop> to <while> to introduce variable-dt
                // 20160623 youiida
                // select by flag
                // 20160810 youiida
                dt  = dynamics->get_dt();
                for (int inter = 0; inter < interval; inter++) {
                    // if (reaction.flag(inter)) {
                        reaction.react(phi, component, dt);
                        energy->updateVertex(component);
                        dynamics->updateIsPolymer(component);
                    // }
                    energy->update(component, phi, mu);
                    dynamics->evolve(phi, mu, component);
                }
            }

            // write output udf/bdf
            present = clock();
            io.write(out_uobj,
                     record,
                     energy->get(component, phi),
                     component,
                     phi,
                     dynamics,
                     timer.lapTime(present),
                     timer.elapseTime(present),
                     timer.estTime(present, record, total_records));

            // get information from dynamics 20160623 youiida
            // Mainly, to detect numerical divergence of the density field
            if (dynamics->logger(record, phi)) {
                io.logger(dynamics->passError());
                exit(1);
            }
        }

        // normal end
        delete dynamics;
        delete energy;
        io.logger(timer.finish());
    } catch(...) {
        cerr << "Undefined error";
        return 2;
    }

    return 0;
}

void small_dt_warning(UDFManager & uobj,
                      DataIO & io,
                      FreeEnergy* & energy,
                      Dynamics* & dynamics,
                      vector<Component> & component,
                      vector<ScalarField>    & phi,
                      int record, double dt) {
    string msg = "small dt:" + dtos(dt);
    io.logger(msg);
    clock_t present = clock();
    io.write(uobj,
             record,
             energy->get(component, phi),
             component,
             phi,
             dynamics,
             "temporary",
             "temporary",
             "temporary",
             "temporary");
    return;
}

