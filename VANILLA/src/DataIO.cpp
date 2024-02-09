#include "DataIO.h"

void DataIO::setup(UDFManager & uobj,
                   vector<Component> & comps,
                   vector<ScalarField> & phi,
                   Dynamics* & dynamics,
                   FreeEnergy* & energy,
                   Reaction & reaction,
                   vector<ScalarField> & mu,
                   int omp_threads) {
    try {
        readSimulationConditions( uobj );
        
        readComponentProperties( uobj, comps );
        
        readDynamics( uobj, dynamics, comps );
        
        readFreeEnergyType( uobj, energy, comps );

        readReactions(uobj, reaction);

        readInitialConditions( uobj, phi, comps );
        

        mu.resize(M);
        for (int i = 0; i < M; i++) {
            mu[i].resize(mesh);
        }
        if (energy->get_name() == "Random_Phase_Approximation") {
            for (int p = 0; p < M; p++) {
                phi[p].fftSetup();
            }
        }

        if (start_condition == "Restart") {
            Location cursor("output.components[]");
            for (int p = 0; p < M; p++) {
                cursor.next();
                readVectorField(uobj, cursor.sub("velocity[]"), dynamics->velocity[p]);
                readStressField(uobj, cursor.sub("shear_stress[]"), dynamics->shear[p]);
                readField(uobj, cursor.sub("bulk_stress[]"), dynamics->bulk[p]);
            }
        }

    } catch (PFException & e) {
        cerr << e.GetErrorMessage() << endl;
        exit(1);
    }

#ifdef _OPENMP
    logger("========== OpenMP ==========");
    logger("Number of Threads / Max Threads ... " + itos(omp_threads) + "/" + itos(omp_get_num_procs())+"\n");
#endif
    
    logger("========================= Successful Initialization =========================\n\n");
    logger(" [Record]        [Lap]            [Elapse]            [Remain]   [Free Energy]");
    return;
}

void DataIO::readInitialConditions(UDFManager & uobj,
                                   vector<ScalarField> & phi,
                                   vector<Component> & comps) {
    Location    cursor("initial_conditions");
    string      type;
    uobj.get(cursor.sub("density_initialize"), type);

    phi.resize(M);
    for (int i = 0; i < M; i++) {
        phi[i].resize(mesh);
    }

    string  seed;
    double  sigma;

    if (type == "random") {
        Location    rcur = cursor.sub("random");

        uobj.get(rcur.sub("random_seed"), seed);
        uobj.get(rcur.sub("random_sigma"), sigma);

        randomInit(comps, phi, seed, sigma);
        start_condition = "Initial";

    } else if (type == "restart") {
        uobj.get("initial_conditions.restart.restart_record", init_step);
        uobj.jump(init_step);
        uobj.get("output.time", init_time);
        
        Location    cursor("output.components[]");
        for (int p = 0; p < M; p++) {
            cursor.next();
            readField(uobj, cursor.sub("density[]"), phi[p]);
        }

        start_condition = "Restart";
    } else {
        warning("readInitialConditions"+type);
        exit(1);
    }

    string head = "========== Density Field Initialization ==========\n";
    string mid = "Initialize by ... " + type + "\n";
    if (type == "random") {
        mid += "Random seed ..... " + seed + "\n";
    }
    logger(head + mid);

    return;
}

void DataIO::randomInit(vector<Component> & comps,
                        vector<ScalarField> & phi,
                        string seed,
                        double sigma) {
    //setup random number
    GaussianRandomNumber    gr;

    if (seed == "constant") {
        gr.setSeed( 12345 );
    } else {
        gr.setSeed( (unsigned)time(NULL) );
    }

    double    min        = 1.0;
    int        maxid    = 0;
    double    max        = -1.0;

    for (int i = 0; i < M; i++) {
        if (min > comps[i].vf) min = comps[i].vf;
        if (max < comps[i].vf) {
            maxid    = i;
            max        = comps[i].vf;
        }
    }
    gr.setSigma( float(min*sigma) );

    
    //initialize
    double    sum, mean, vf;
    int        size = mesh.nx*mesh.ny*mesh.nz;


    for (int i = 0; i < M; i++) {
        vf = comps[i].vf;
        if (i != maxid) {
            sum = 0.0;
            for (int ix = 0; ix < size; ix++) {
                //phi[i][ix] =    comps[i].vf *(1+0.1*sin(ix*PI2/ mesh.nx));
                phi[i][ix]   =  comps[i].vf + gr();
                sum            +=    phi[i][ix];
            }
            mean = sum / double(size) - vf;
            for (int ix = 0; ix < size; ix++) {
                phi[i][ix] -= mean;
                if (phi[i][ix]>1.0 || phi[i][ix] < 0.0) {
                    warning("too large sigma!?");
                    exit(1);
                }
            }
        }
    }
    for (int ix = 0; ix < size; ix++) {
        sum = 0.0;
        for (int i = 0; i < M; i++) {
            if (i != maxid) {
                sum += phi[i][ix];
            }
        }
        phi[maxid][ix] = 1.0 - sum;
        if (phi[maxid][ix]>1.0 || phi[maxid][ix] < 0.0) {
            warning("too large sigma!?");
            exit(1);
        }
    }

    return;
}


void DataIO::readFreeEnergyType(UDFManager & uobj,
                                FreeEnergy* & energy,
                                vector<Component> & comps) {
    Location    cursor("free_energy");
    string      model;
    uobj.get( cursor.sub("model"), model );

    string head = "========== Free Energy ==========\n";
    string mid  =  "Free Energy Type ... " + model + "\n";
            

    if (model == "Flory_Huggins_deGennes") {
        energy = new FHdG(mesh, M, chi);
    } else if(model == "Random_Phase_Approximation") {
        //double    kappa, alpha;
        //uobj.get(cursor.sub("Random_Phase_Approximation.compressibility"), kappa);
        //alpha = 1.0 / kappa;
        energy = new RPA(mesh, M, /*alpha,*/ comps, chi);

        //mid += "Compressibility .... " + dtos(kappa) + "\n";
    } else {
        warning("readFreeEnergyTpe:" + model);
        exit(1);
    }

    logger(head + mid);
    return;
}

void DataIO::readComponentProperties(UDFManager & uobj,
                                     vector<Component> & comps) {
    Location        cursor("component_properties");

    //read segment
    Location        scur = cursor.sub("segment[]");
    int             segs = uobj.size(scur);
    vector<Segment> segments(segs);
    string          name;
    double          size, zeta, gb, gs, tb, ts, tension;


    for (int i = 0; i < segs; i++) {
        scur.next();
        uobj.get(scur.sub("name"), name);
        uobj.get(scur.sub("size"), size);
        uobj.get(scur.sub("bulk_modulus"), gb);
        uobj.get(scur.sub("shear_modulus"), gs);
        uobj.get(scur.sub("bulk_relaxation_time"), tb);
        uobj.get(scur.sub("shear_relaxation_time"), ts);
        uobj.get(scur.sub("friction"), zeta);
        uobj.get(scur.sub("interfacial_tension"), tension);
        segments[i].set(name, size, gb, gs, tb, ts, zeta, tension);
        logger( segments[i].output(i) );
    }

    //read components
    Location    ccur = cursor.sub("components[]");
    string      segname;
    Segment     seg;
    double      N, vf;
    M = uobj.size(ccur);
    comps.resize(M);

    for (int i = 0; i < M; i++) {
        ccur.next();
        uobj.get(ccur.sub("name"), name);
        uobj.get(ccur.sub("volume_fraction"), vf);
        uobj.get(ccur.sub("segment"), segname);
        uobj.get(ccur.sub("degree_of_polymerization"), N);
        seg = find(segname, segments);
        comps[i].set(N, vf, name, seg);
        logger( comps[i].output(i) );
    }


    //read chi-parameter
    Location    chicur = cursor.sub("chi_parameter[]");
    chi.resize(M*M);
    for (int i = 0; i < M*M; i++) {
        chi[i] = 0.0;
    }
    vector<int>   iv1, iv2;
    string        s1, s2;
    double        value;

    string head = "========== Chi parameter ==========\n";

    for (int i = 0; i < (int)uobj.size(chicur); i++) {
        chicur.next();
        uobj.get(chicur.sub("segment_name1"), s1);
        uobj.get(chicur.sub("segment_name2"), s2);
        uobj.get(chicur.sub("parameter"), value);
        iv1 = find(s1, comps);
        iv2 = find(s2, comps);
        for (int i1 = 0; i1 < (int)iv1.size(); i1++) {
            for (int i2 = 0; i2 < (int)iv2.size(); i2++) {
                chi[iv1[i1]*M + iv2[i2]] = value;
                chi[iv2[i2]*M + iv1[i1]] = value;
                head += "chi(" + itos(iv1[i1]) + "," + itos(iv2[i2]) + ") ... " + dtos(value) +"\n";
                head += "chi(" + itos(iv2[i2]) + "," + itos(iv1[i1]) + ") ... " + dtos(value) + "\n";
            }
        }
    }
    logger(head);
    return;
}

void DataIO::readSimulationConditions(UDFManager & uobj) {
    
    Location    cursor("simulation_conditions");

    //mesh_conditions
    Location    mcur = cursor.sub("mesh_conditions");
    int         nx, ny, nz;
    double      dx, dy, dz;
    uobj.get(mcur.sub("number_of_mesh.Nx"), nx);
    uobj.get(mcur.sub("number_of_mesh.Ny"), ny);
    uobj.get(mcur.sub("number_of_mesh.Nz"), nz);
    uobj.get(mcur.sub("mesh_size.dx"), dx);
    uobj.get(mcur.sub("mesh_size.dy"), dy);
    uobj.get(mcur.sub("mesh_size.dz"), dz);
    mesh.dx = dx;
    mesh.dy = dy;
    mesh.dz = dz;
    mesh.nx = nx;
    mesh.ny = ny;
    mesh.nz = nz;

    //total_records
    uobj.get(cursor.sub("total_records"), totalRecords);

    //output_interval_steps
    uobj.get(cursor.sub("output_interval_steps"), interval);

    //output_data_flag
    string vstr,sstr,dtstr;
    uobj.get(cursor.sub("output_data_flag.velocity_field"), vstr);
    uobj.get(cursor.sub("output_data_flag.stress_field"), sstr);
    uobj.get(cursor.sub("variable_dt"), dtstr);
    if (vstr == "yes") {
        vector_flag = true;
    }
    if (sstr == "yes") {
        stress_flag = true;
    }
    if (dtstr == "yes") {
        variable_dt_flag = true;
    }

    string head  = "========== Simulation Conditions ==========\n";
    string  mid  = "Mesh(nx,ny,nz) ........... (" + itos(mesh.nx) + "," + itos(mesh.ny) + "," + itos(mesh.nz) + ")\n";
            mid += "Mesh(dx,dy,dz) ........... (" + dtos(mesh.dx) + "," + dtos(mesh.dy) + "," + dtos(mesh.dz) + ")\n";
            mid += "Total Records ............ " + itos(totalRecords) + "\n";
            mid += "Output Interaval Steps ... " + itos(interval) + "\n";
            mid += "Variable dt .............. " + dtstr + "\n";
            mid += "Output Vector Fields ..... " + vstr + "\n";
            mid += "Output Stress Fields ..... " + sstr + "\n";
            logger( head + mid );
    return;
}

void DataIO::readDynamics(UDFManager & uobj,
                            Dynamics* & dynamics,
                            vector<Component> & comps) {
    Location    cursor("dynamics");
    string      head = "========== Dynamics ==========\n";
    string      model, mid;
    uobj.get(cursor.sub("model"), model);

    if (model == "time_depend_Ginzburg_Landau") { 
        uobj.get(cursor.sub("time_depend_Ginzburg_Landau.delta_t"), dt);
        dynamics = new TDGL(mesh, dt, M);
        stress_flag = false;
        mid =  "Dynamics Type ... " + dynamics->getName() + "\n";
        mid += "Delta t ......... " + dtos(dt) + "\n";
    } else if(model == "visco_elastic"){
        Location    vcur = cursor.sub("visco_elastic");
        double      Ne, phic, converge, omega;
        uobj.get(vcur.sub("delta_t"), dt);
        uobj.get(vcur.sub("Ne"), Ne);
        uobj.get(vcur.sub("converge"), converge);
        uobj.get(vcur.sub("critical_phi"), phic);
        uobj.get(vcur.sub("omega_SOR"), omega);
        uobj.get(vcur.sub("model"), model);

        dynamics = new ViscoElastic(mesh, dt, M, Ne, phic, converge, comps, model, omega);
        for (int i = 0; i < M; i++) {
            comps[i].set_mpzeta(Ne);
        }

        mid =  "Dynamics Type ... " + dynamics->getName() + "\n";
        mid += "Delta t ......... " + dtos(dt) + "\n";
        mid += "Ne .............. " + dtos(Ne) + "\n";
        mid += "Critical Phi .... " + dtos(phic) + "\n";
        mid += "Solver Method ... " + model + "\n";
        if (model == "MAC-SOR") {
            mid += "Omega SOR ....... " + dtos(omega) + "\n";
        }
        mid += "Convergence ..... " + dtos(converge) + "\n";


    //} else if(model =="diffusion_equation"){
    //    uobj.get(cursor.sub("diffusion_equation.delta_t"), dt);
    //    dynamics = new Diffusion(mesh, dt, M);
    //    mid = "Dynamics Type ........... " + dynamics->get_name() + "\n";
    } else {
        warning("readDynamics:"+model);
        exit(1);
    }

    logger( head + mid );
    return;
}

void DataIO::readReactions(UDFManager & uobj,
                            Reaction & reaction) {
    Location    cursor("reaction.reaction[]");
    int         num_reactions = uobj.size(cursor);

    if(num_reactions>0){
        int           rid1, rid2, pid, rinter, size;
        double        k;
        string        model;

        size = mesh.nx*mesh.ny*mesh.nz;

        uobj.get(Location("reaction.interval"), rinter);
        rinter = 1;
        if (rinter <= 0 || rinter >= interval) {
            warning("invalid reaction interval" + itos(rinter));
        }
        reaction.set(num_reactions, rinter);

        for (int i = 0; i < num_reactions; i++) {
            cursor.next();
            uobj.get(cursor.sub("model"), model);
            Location cur(cursor.sub(model));
            if (model == "first_order") {
                uobj.get(cur.sub("reactant"), rid1);
                uobj.get(cur.sub("product"), pid);
                uobj.get(cur.sub("reaction_rate"), k);
                reaction.reactions[i] = new FirstOrder(rid1, pid, k, rinter, size);
            } else if (model == "second_order") {
                uobj.get(cur.sub("reactant1"), rid1);
                uobj.get(cur.sub("reactant2"), rid2);
                uobj.get(cur.sub("product"), pid);
                uobj.get(cur.sub("reaction_rate"), k);
                reaction.reactions[i] = new SecondOrder(rid1, rid2, pid, k, rinter, size);
            } else if (model == "chain_growth") {
                uobj.get(cur.sub("reactant"), rid1);
                uobj.get(cur.sub("product"), pid);
                uobj.get(cur.sub("reaction_rate"), k);
                reaction.reactions[i] = new ChainGrowth(rid1, pid, k, rinter, size);
            } else {
                warning("readReactions:" + model);
            }
        }
    }
    logger( reaction.output() );
    return;
}


void DataIO::write(UDFManager & uobj,
                   const int record,
                   const double free_energy,
                   vector<Component> & comps,
                   vector<ScalarField> & phi,
                   Dynamics* & dynamics,
                   const string & lap_time,
                   const string & elapse_time,
                   const string & remaining_time,
                   const string & comment) {
    if (record==totalRecords) {
        vector_flag = true;
        stress_flag = true;
        if(dynamics->getName()!="Visco Elastic")
            stress_flag = false;
    }
    uobj.newRecord( comment );
    Location    out( "output" );
    uobj.put( out.sub("steps"), record + init_step );
    uobj.put( out.sub("time"), record*interval*dt + init_time );
    uobj.put( out.sub("free_energy"), free_energy );

    Location    cursor = out.sub("components[]");
    for (int i = 0; i < M; i++) {
        cursor.next();
        uobj.put( cursor.sub("name"), comps[i].get_name() );
        uobj.put( cursor.sub("volume_fraction"), comps[i].vf );
        putField( uobj, cursor.sub("density[]"), phi[i] );
        if (vector_flag) {
            putVectorField( uobj, cursor.sub("velocity[]"), dynamics->velocity[i] );
        }
        if (stress_flag) {
            putField(uobj, cursor.sub("bulk_stress[]"), dynamics->bulk[i]);
            putStressField(uobj, cursor.sub("shear_stress[]"), dynamics->shear[i]);
        }
    }

    logger(align(progress(record, totalRecords), 9)  + "  " +
           align(lap_time, 11)                       + "  " +
           align(elapse_time, 18)                    + "  " +
           align(remaining_time, 18)                 + "  " +
           align(dtos(free_energy), 14)                      );

    uobj.write();
    return;
}

template <class T>
void DataIO::putField(UDFManager & uobj,
                      Location cursor,
                      const T & f) {
    for (int i = 0; i < f.size(); i++) {
        cursor.next();
        uobj.put(cursor, f[i]);
    }
    return;
}

template <class T>
bool DataIO::readField(UDFManager & uobj,
                       Location cursor,
                       T & f) {
    for (int i = 0; i < f.size(); i++) {
        cursor.next();
        uobj.get(cursor, f[i]);
    }
    return true;
}

bool DataIO::readVectorField(UDFManager & uobj,
                             Location cursor,
                             VectorField & vf) {
    for (int i = 0; i < vf.size(); i++) {
        cursor.next();
        uobj.get(cursor.sub("x"), vf[i].x);
        uobj.get(cursor.sub("y"), vf[i].y);
        uobj.get(cursor.sub("z"), vf[i].z);
    }
    return true;
}
void DataIO::putVectorField(UDFManager & uobj,
                            Location cursor,
                            const VectorField & vf) {
    for (int i = 0; i < vf.size(); i++) {
        cursor.next();
        uobj.put(cursor.sub("x"), vf[i].x);
        uobj.put(cursor.sub("y"), vf[i].y);
        uobj.put(cursor.sub("z"), vf[i].z);
    }
    return;
}
bool DataIO::readStressField(UDFManager & uobj,
                             Location cursor,
                             StressField & sf) {
    for (int i = 0; i < sf.size(); i++) {
        cursor.next();
        uobj.get(cursor.sub("xx"), sf[i].xx);
        uobj.get(cursor.sub("xy"), sf[i].xy);
        uobj.get(cursor.sub("xz"), sf[i].xz);
        //uobj.get(cursor.sub("yx"), sf[i].yx);
        uobj.get(cursor.sub("yy"), sf[i].yy);
        uobj.get(cursor.sub("yz"), sf[i].yz);
        //uobj.get(cursor.sub("zx"), sf[i].zx);
        //uobj.get(cursor.sub("zy"), sf[i].zy);
        uobj.get(cursor.sub("zz"), sf[i].zz);
    }
    return true;
}
void DataIO::putStressField(UDFManager & uobj,
                            Location cursor,
                            const StressField & sf) {
    for (int i = 0; i < sf.size(); i++) {
        cursor.next();
        uobj.put(cursor.sub("xx"), sf[i].xx);
        uobj.put(cursor.sub("xy"), sf[i].xy);
        uobj.put(cursor.sub("xz"), sf[i].xz);
        //uobj.put(cursor.sub("yx"), sf[i].yx);
        uobj.put(cursor.sub("yy"), sf[i].yy);
        uobj.put(cursor.sub("yz"), sf[i].yz);
        //uobj.put(cursor.sub("zx"), sf[i].zx);
        //uobj.put(cursor.sub("zy"), sf[i].zy);
        uobj.put(cursor.sub("zz"), sf[i].zz);
    }
    return;
}
