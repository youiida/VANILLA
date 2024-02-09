#ifndef _INC_MY_UDF_IO_H_
#define _INC_MY_UDF_IO_H_

#include "MyTools.h"
#include "udfmanager.h"

void simplifyInput(string & in_name,
                   string & out_name,
                   string & log_name,
                   bool key_out,
                   bool key_log) {
    if (!key_out) {
        out_name = replace2(in_name, "_in.udf", ".udf", "_out.bdf");
        if (out_name.size()<1) {
            out_name = replace(in_name, ".", "_out.bdf");
            log_name = replace(in_name, ".", "_log.dat");
        }
    }
    if (!key_log) {
        log_name = replace2(in_name, "_in.udf", ".udf", "_log.dat");
        if (log_name.size()<1) {
            log_name = replace(in_name, ".", "_log.dat");
        }
    }
    return;
}

void usage() {
    cout << ">> (exe) -I (input_udf) -O (output_udf) -L (log_file) -N (openMP threads)" << endl;
    cout << "or" << endl;
    cout << ">> (exe) (input_udf) [auto output_udf, log_file name]" << endl;
    cout << "options (-I, -O, -L, -N) are exchangeable" << endl;

    exit(2);
    return;
}
void parseArg(const int argc,
              char* argv[],
              string & input_udf,
              string & output_udf,
              string & log_file,
              int & num_threads) {
    string    arg;
    bool      key_in = false;
    bool      key_out = false;
    bool      key_log = false;
    char      at;

    try {
        for (int i = 1; i<argc; i++) {
            arg = string(argv[i]);
            if (arg.at(0) == '-') {
                at = arg.at(1);
                if (at == 'I') {
                    key_in = true;
                    input_udf = string(argv[i + 1]);
                    i++;
                }
                else if (at == 'O') {
                    key_out = true;
                    output_udf = string(argv[i + 1]);
                    i++;
                }
                else if (at == 'L') {
                    key_log = true;
                    log_file = string(argv[i + 1]);
                    i++;
                }
                else if (at == 'N') {
                    num_threads = atoi(argv[i + 1]);
                    i++;
                }
                else {
                    usage();
                }
            }
        }
        if (!key_in && argc > 1) {
            input_udf = string(argv[1]);
            key_in = true;
        }
    } catch (out_of_range& e) {
        cerr << "Do not put space or tab just after '-'" << endl;
        cerr << e.what() << endl;
        usage();
    }
    if (key_in) {
        simplifyInput(input_udf, output_udf, log_file, key_out, key_log);
    } else {
        usage();
    }

    try {
        UDFManager    in_uobj(input_udf);
    } catch (UDFManager::UDFManagerException & e) {
        cerr << e.what() << endl;
        cerr << "Input udf and/or def_udf does not exsist!" << endl;
        exit(2);
    }
    try {
        UDFManager    out_uobj(output_udf, input_udf);

    } catch (UDFManager::UDFManagerException & e) {
        cerr << e.what() << endl;
        cerr << "Can not set output udf file name! " << endl;
        cerr << "(invalid name or input udf is not '.udf' ?)" << endl;
        exit(2);
    }
    return;
}



#endif
