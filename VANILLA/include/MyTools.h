#ifndef _INC_MY_TOOLS_H_
#define _INC_MY_TOOLS_H_
#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <numeric>
#include <complex>

#pragma warning(disable : 4996)  // localtime,asctime
#pragma warning(disable : 4290)  // nothrow
using namespace std;

#ifndef _DCOMPLEX_
#define _DCOMPLEX_
typedef complex<double> dcomplex;
#endif
#if defined(_CONSOLE)
#define finite(x) _finite(x)
#define isnan(x) _isnan(x)
#else
#define INT_MAX std::numeric_limits<int>::max()
#endif

#ifndef PI
#define PI  (3.14159265358979)
#define PI2 (6.28318530717959)
#endif

int getOrder(const int integer);

/* string handle */
list<string>    split(string str, const string & separator);
string          getHead(const string & str, const string & separator);
string          getTail(const string & str, const string & separator);
string          eraseHead(const string & str, const string & separator);
string          eraseTail(const string & str, const string & separator);
string          sfill(const int integer, const int max);
string          zfill(const int integer, const int max);
string          cfill(const int integer, const int max, const string & c);
string          cfill(string & text, const string & c, const int size);
string          itobyte(const int byte);
string          itobyte(const int byte, const string & unit);
string          replace(const string & target_str, const string & tar_word, const string & rep_word);
string          replace2(const string & target_str, const string & tar1,     const string & tar2, const string & rep);
string          align(const string & str, const int order);
string          progress(const int present, const int end);
/* read file */
void    readtext(vector< int            > & data, fstream & fobj, const char sep);
void    readtext(vector< double         > & data, fstream & fobj, const char sep);
void    readtext(vector< string         > & data, fstream & fobj, const char sep);
void    readtext(vector< vector<int>    > & data, fstream & fobj, const char sep);
void    readtext(vector< vector<double> > & data, fstream & fobj, const char sep);
void    readtext(vector< vector<string> > & data, fstream & fobj, const char sep);

/* convert type */
string    dtos(const double x);
string    btos(const bool b);
string    itos(const int i);

bool      stob(const string & str);

#ifndef VCC
double    stod(const string & str);
int       stoi(const string & str);
#endif

/* warning caution */
void    warning(const string & msg);
void    caution(const string & msg);

/* error */
template<class T>
double  error(const T before, const T after);
bool    isError(const double value);
bool    isInRange(const double data, const double upper, const double lower);
bool    isShouldBe(const double data, const double ideal, const double tolerance);

/* math */
double  gaussWindow(const double x, const double x0);
//  template<class T>
//  double    step(const T var1, const T var2);
template<class T>
double step(const T var1, const T var2) {
    if(var1 > var2) {
        return 1.0;
    } else {
        return 0.0;
    }
}
double ddelta(int i1, int i2);

#include "Vector3d.h"
Vector3d vdelta(int i1, int i2);
#ifdef _INC_VECTOR3D_H_
bool setInvMatrix(vector<Vector3d> & target, vector<Vector3d> & inv, vector<double> & x, vector<int> & indx, int size);
#endif
bool setInvMatrix(vector<double>   & target, vector<double>   & inv, vector<double> & x, vector<int> & indx, double & det, int size);

void lubksb(vector<double> & a, vector<int> & indx, vector<double> & b, int n);
bool ludcmp(vector<double> & a, vector<int> & indx, vector<double> & vv, double & d, int n);

#ifndef _INC_OPT_HISTOGRAM_H_
#define _INC_OPT_HISTOGRAM_H_
/* opt_histogram */
class OptHistogram {
private:
    double    width, max, min;
    int       bin;

    vector<int> get_hist(vector<double> & data,
                         double max,
                         double min,
                         int bins) {
        vector<int>    temp_hist;
        int            ir;
        temp_hist.resize(bins, 0);
        double delta = (max - min) / (static_cast<double>(bins));
        for (int j = 0; j < data.size(); j++) {
            ir = static_cast<int>((data[j] - min) / delta);
            if (ir == (bins))ir -= 1;
            temp_hist[ir] += 1;
        }
        return temp_hist;
    }

public:
    OptHistogram() {}
    ~OptHistogram() {}
    void make_hist(vector<double> & data) {
        int min_bin = 2;
        int max_bin = static_cast<int>(data.size()) / 2;
        sort(data.begin(), data.end());
        max = data.back();
        min = data.front();

        vector<double>    D;  // bin sizes
        vector<double>    C;  // cost parameter

        for (int i = min_bin; i <= max_bin; i++) {
            D.push_back((max - min) / static_cast<double>(i));
            C.push_back(0.0);
        }

        double sum, sum2, mean, dev, delta;

        for (int i = min_bin; i <= max_bin; i++) {
            hist = get_hist(data, max, min, i);

            sum = 0.0;
            sum2 = 0.0;
            for (int j = 0; j < static_cast<int>(hist.size()); j++) {
                sum += hist[j];
                sum2 += hist[j] * hist[j];
            }
            mean = sum / static_cast<double>(i);
            dev = sum2 / static_cast<double>(i) - mean*mean;
            C[i - min_bin] = (2.0*mean - dev) / (delta*delta);
        }

        double temp = C[0];
        int idx = 0;
        for (int i = 1; i < static_cast<int>(C.size()); i++) {
            if (temp > C[i]) {
                temp = C[i];
                idx = i;
            }
        }

        hist = get_hist(data, max, min, idx);
        hist_val.resize(hist.size());
        for (int i = 0; i < static_cast<int>(hist_val.size()); i++) {
            hist_val[i] = min + (0.5 + i)*delta;
        }
    }

    vector<int>       hist;
    vector<double>    hist_val;  // [min+0.5*delta, min+1.5*delta, ... max-0.5*delta]
    void output(string name,
                string norm = "yes") {
        fstream out;
        out.open(name.c_str(), ios::out);
        string type = getTail(name, ".");
        string sep;
        if (type == "csv") {
            sep = ",";
        } else {
            sep = "\t";
        }
        double count = accumulate(hist.begin(), hist.end(), 0.0);


        out << "value" + sep + "population\n";
        for (int i = 0; i < static_cast<int>(hist.size()); i++) {
            out << hist_val[i] << sep;
            if (norm == "yes") {
                out << static_cast<double>(hist[i]) / count << endl;
            }
        }
        out.close();
        return;
    }
};
#endif

#ifndef _INC_TIMER_H_
#define _INC_TIMER_H_

class Timer{
private:
    clock_t    timer_start, prev;
	string itotime(int clock) {
		int hour = clock / 3600;
		if (hour >= 24) {
			int day = hour / 24;
			hour -= day * 24;
			return sfill(day, 10) + "[d] " + zfill(hour, 24) + "[h]";
		}
		int       sec = clock % 60;
		int       min = (clock % 3600) / 60;
		string    time;

		if (hour != 0) {
			time = sfill(hour, 100) + "[h] " + zfill(min, 60) + "[m] " + zfill(sec, 60) + "[s]";
		} else {
			if (min != 0) {
				time = sfill(min, 60) + "[m] " + zfill(sec, 60) + "[s]";
			} else {
				time = sfill(sec, 60) + "[s]";
			}
        }
        return time;
    }
    string getDate(void) {
        time_t      rawtime;
        time(&rawtime);
        struct tm * timeinfo = localtime(&rawtime);
        string      str      = asctime(timeinfo);
        return      str.substr(0, str.size()-1);
    }
public:
    Timer() {timer_start =  clock(); prev = clock(); }
    ~Timer() {}
    string lapTime(const clock_t present) {
        int clocks = static_cast<int>(present - prev) / static_cast<int>CLOCKS_PER_SEC;
        prev = present;
        return itotime(clocks);
    }
    string elapseTime(const clock_t present) {
        int clocks = static_cast<int>(present - timer_start) / (static_cast<int>CLOCKS_PER_SEC;
        return itotime(clocks);
    }
    string estTime(const clock_t present, 
                   const int record, 
                   const int total) {
        int clocks = static_cast<int>(present - timer_start) / static_cast<int>CLOCKS_PER_SEC;
            clocks *= (total - record);
        record !=0 ? clocks /= record : clocks = 0;
        return itotime(clocks);
    }
    string finish(void){
        return "The program executed successfully at " + getDate() + "\nexecution time:" + elapseTime(clock());
    }
};
#endif

#ifndef    _INC_RANDOM_H_
#define    _INC_RANDOM_H_

#define NTAB 32
class UniformRandomNumber {
  public:
     UniformRandomNumber();
     float operator() ();
     UniformRandomNumber&  setMax(float max_ = 1.0);
     UniformRandomNumber&  setMin(float min_ = 0.0);
     UniformRandomNumber&  setSeed(int idum_ = 0);
  private:
     float max;
     float min;
     float a, b;
     int idum;
     int iy;
     int iv[NTAB];
     static int count;
};

class GaussianRandomNumber {
public:
    GaussianRandomNumber();
    float operator() ();
    GaussianRandomNumber&  setSigma(float sigma_ = 1.0);
    GaussianRandomNumber&  setMean(float mean_ = 0.0);
    GaussianRandomNumber&  setSeed(int idum_ = 0);
private:
    float sigma;
    float mean;
    UniformRandomNumber rand;
    int iset;
    float gset;
};
#endif

#endif
