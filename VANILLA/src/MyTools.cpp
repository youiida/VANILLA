#include "MyTools.h"

int getOrder(const int integer){
    stringstream stream;
    int absInt = abs(integer);
    stream << absInt;
    return (int)(stream.str()).size();
}

list<string> split(string str, 
                   const string & separator){
    list<string> result;
    int cutAt;
    while( (cutAt = (int)str.find_first_of(separator)) != str.npos ){
        if(cutAt > 0)
            result.push_back(str.substr(0, cutAt));
        str = str.substr(cutAt + 1);
    }
    if(str.length() > 0)
        result.push_back(str);
    return result;
}

string replace(const string & target_str, 
               const string & tar_word, 
               const string & rep_word){
    string    rs;
    int        idx    =    (int)target_str.find(tar_word);
    if(idx!=-1){
        rs    =    target_str.substr(0,idx) + rep_word + target_str.substr(idx+tar_word.size());
        rs    =    replace(rs, tar_word, rep_word);
    }else{
        rs    =    target_str;
    }
    return rs;
}

string    replace2(const string & target_str,
                   const string & tar1,
                   const string & tar2,
                   const string & rep) {
    string rs = replace(target_str,
        tar1,
        rep);
    if (rs == target_str)
        rs = replace(target_str,
            tar2,
            rep);
    return rs;
}

string getHead(const string & str, 
               const string & separator){
    int    idx    =    (int)str.find(separator);
    return str.substr(0,idx);
}

string getTail(const string & str, 
               const string & separator){
    int    idx    =    (int)str.rfind(separator);
    return str.substr(idx+separator.size());
}

string eraseTail(const string & str, 
                 const string & separator){
    int    idx    =    (int)str.rfind(separator);
    return str.substr(0,idx);
}

string eraseHead(const string & str, 
                 const string & separator){
    int    idx    =    (int)str.find(separator);
    return str.substr(idx);    
}

string itobyte(const int byte){
    string unit;

    if(byte>1024){
        int    KB        = (int)(double(byte) / 1024.0);
        int    perKB     = (int)(100.0*double((byte % 1024))/1024.0);
        if(KB>1024.0){
            int    MB       = (int)(KB / 1024.0);
            int    perMB    = (int)(100.0*double((KB % 1024))/1024.0);
            if(MB>1024.0){
                int    GB       = (int)(MB / 1024.0);
                int    perGB    = (int)(100.0*double((MB % 1024))/1024.0);
                unit = zfill(GB, 1) + "." + zfill(perGB, 1) + " GB";
            }else
                unit = zfill(MB, 1) + "." + zfill(perMB, 1) + " MB";
        }else
            unit = zfill(KB, 1) + "." + zfill(perKB, 1) + " KB";
    }else
        unit = zfill(byte,1) + " B";
    return unit;
}

string itobyte(const int byte, 
               const string & prefix){
    string unit;
    if(prefix=="K"){
        if(byte>1024.0){
            int    MB       = (int)(byte / 1024.0);
            int    perMB    = (int)(100.0*double((byte % 1024))/1024.0);
            if(MB>1024.0){
                int    GB       = (int)(MB / 1024.0);
                int    perGB    = (int)(100.0*double((MB % 1024))/1024.0);
                unit = zfill(GB, 1) + "." + zfill(perGB, 1) + " GB";
            }else
                unit = zfill(MB, 1) + "." + zfill(perMB, 1) + " MB";
        }else
            unit = zfill(byte, 1) + " KB";
    }else if(prefix=="M"){
            if(byte>1024.0){
                int    GB       = (int)(byte / 1024.0);
                int    perGB    = (int)(100.0*double((byte % 1024))/1024.0);
                unit = zfill(GB, 1) + "." + zfill(perGB, 1) + " GB";
            }else
                unit = zfill(byte, 1) + " MB";
    }else{
        warning("invalid itobypte prefix");
    }
    return unit;
}


string zfill(const int integer, 
             const int max){
    return cfill(integer, max, "0");
}

string sfill(const int integer, 
             const int max){
    return cfill(integer, max, " ");
}
string progress(const int present, 
                const int end){
    return sfill(present,end)+"/"+sfill(end,end);
}

string cfill(const int integer, 
             const int max, 
             const string & c){
    stringstream stream;
    string       str;
    int          maxOrder = getOrder(max);
    int          absInt = abs(integer);
    stream << absInt;
    int          order = (int)(stream.str()).size();
    int          i;
    for(i=0;i<(maxOrder - order);i++)
        str += c;
    if(integer>=0)    return( str + stream.str() );
    else              return( "-" + str + stream.str() );
}

string cfill(string & text, 
             const string & c, 
             const int size){
    int textSize = (int)text.size();
    int i;
    for(i=0;i<(size-textSize);i++)    text += c;
    return text;
}


string align(const string & str, 
             const int order){
    int        size = (int)str.size();
    string     rs;
    int        i;
    for(i=0;i<(order-size);i++)    rs += " ";
    return rs+str;
}

string dtos(const double x){
    ostringstream oss;
    oss << x;
    return oss.str();
}

string btos(const bool b){
    if(b)    return "Yes";
    else     return "No";
}

string itos(const int i){
    stringstream sst;
    sst << i;
    return sst.str();
}

bool stob(const string & str){
    if(str=="YES" || str=="Yes" || str=="yes" || str=="Y" || str=="y")    return true;
    else                                                                return false;
}

#ifndef VCC
double stod(const string &str){
    return atof(str.c_str());
}
int stoi(const string & str){
    return (int)atof(str.c_str());
}
#endif


void readtext(vector<int> & data, 
              fstream & fobj, 
              char sep){
    string line, token;
    while(getline(fobj, line)){
        istringstream stream(line);
        while( getline(stream, token, sep) ){
            data.push_back((int)atof(token.c_str()));
        }
    }
    fobj.close();
    return;
};
void readtext(vector<double> & data, 
              fstream & fobj, 
              char sep){
    string line, token;
    while(getline(fobj, line)){
        istringstream stream(line);
        while( getline(stream, token, sep) ){
            data.push_back(atof(token.c_str()));
        }
    }
    fobj.close();
    return;
};
void readtext(vector<string> & data, 
              fstream & fobj, 
              char sep){
    string line, token;
    while(getline(fobj, line)){
        istringstream stream(line);
        while( getline(stream, token, sep) ){
            data.push_back(token);
        }
    }
    fobj.close();
    return;
};
void readtext(vector< vector<int> > & data, 
              fstream & fobj, 
              char sep){
    string         line, token;
    vector<int>    temp;
    while(getline(fobj, line)){
        istringstream stream(line);
        while( getline(stream, token, sep) ){
            temp.push_back((int)atof(token.c_str()));
        }
        data.push_back(temp);
        temp.clear();
    }
    fobj.close();
    return;
};
void readtext(vector< vector<double> > & data, 
              fstream & fobj, 
              char sep){
    string            line, token;
    vector<double>    temp;
    while(getline(fobj, line)){
        istringstream stream(line);
        while( getline(stream, token, sep) ){
            temp.push_back(atof(token.c_str()));
        }
        data.push_back(temp);
        temp.clear();
    }
    fobj.close();
    return;
};
void readtext(vector< vector<string> > & data, 
              fstream & fobj, 
              char sep){
    string            line, token;
    vector<string>    temp;
    while(getline(fobj, line)){
        istringstream stream(line);
        while( getline(stream, token, sep) ){
            temp.push_back(token);
        }
        data.push_back(temp);
        temp.clear();
    }
    fobj.close();
    return;
};

void warning(const string & msg){
    cerr << endl << endl << endl;
    cerr << " !! warning !! " << msg << endl;
    //throw(msg);
    exit(1);
}

void caution(const string & msg){
    cerr << endl;
    cerr << " ## caution ## " << msg << endl;
    cerr << endl;
}

//error
template <class T>
double error(const T before, 
             const T after){
    return  fabs(before-after) / fabs(before);
}

bool isError(const double value){
    if(isnan(value) || !isfinite(value) )    return true;
    else                                    return false;
}

bool isInRange(const double data, 
               const double upper, 
               const double lower){
    if( data < lower || data > upper )    return false;
    else                                return true;
}
bool isShouldBe(const double data, 
                const double ideal, 
                const double tolerance){
    if( data < (ideal-tolerance) || data > (ideal+tolerance) )    return false;
    else                                                        return true;
}

//math
double gaussWindow(const double x, 
                   const double x0){
    return( 2.0 * sqrt( log(2.0) / PI ) * exp( -x * x * log(2.0) * 4.0 / (x0*x0) ) / x0 );
}

double    ddelta(int i1,
                 int i2) {
    return (i1 == i2 ? 1.0 : 0.0);
    //condition ? true : false
};


//LU-decomposition
#ifdef _INC_VECTOR3D_H_
Vector3d    vdelta(int i1,
                   int i2) {
    return (i1 == i2 ? Vector3d(1.0,1.0,1.0) : Vector3d(0.0,0.0,0.0));
    //condition ? true : false
};

bool    setInvMatrix(vector<Vector3d> & target,
                     vector<Vector3d> & inv,
                     vector<double> & x,
                     vector<int>    & indx,
                     int size) {
    vector<double>    tx(target.size());
    vector<double>    ty(target.size());
    vector<double>    tz(target.size());
    vector<double>    ix(inv.size());
    vector<double>    iy(inv.size());
    vector<double>    iz(inv.size());
    double            det;
    for (int i = 0; i < target.size(); i++) {
        tx[i] = target[i].x;
        ty[i] = target[i].y;
        tz[i] = target[i].z;
    }
    setInvMatrix(tx, ix, x, indx, det, size);
    setInvMatrix(ty, iy, x, indx, det, size);
    setInvMatrix(tz, iz, x, indx, det, size);
    for (int i = 0; i < target.size(); i++) {
        inv[i] = Vector3d(ix[i], iy[i], iz[i]);
    }
    return true;
}
#endif

bool    setInvMatrix(vector<double> & target,
                     vector<double> & inv,
                     vector<double> & x,
                     vector<int> & indx,
                     double & det,
                     int size) {
    if (!ludcmp(target, indx, x, det, size))    return false;
    int n1, n2;
    for (n2 = 0; n2<size; n2++) {
        det *= target[n2*size+n2];
        for (n1 = 0; n1<size; n1++)
            x[n1] = 0.0;
        x[n2] = 1.0;
        lubksb(target, indx, x, size);
        for (n1 = 0; n1<size; n1++)
            inv[n1*size+n2] = x[n1];
    }
    return true;
}
void    lubksb(vector<double> & a,
               vector<int> & indx,
               vector<double> & b,
               int n) {
    int       i, ii = 0, ip, j;
    double    sum;

    for (i = 0; i<n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii != 0)
            for (j = ii - 1; j<i; j++)
                sum -= a[n*i + j] * b[j];
        else if (sum != 0.0)
            ii = i + 1;
        b[i] = sum;
    }
    for (i = n - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j<n; j++)
            sum -= a[n*i + j] * b[j];
        b[i] = sum / a[n*i + i];
    }
    return;
}
bool    ludcmp(vector<double> & a,
               vector<int> & indx,
               vector<double> & vv,
               double & d,
               int n) {
    const double    TINY = 1.0e-20;
    int             i, imax, j, k;
    double          big, dum, sum, temp;
    imax = -INT_MAX;
    d    =  1.0;
    for (i = 0; i<n; i++) {
        big = 0.0;
        for (j = 0; j<n; j++)
            if ((temp = fabs(a[n*i + j])) > big)
                big = temp;
        if (big == 0.0)
            return false;
        vv[i] = 1.0 / big;
    }
    for (j = 0; j<n; j++) {
        for (i = 0; i<j; i++) {
            sum = a[n*i + j];
            for (k = 0; k<i; k++)
                sum -= a[n*i + k] * a[n*k + j];
            a[n*i + j] = sum;
        }
        big = 0.0;
        for (i = j; i<n; i++) {
            sum = a[n*i + j];
            for (k = 0; k<j; k++)
                sum -= a[n*i + k] * a[n*k + j];
            a[n*i + j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 0; k<n; k++) {
                dum = a[n*imax + k];
                a[n*imax + k] = a[n*j + k];
                a[n*j + k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[n*j + j] == 0.0)
            a[n*j + j] = TINY;
        if (j != n - 1) {
            dum = 1.0 / (a[n*j + j]);
            for (i = j + 1; i<n; i++)
                a[n*i + j] *= dum;
        }
    }
    return true;
}




#ifdef _INC_RANDOM_H_

#include <math.h>
#include <time.h>
//#include "Random.h"

/**********************************/
#define DEFAULT 0
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define NTAB 32
#define IR 2836
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

//UniformRandomNumber
UniformRandomNumber::UniformRandomNumber(){
   setSeed(0);
   max=1.0; min=0.0;
   a=1.0; b=0.0;
}

UniformRandomNumber&
UniformRandomNumber::setMax(float max_){
   max=max_; 
   a=max-min;
   return *this;
}

UniformRandomNumber&
UniformRandomNumber::setMin(float min_){
   min=min_; 
   a=max-min;
   b=min;
   return *this;
}

UniformRandomNumber&
UniformRandomNumber::setSeed(long idum_){
        int j;
        long k;
        idum = idum_;
        if (idum<0) idum=-idum;
        if (idum == 0) idum= long( time(0) ) + count++;
        for (j=NTAB+7;j>=0;j--) {
            k=idum/IQ;
            idum=IA*(idum-k*IQ)-IR*k;
            if (idum < 0) idum += IM;
            if (j < NTAB) iv[j] = idum;
        }
        iy = iv[0];
         k = idum/IQ;
        idum = IA*(idum-k*IQ)-IR*k;
        if (idum < 0) idum += IM;
        j = static_cast<int>(iy)/ static_cast<int>(NDIV);
        iy = iv[j];
        iv[j] = idum;
        return *this;
}

float 
UniformRandomNumber::operator() (){
        int j;
        long k;
        float temp;
        k = idum/IQ;
        idum = IA*(idum-k*IQ)-IR*k;
        if (idum < 0) idum += IM;
        j = static_cast<int>(iy)/ static_cast<int>(NDIV);
        iy = iv[j];
        iv[j] = idum;
        if ((temp = float(AM*iy)) > RNMX) return float(RNMX*a+b);
        else return temp*a+b;
}

GaussianRandomNumber::GaussianRandomNumber()
{  rand.setSeed(0);
   sigma=1.0; mean=0.0;
   iset=0;
}

float
GaussianRandomNumber::operator() () {
    float fac,rsq,v1,v2;
    if  (iset == 0) {
        do {
            v1 = float(2.0*rand()-1.0);
            v2 = float(2.0*rand()-1.0);
            rsq = v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=float(sqrt(-2.0*log(rsq)/rsq));
        gset = v1*fac;
        iset = 1;
        return v2*fac*sigma + mean ;
    } else {
        iset = 0;
        return gset*sigma + mean ;
    }
}

GaussianRandomNumber&  
GaussianRandomNumber::setSigma(float sigma_)
           {sigma=sigma_;return *this; }

GaussianRandomNumber&  
GaussianRandomNumber::setMean(float mean_)
           {mean=mean_;return *this; }

GaussianRandomNumber&  
GaussianRandomNumber::setSeed(int idum_)
           {rand.setSeed(idum_);return *this; }

int UniformRandomNumber::count = 0;
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#endif


