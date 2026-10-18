// Portable inference only. Feature extraction and search are separate costs.
#include <cmath>
#include <algorithm>
extern "C" void pcctp_operator_scores(const double* x,const double* mean,const double* scale,
 const double* probability,const double* timing,int context,double* p,double* logtime) {
  const int features=context==0?0:context==1?2:10;
  double z[10];
  for(int j=0;j<features;++j)z[j]=std::max(-5.,std::min(5.,(x[j]-mean[j])/scale[j]));
  for(int move=0;move<15;++move) {
    double a=probability[move],b=timing[move];
    for(int j=0;j<features;++j) {
      a+=z[j]*probability[15+move*features+j];
      b+=z[j]*timing[15+move*features+j];
    }
    p[move]=a>=0?1./(1.+std::exp(-a)):std::exp(a)/(1.+std::exp(a));
    logtime[move]=b;
  }
}
