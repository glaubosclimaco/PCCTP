#include <cmath>
#include <vector>
#include <algorithm>
// Costs/coverage refer only to visit vertices; static summaries cached by caller.
extern "C" void pcctp_state_features(int n,int length,const int* route,
 const unsigned char* mandatory,const double* prizes,double quota,const double* cost,
 int demands,const unsigned char* cover,double mean_cost,double total_prize,double* out) {
  double edge_sum=0.,squares=0.,prize=0.,saving=0.;int required=0,tight=0;
  for(int v=0;v<n;++v)required+=mandatory[v]!=0;
  for(int i=0;i<length;++i) {
    int v=route[i],a=route[(i+length-1)%length],b=route[(i+1)%length];
    double edge=cost[a*n+v];edge_sum+=edge;squares+=edge*edge;prize+=prizes[v];
    if(!mandatory[v])saving+=std::max(0.,cost[a*n+v]+cost[v*n+b]-cost[a*n+b]);
  }
  for(int w=0;w<demands;++w) {
    int count=0;
    for(int i=0;i<length;++i)count+=cover[w*n+route[i]]!=0;
    tight+=count==1;
  }
  double average=edge_sum/length;
  out[0]=std::log1p(n);out[1]=std::log1p(length);out[2]=double(length)/n;
  out[3]=double(required)/n;out[4]=quota/std::max(1.,total_prize);
  out[5]=prize/std::max(1.,quota);out[6]=double(tight)/std::max(1,demands);
  out[7]=average/std::max(1.,mean_cost);
  out[8]=std::sqrt(std::max(0.,squares/length-average*average))/std::max(1.,average);
  out[9]=saving/std::max(1.,edge_sum);
}
