#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include "schedule_models.hpp"
extern "C" void pcctp_state_features(int,int,const int*,const unsigned char*,const double*,double,
 const double*,int,const unsigned char*,double,double,double*);
extern "C" void pcctp_operator_scores(const double*,const double*,const double*,const double*,const double*,int,double*,double*);
static double scheduler_rewards[15]={},scheduler_times[15]={};
static int scheduler_counts[15]={};
static void schedule_update(int move,double gain,double seconds) {
  scheduler_counts[move]++;scheduler_rewards[move]+=(gain>0);
  scheduler_times[move]+=seconds;
}
static void schedule_order(Vertice* vertices,Matriz* C,int* successors,int quota,int raw_n,int* order) {
  const char* policy=std::getenv("PCCTP_OPERATOR_POLICY");
  if(!policy || !std::strcmp(policy,"original"))return;
  // Simple controls avoid state features they do not use.
  if(!std::strcmp(policy,"static") || !std::strcmp(policy,"adaptive") || !std::strcmp(policy,"size")) {
    double x[10]={},p[15],lt[15],score[15];int mode=!std::strcmp(policy,"size")?1:0;
    if(mode==1) {
      int visits=0,length=0;
      for(int i=0;i<raw_n;++i){visits+=std::strcmp(vertices[i].tipo,"w")!=0;length+=successors[i]>=0;}
      x[0]=std::log1p(visits);x[1]=std::log1p(length);
    }
    pcctp_operator_scores(x,model_mean[mode],model_scale[mode],model_probability[mode],model_timing[mode],mode,p,lt);
    for(int m=0;m<15;++m)score[m]=!std::strcmp(policy,"adaptive")?
      (scheduler_rewards[m]+1.)/(scheduler_times[m]+.001):p[m]/std::exp(std::max(-20.,std::min(20.,lt[m])));
    std::stable_sort(order,order+15,[&](int a,int b){return score[a]>score[b];});return;
  }
  static bool initialized=false;static int n,demands;
  static double mean,total;static std::vector<int> map,ids;
  static std::vector<double> cost,prizes;
  static std::vector<unsigned char> mandatory,cover;
  if(!initialized) {
    map.assign(raw_n,-1);demands=0;total=0.;mean=0.;
    for(int v=0;v<raw_n;++v)if(std::strcmp(vertices[v].tipo,"w")) {map[v]=ids.size();ids.push_back(v);}else demands++;
    n=ids.size();cost.resize(n*n);prizes.resize(n);mandatory.resize(n);cover.assign(demands*n,0);
    for(int i=0;i<n;++i) {
      mandatory[i]=!std::strcmp(vertices[ids[i]].tipo,"t");prizes[i]=vertices[ids[i]].premio;total+=prizes[i];
      for(int j=0;j<n;++j){cost[i*n+j]=get(C,ids[i],ids[j]);mean+=cost[i*n+j];}
    }
    mean/=n*n;int w=0;
    for(int v=0;v<raw_n;++v)if(!std::strcmp(vertices[v].tipo,"w")) {
      for(int j=0;j<vertices[v].qtd_coberto_por;++j) {
        int index=map[vertices[v].coberto_por[j]];if(index>=0)cover[w*n+index]=1;
      }w++;
    }initialized=true;
  }
  int root=0;while(root<raw_n && successors[root]<0)root++;
  std::vector<int> route;int v=root;
  do {route.push_back(map[v]);v=successors[v];}while(v!=root);
  double x[10],p[15],logtime[15],unused[15],score[15];
  pcctp_state_features(n,route.size(),route.data(),mandatory.data(),prizes.data(),quota,cost.data(),
    !std::strcmp(policy,"generic")?0:demands,cover.data(),mean,total,x);
  const char* verify=std::getenv("PCCTP_SCHEDULE_VERIFY");static int exported=0;
  if(verify && exported++<20) {
    std::ofstream out(verify,std::ios::app);out<<std::setprecision(17)<<"{\"route\":[";
    for(size_t i=0;i<route.size();++i)out<<(i?",":"")<<ids[route[i]];
    out<<"],\"features\":[";
    for(int i=0;i<10;++i)out<<(i?",":"")<<x[i];out<<"]}\n";
  }
  int mode=!std::strcmp(policy,"static")?0:!std::strcmp(policy,"size")?1:!std::strcmp(policy,"generic")?2:3;
  pcctp_operator_scores(x,model_mean[mode],model_scale[mode],model_probability[mode],model_timing[mode],mode,p,logtime);
  // Generic and full share the same generic runtime prediction.
  if(mode==3)pcctp_operator_scores(x,model_mean[2],model_scale[2],model_probability[2],model_timing[2],2,unused,logtime);
  for(int m=0;m<15;++m) {
    score[m]=p[m]/std::exp(std::max(-20.,std::min(20.,logtime[m])));
    if(!std::strcmp(policy,"adaptive"))
      score[m]=(scheduler_rewards[m]+1.)/(scheduler_times[m]+.001);
  }
  std::stable_sort(order,order+15,[&](int a,int b){return score[a]>score[b];});
}
