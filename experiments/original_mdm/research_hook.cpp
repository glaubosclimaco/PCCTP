#include "research_hook.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
int n,seed,support,stability,calls=0;
int observed_iterations=0;
bool final_export=false;
std::string instance,export_path;
std::chrono::steady_clock::time_point began;
struct Snapshot { int iteration; double cost; std::vector<int> route; };
std::vector<Snapshot> snapshots;
std::vector<int> decode(int* successors) {
  int root=0; while(root<n && successors[root]<0) ++root;
  if(root==n) throw std::runtime_error("Empty legacy route");
  int v=root; std::vector<int> route; std::set<int> seen;
  do {
    if(v<0 || v>=n || seen.count(v)) throw std::runtime_error("Malformed legacy cycle");
    route.push_back(v); seen.insert(v); v=successors[v];
  } while(v!=root);
  for(int i=0;i<n;i++) if(successors[i]>=0 && !seen.count(i)) throw std::runtime_error("Multiple legacy cycles");
  return route;
}
void write_route(std::ostream& out,const std::vector<int>& route) {
  out << "["; for(size_t j=0;j<route.size();j++) out << (j ? ",":"") << route[j]; out << "]";
}
}
void research_initialize(const char* path,int size,int s,int sup,int stable,const char* output) {
  instance=path; n=size; seed=s; support=sup; stability=stable;
  export_path=output; observed_iterations=0; final_export=false;
  calls=0; snapshots.clear(); began=std::chrono::steady_clock::now();
}
void research_snapshot(int iteration,double value,int* route) {
  observed_iterations=iteration;
  if(iteration==1 || iteration==5 || iteration==10 || iteration==20 || iteration==35 || iteration==70) {
    snapshots.push_back({iteration,value,decode(route)});
    std::string temporary=export_path+".checkpoint.tmp";
    research_export(temporary.c_str(),value,route);
    if(std::rename(temporary.c_str(),(export_path+".checkpoint").c_str())!=0)
      throw std::runtime_error("Cannot publish legacy checkpoint");
  }
}
void research_mining() { ++calls; }
void research_export(const char* path,double value,int* route) {
  if(export_path==path) final_export=true;
  std::ofstream out(path);
  if(!out) throw std::runtime_error("Cannot write legacy route export");
  out << std::setprecision(16) << "{\"method\":\"original_MDM_legacy_instrumented\",\"seed\":" << seed
      << ",\"support_absolute\":" << support << ",\"stability_percent\":" << stability
      << ",\"mining_calls\":" << calls << ",\"iterations\":" << observed_iterations
      << ",\"completed\":" << (final_export ? "true":"false") << ",\"cost\":" << value
      << ",\"seconds\":" << std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()
      << ",\"route\":";
  write_route(out,decode(route)); out << ",\"snapshots\":[";
  for(size_t j=0;j<snapshots.size();j++) {
    auto& s=snapshots[j]; out << (j ? ",":"") << "{\"iteration\":" << s.iteration << ",\"cost\":" << s.cost << ",\"route\":";
    write_route(out,s.route); out << "}";
  }
  out << "]}\n";
}
