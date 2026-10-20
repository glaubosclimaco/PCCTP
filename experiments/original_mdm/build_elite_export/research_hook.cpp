#include "research_hook.hpp"
#include <chrono>
#include <time.h>
#include <algorithm>
#include <cstdlib>
#include <limits>
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
double began=0.;
struct Snapshot { int iteration; double cost; std::vector<int> route; double seconds; };
std::vector<Snapshot> snapshots;
double budget=0.,startup_debit=0.;bool early_start=false;
struct EliteCapture { std::vector<std::vector<int>> routes; std::vector<double> costs; std::string patterns; };
std::vector<EliteCapture> elite_captures;
EliteCapture pending_elite;
bool elite_pending=false;
std::vector<Snapshot> progress;
double seconds_now() {return startup_debit+research_monotonic_seconds()-began;}
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
double research_monotonic_seconds() {
  struct timespec value;
  if(clock_gettime(CLOCK_MONOTONIC,&value)!=0)throw std::runtime_error("Monotonic clock unavailable");
  return double(value.tv_sec)+double(value.tv_nsec)*1e-9;
}
void research_process_start() {
#ifdef PCCTP_TIMED
  began=research_monotonic_seconds();early_start=true;
  const char* b=std::getenv("PCCTP_BUDGET");budget=b?std::atof(b):0.;
  const char* parent=std::getenv("PCCTP_STARTED_MONOTONIC");
  if(parent)startup_debit=std::max(0.,began-std::atof(parent));
#endif
}
bool research_budget_active(){return budget>0.;}
bool research_expired(){return budget>0. && seconds_now()>=budget;}
void research_initialize(const char* path,int size,int s,int sup,int stable,const char* output) {
  instance=path; n=size; seed=s; support=sup; stability=stable;
  export_path=output; observed_iterations=0; final_export=false;
  calls=0; snapshots.clear();progress.clear();
  elite_captures.clear(); elite_pending=false;
  if(!early_start)began=research_monotonic_seconds();
}
void research_progress(double value,int* route) {
  if(!research_budget_active())return;
  double elapsed=seconds_now();
  if(elapsed>budget || (!progress.empty() && value>=progress.back().cost))return;
  progress.push_back({observed_iterations,value,decode(route),elapsed});
  std::string temporary=export_path+".checkpoint.tmp";
  research_export(temporary.c_str(),value,route);
  if(std::rename(temporary.c_str(),(export_path+".checkpoint").c_str())!=0)
    throw std::runtime_error("Cannot publish timed progress");
}
void research_snapshot(int iteration,double value,int* route) {
  observed_iterations=iteration;
  if(iteration==1 || iteration==5 || iteration==10 || iteration==20 || iteration==35 || iteration==70) {
    snapshots.push_back({iteration,value,decode(route),
      seconds_now()});
    std::string temporary=export_path+".checkpoint.tmp";
    research_export(temporary.c_str(),value,route);
    if(std::rename(temporary.c_str(),(export_path+".checkpoint").c_str())!=0)
      throw std::runtime_error("Cannot publish legacy checkpoint");
  }
}
void research_mining() { ++calls; }
static std::string json_escape(const std::string& text) {
  std::string out;
  for(unsigned char c: text) {
    if(c=='\\' || c=='"') { out.push_back('\\'); out.push_back(char(c)); }
    else if(c=='\n') out += "\\n";
    else if(c=='\r') continue;
    else out.push_back(char(c));
  }
  return out;
}
void research_elite(int rows, int vertices, const float* matrix, int columns) {
  if(vertices!=n || columns<n+2 || rows<1) throw std::runtime_error("Unexpected elite matrix");
  pending_elite.routes.clear(); pending_elite.costs.clear(); pending_elite.patterns.clear();
  for(int i=0;i<rows;i++) {
    std::vector<int> successors(n);
    for(int j=0;j<n;j++) successors[j]=static_cast<int>(matrix[(size_t)i*columns+j]);
    pending_elite.routes.push_back(decode(successors.data()));
    pending_elite.costs.push_back(matrix[(size_t)i*columns+n]);
  }
  elite_pending=true;
}
void research_capture_patterns() {
  if(!elite_pending) throw std::runtime_error("Pattern capture without an elite bank");
  std::ifstream in("padroes.txt");
  if(!in) throw std::runtime_error("Mined pattern file is missing");
  pending_elite.patterns.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  elite_captures.push_back(pending_elite);
  elite_pending=false;
  std::ofstream out((export_path+".elites.jsonl").c_str(), std::ios::app);
  if(!out) throw std::runtime_error("Cannot write elite bank");
  const EliteCapture& bank=elite_captures.back();
  out << std::setprecision(16) << "{\"mining_index\":" << (elite_captures.size()-1)
      << ",\"n\":" << n << ",\"routes\":[";
  for(size_t r=0;r<bank.routes.size();r++) {
    out << (r?",":""); write_route(out, bank.routes[r]);
  }
  out << "],\"costs\":[";
  for(size_t r=0;r<bank.costs.size();r++) out << (r?",":"") << bank.costs[r];
  out << "],\"patterns_text\":\"" << json_escape(bank.patterns) << "\"}\n";
}
void research_move(int move,double before_cost,double after_cost,double seconds,int* before,int* after) {
  static int recorded=0;
  if(recorded>=2000) return;
  std::ofstream out((export_path+".moves.jsonl").c_str(),std::ios::app);
  if(!out) throw std::runtime_error("Cannot write move profile");
  out << std::setprecision(16) << "{\"call\":" << recorded++ << ",\"move\":" << move
      << ",\"before_cost\":" << before_cost << ",\"after_cost\":" << after_cost
      << ",\"operator_seconds\":" << seconds << ",\"before_route\":";
  write_route(out,decode(before));out << ",\"after_route\":";
  write_route(out,decode(after));out << "}\n";
}
void research_export(const char* path,double value,int* route) {
  if(export_path==path) final_export=true;
  std::ofstream out(path);
  if(!out) throw std::runtime_error("Cannot write legacy route export");
  const char* method="original_MDM_legacy_instrumented";
#ifdef PCCTP_DECODER_INITIALIZED
  method="original_MDM_legacy_decoder_index_initialized";
#endif
#ifdef PCCTP_ELITE_BOUNDED
  method="original_MDM_legacy_decoder_index_initialized_elite_count_bounded";
#endif
  out << std::setprecision(16) << "{\"method\":\"" << method << "\",\"seed\":" << seed
      << ",\"clock_source\":\"clock_gettime(CLOCK_MONOTONIC)\""
      << ",\"support_absolute\":" << support << ",\"stability_percent\":" << stability
      << ",\"mining_calls\":" << calls << ",\"iterations\":" << observed_iterations
      << ",\"completed\":" << (final_export ? "true":"false") << ",\"cost\":" << value
      << ",\"seconds\":" << seconds_now()
      << ",\"route\":";
  write_route(out,decode(route)); out << ",\"snapshots\":[";
  for(size_t j=0;j<snapshots.size();j++) {
    auto& s=snapshots[j]; out << (j ? ",":"") << "{\"iteration\":" << s.iteration << ",\"cost\":" << s.cost
      << ",\"seconds\":" << s.seconds << ",\"route\":";
    write_route(out,s.route); out << "}";
  }
  out << "]";
#ifdef PCCTP_TIMED
  out << ",\"budget_seconds\":" << budget << ",\"startup_debit_seconds\":" << startup_debit
      << ",\"feasible_before_deadline\":" << (progress.empty()?"false":"true") << ",\"deadline_cost\":";
  if(progress.empty())out<<"null";else out<<progress.back().cost;
  out<<",\"deadline_route\":";if(progress.empty())out<<"[]";else write_route(out,progress.back().route);
  out<<",\"anytime\":[";
  for(size_t i=0;i<progress.size();++i) {
    auto& s=progress[i];out<<(i?",":"")<<"{\"seconds\":"<<s.seconds<<",\"cost\":"<<s.cost<<",\"route\":";
    write_route(out,s.route);out<<"}";
  }out<<"]";
#endif
  out << "}\n";
}
