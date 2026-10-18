// Exact Held--Karp repair on fixed path components, with optional supernodes.
// Finite nonnegative transition costs; optional-selection feasibility is supplied
// separately. This is an experimental backend, not a novelty claim for DP.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

extern "C" int pcctp_path_dp(int n,int root,const double* edge,const int* orientations,
  const unsigned char* feasible,double upper,double seconds,int* order,int* length,
  double* value,std::uint64_t* processed) {
  if(n<1 || n>20 || root<0 || root>=n) return 2;
  const int width=2*n; const std::uint32_t total=1u<<n;
  const double inf=std::numeric_limits<double>::infinity();
  auto began=std::chrono::steady_clock::now();
  std::vector<double> dp(std::size_t(total)*width,inf);
  std::vector<unsigned char> parent(std::size_t(total)*width,255);
  const std::uint32_t start=1u<<root;
  dp[std::size_t(start)*width+2*root]=0.;
  double best=upper; std::uint32_t bestmask=0; int bestlast=-1; int status=0;
  *processed=0; *length=0;
  std::uint32_t scanned=0;
  for(std::uint32_t mask=start;mask<total;++mask) {
    if(!(mask&start)) continue;
    if((++scanned&255u)==0 && std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()>=seconds) {
      status=1; break;
    }
    for(int last=0;last<n;++last) if(mask&(1u<<last)) {
      for(int orientation=0;orientation<orientations[last];++orientation) {
        const int state=2*last+orientation;
        double base=dp[std::size_t(mask)*width+state];
        if(!std::isfinite(base) || base>=best) continue;
        ++*processed;
        if(feasible[mask]) {
          const double closed=base+edge[state*width+2*root];
          if(closed<best) { best=closed;bestmask=mask;bestlast=state; }
        }
        for(int next=0;next<n;++next) if(!(mask&(1u<<next))) {
          auto target=mask|(1u<<next);
          for(int o=0;o<orientations[next];++o) {
            int nextstate=2*next+o;
            double candidate=base+edge[state*width+nextstate];
            auto index=std::size_t(target)*width+nextstate;
            if(candidate<best && candidate<dp[index]) { dp[index]=candidate; parent[index]=state; }
          }
        }
      }
    }
  }
  *value=best;
  if(bestlast>=0) {
    std::vector<int> reversed; auto mask=bestmask; int state=bestlast;
    while(true) {
      reversed.push_back(state);
      if(mask==start && state==2*root) break;
      int previous=parent[std::size_t(mask)*width+state];
      if(previous==255 || reversed.size()>std::size_t(n)) return 3;
      mask^=1u<<(state/2); state=previous;
    }
    std::reverse(reversed.begin(),reversed.end());
    *length=int(reversed.size());std::copy(reversed.begin(),reversed.end(),order);
  }
  return status;
}
