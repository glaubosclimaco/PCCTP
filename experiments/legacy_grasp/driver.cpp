#include "Graph.h"
#include "Utils.h"
#include "mineracao.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <stdexcept>
#include "patterns.hpp"

extern int qtd_de_padroes;

void construir_solucao_inicial(std::vector<Vertex>, Vertice*, Matriz*, float, int*,
    int, int, float, int, int, int, DMPadrao*, int*&);
void busca_local(Vertice*, Matriz*, int*, int*, int, int, int);

int main(int argc, char** argv) {
  if (argc != 6 && argc != 7) { std::cerr << "instance seed iterations output max_seconds [plain|mdm_reconstruction]\n"; return 2; }
  try {
    std::ifstream f(argv[1]); int n, quota, known; double radius;
    if (!(f >> n >> quota >> radius >> known) || n < 3 || n > 500)
      throw std::runtime_error("Invalid instance header");
    std::vector<int> types(n), prizes(n);
    for (int& x : types) f >> x;
    for (int& x : prizes) f >> x;
    Matriz* C = create_matriz(n,n);
    std::vector<std::vector<double>> costs(n,std::vector<double>(n));
    for (int i=0;i<n;i++) for(int j=0;j<n;j++) { f >> costs[i][j]; setMatriz(C,i,j,costs[i][j]); }
    if (!f) throw std::runtime_error("Truncated instance");
    Vertice* vertices = (Vertice*)calloc(n,sizeof(Vertice));
    int nt=0;
    for (int i=0;i<n;i++) {
      vertices[i].id=i; vertices[i].premio=prizes[i];
      strcpy(vertices[i].tipo, types[i]==0 ? "v/t" : types[i]==1 ? "t" : "w");
      nt += types[i]==1;
      std::vector<int> a,b;
      for(int j=0;j<n;j++) if(costs[i][j]<=radius) {
        if(types[i]!=2 && types[j]==2) a.push_back(j);
        if(types[i]==2 && types[j]!=2) b.push_back(j);
      }
      vertices[i].qtd_cobre=a.size(); vertices[i].qtd_coberto_por=b.size();
      vertices[i].cobre=(int*)calloc(std::max(size_t(1),a.size()),sizeof(int));
      vertices[i].coberto_por=(int*)calloc(std::max(size_t(1),b.size()),sizeof(int));
      std::copy(a.begin(),a.end(),vertices[i].cobre);
      std::copy(b.begin(),b.end(),vertices[i].coberto_por);
    }
    if(nt<3) throw std::runtime_error("Legacy GENIUS requires >=3 mandatory vertices");
    auto decode = [&](int* successor) {
      int start=0; while(start<n && successor[start]<0) ++start;
      if(start==n) throw std::runtime_error("Empty successor route");
      std::vector<int> route; std::set<int> seen;
      int v=start;
      do {
        if(v<0 || v>=n || types[v]==2 || seen.count(v)) throw std::runtime_error("Invalid successor cycle");
        route.push_back(v); seen.insert(v); v=successor[v];
      } while(v!=start);
      if(route.size()<3) throw std::runtime_error("Short cycle");
      int collected=0;
      for(int i=0;i<n;i++) {
        if(types[i]==1 && !seen.count(i)) throw std::runtime_error("Missing mandatory");
        if(seen.count(i)) collected+=prizes[i];
        if(types[i]==2) { bool covered=false; for(int u : route) covered |= costs[u][i]<=radius;
          if(!covered) throw std::runtime_error("Uncovered demand"); }
        if(successor[i]>=0 && !seen.count(i)) throw std::runtime_error("Disconnected successor cycle");
      }
      if(collected<quota) throw std::runtime_error("Insufficient prize");
      return route;
    };
    int seed=std::stoi(argv[2]), iterations=std::stoi(argv[3]);
    double max_seconds=std::stod(argv[5]);
    bool use_mining=argc==7 && std::string(argv[6])=="mdm_reconstruction";
    if(argc==7 && !use_mining && std::string(argv[6])!="plain") throw std::runtime_error("Unknown method");
    init_genrand(seed);
    std::vector<int> best; double bestcost=1e100;
    int *route=(int*)malloc(n*sizeof(int)), *neighbor=(int*)malloc(n*sizeof(int));
    int* sub=new int[4]();
    auto began=std::chrono::steady_clock::now(); int done=0;
    std::vector<std::pair<int,double>> trajectory;
    struct Snapshot { int iteration; double cost; std::vector<int> route; };
    std::vector<Snapshot> snapshots;
    std::vector<std::pair<double,Itemset>> elite;
    std::vector<OwnedPattern> patterns;
    int stable=0, mining_calls=0, mined_constructions=0, cursor=0; bool dirty=false;
    freopen("/dev/null","w",stdout);
    for(int it=0;it<iterations;it++) {
      DMPadrao* pattern=nullptr;
      if(use_mining && !patterns.empty()) {
        auto& owned=patterns[cursor++%patterns.size()];
        pattern=&owned.p; qtd_de_padroes=1; ++mined_constructions;
      }
      construir_solucao_inicial({},vertices,C,.6f,route,nt,quota,radius,n,1,pattern ? 1:0,pattern,sub);
      decode(route);
      busca_local(vertices,C,route,neighbor,radius,quota,n);
      auto solution=decode(neighbor);
      double value=0; for(size_t j=0;j<solution.size();j++) value+=costs[solution[j]][solution[(j+1)%solution.size()]];
      if(value<bestcost) { bestcost=value; best=solution; trajectory.push_back({it+1,value}); }
      if(use_mining) {
        auto tx=transaction(solution,n); bool duplicate=false;
        for(auto& member:elite) duplicate |= member.second==tx;
        bool changed=false;
        if(!duplicate && (elite.size()<10 || value<elite.back().first)) {
          elite.push_back({value,tx});
          std::sort(elite.begin(),elite.end());
          if(elite.size()>10) elite.pop_back();
          changed=true; dirty=true;
        }
        stable=changed ? 0:stable+1;
        if(stable>=5 && dirty && elite.size()>=2) {
          std::vector<Itemset> transactions; for(auto& member:elite) transactions.push_back(member.second);
          patterns=mine_patterns(transactions,n); ++mining_calls;
          cursor=0; stable=0; dirty=false;
        }
      }
      ++done;
      if(done==1 || done==5 || done==10 || done==20 || done==35 || done==70)
        snapshots.push_back({done,bestcost,best});
      double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
      if(max_seconds>0 && elapsed>=max_seconds) break;
    }
    std::ofstream out(argv[4]);
    out << std::setprecision(16) << "{\"method\":\"" << (use_mining ? "MDM_reconstruction_not_original_publication_binary":"legacy_GENIUS_RVND_not_MDM") << "\",\"seed\":" << seed
        << ",\"iterations\":" << done << ",\"cost\":" << bestcost << ",\"seconds\":"
        << std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count() << ",\"route\":[";
    for(size_t i=0;i<best.size();i++) out << (i ? "," : "") << best[i];
    out << "],\"improvements\":[";
    for(size_t i=0;i<trajectory.size();i++) out << (i ? "," : "") << "[" << trajectory[i].first << "," << trajectory[i].second << "]";
    out << "],\"mining_calls\":" << mining_calls << ",\"mined_constructions\":" << mined_constructions
        << ",\"elite_capacity\":10,\"support\":0.2,\"stability_iterations\":5,\"snapshots\":[";
    for(size_t i=0;i<snapshots.size();i++) {
      out << (i ? ",":"") << "{\"iteration\":" << snapshots[i].iteration << ",\"cost\":" << snapshots[i].cost << ",\"route\":[";
      for(size_t j=0;j<snapshots[i].route.size();j++) out << (j ? ",":"") << snapshots[i].route[j];
      out << "]}";
    }
    out << "]}\n";
    return 0;
  } catch(const std::exception& e) { std::cerr << e.what() << "\n"; return 1; }
}
