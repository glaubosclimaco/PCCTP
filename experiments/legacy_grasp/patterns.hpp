// Exact maximal frequent undirected edge sets by support-subset intersections.
// At most ten elite transactions; no external mining executable is required.
#pragma once
#include "mineracao.h"
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

using Itemset = std::set<int>;
inline Itemset transaction(const std::vector<int>& route, int n) {
  Itemset result;
  for(size_t j=0;j<route.size();j++) {
    int a=route[j], b=route[(j+1)%route.size()];
    if(a>b) std::swap(a,b);
    result.insert(a*n+b);
  }
  return result;
}

struct OwnedPattern {
  DMPadrao p{};
  std::vector<Aresta> arcs;
  std::vector<ComponenteConexa> components;
  std::vector<std::vector<int>> ids;
  void bind() {
    p.id=0; p.tamanho=arcs.size(); p.arestas=arcs.data();
    p.qtd_cc=components.size(); p.cc_atual=p.qtd_cc-1;
    for(size_t i=0;i<ids.size();i++) {
      components[i].tamanho=ids[i].size(); components[i].aresta_id=ids[i].data();
    }
    p.ccs=components.data();
  }
};

inline std::vector<OwnedPattern> mine_patterns(const std::vector<Itemset>& transactions, int n) {
  const int count=transactions.size();
  if(count<2 || count>10) return {};
  const int minimum=std::max(2,(int)std::ceil(.2*count));
  std::set<Itemset> intersections;
  for(unsigned mask=1;mask<(1u<<count);mask++) {
    if(__builtin_popcount(mask)!=minimum) continue;
    Itemset common; bool first=true;
    for(int j=0;j<count;j++) if(mask&(1u<<j)) {
      if(first) { common=transactions[j]; first=false; }
      else {
        Itemset next;
        std::set_intersection(common.begin(),common.end(),transactions[j].begin(),transactions[j].end(),
                              std::inserter(next,next.begin()));
        common.swap(next);
      }
    }
    if(!common.empty()) intersections.insert(common);
  }
  std::vector<Itemset> maximal;
  for(const auto& a : intersections) {
    bool dominated=false;
    for(const auto& b : intersections) if(b.size()>a.size() && std::includes(b.begin(),b.end(),a.begin(),a.end())) {
      dominated=true; break;
    }
    if(!dominated) maximal.push_back(a);
  }
  std::sort(maximal.begin(),maximal.end(),[](const Itemset& a,const Itemset& b) {
    return a.size()!=b.size() ? a.size()>b.size() : a<b;
  });
  std::vector<OwnedPattern> result;
  for(const auto& pattern : maximal) {
    if(result.size()>=10) break;
    std::map<int,std::vector<int>> adj;
    for(int edge : pattern) { int a=edge/n,b=edge%n; adj[a].push_back(b); adj[b].push_back(a); }
    bool degree_ok=true; for(auto& entry:adj) degree_ok &= entry.second.size()<=2;
    if(!degree_ok) continue;
    std::set<int> unused; for(const auto& entry:adj) unused.insert(entry.first);
    std::vector<std::vector<int>> paths;
    while(!unused.empty()) {
      int start=*unused.begin();
      // Find the connected component first; orient its path from an endpoint.
      std::set<int> component; std::vector<int> queue{start};
      for(size_t q=0;q<queue.size();q++) if(component.insert(queue[q]).second)
        for(int v:adj[queue[q]]) if(!component.count(v)) queue.push_back(v);
      bool endpoint=false;
      for(int v:component) if(adj[v].size()==1) { start=v; endpoint=true; break; }
      for(int v:component) unused.erase(v);
      if(!endpoint) continue; // closed cycles cannot safely seed the legacy path constructor
      std::vector<int> path{start}; int prev=-1, current=start;
      while(true) {
        int next=-1; for(int v:adj[current]) if(v!=prev) { next=v; break; }
        if(next<0) break;
        path.push_back(next); prev=current; current=next;
      }
      paths.push_back(path);
    }
    std::sort(paths.begin(),paths.end(),[](const std::vector<int>& a,const std::vector<int>& b) {
      return a.size()!=b.size() ? a.size()<b.size() : a<b;
    });
    if(paths.empty() || paths.back().size()<3) continue;
    OwnedPattern owned;
    int support=0; for(const auto& tx:transactions) support+=std::includes(tx.begin(),tx.end(),pattern.begin(),pattern.end());
    owned.p.suporte=support;
    for(const auto& path:paths) {
      owned.ids.emplace_back();
      for(size_t j=1;j<path.size();j++) {
        int id=owned.arcs.size();
        Aresta a{}; a.id=id; a.vertice_origem=path[j-1]; a.vertice_destino=path[j];
        owned.arcs.push_back(a); owned.ids.back().push_back(id);
      }
    }
    owned.components.resize(owned.ids.size()); result.push_back(std::move(owned));
  }
  for(auto& pattern:result) pattern.bind();
  return result;
}
