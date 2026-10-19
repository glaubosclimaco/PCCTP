// Symmetric integer TSP first-improvement 2-opt; standalone research control.
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
int main(){
 int n;double seconds;if(!(std::cin>>n>>seconds)||n<3)return 2;
 std::vector<long long> c(n*n);for(auto &x:c)if(!(std::cin>>x))return 2;
 std::vector<int> r(n);for(auto &v:r)if(!(std::cin>>v)||v<0||v>=n)return 2;
 auto sorted=r;std::sort(sorted.begin(),sorted.end());for(int i=0;i<n;i++)if(sorted[i]!=i)return 2;
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(c[i*n+j]!=c[j*n+i])return 2;
 auto end=std::chrono::steady_clock::now()+std::chrono::duration<double>(seconds);
 bool changed=true,complete=false;
 while(changed && std::chrono::steady_clock::now()<end){
  changed=false;
  for(int i=0;i<n-1 && !changed;i++)for(int j=i+2;j<n;j++){
   if(std::chrono::steady_clock::now()>=end)goto done;
   int a=r[i],b=r[i+1],d=r[j],e=r[(j+1)%n];
   auto delta=c[a*n+d]+c[b*n+e]-c[a*n+b]-c[d*n+e];
   if(delta<0){std::reverse(r.begin()+i+1,r.begin()+j+1);changed=true;break;}
  }
  if(!changed)complete=true;
 }
 done:
 long long cost=0;for(int i=0;i<n;i++)cost+=c[r[i]*n+r[(i+1)%n]];
 std::cout<<cost<<' '<<complete<<'\n';for(int v:r)std::cout<<v<<' ';std::cout<<'\n';
}
