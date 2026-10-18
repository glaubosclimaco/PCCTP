"""Extract actual encoder and expose reads beyond occupied elite rows."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def main():
    source=(ROOT/'mdm-grasp-genius-vnd/mineracao.cpp').read_text()
    encoder=source[source.index('void encode('):source.index('int decode(')]
    prefix='''
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
using namespace std;
struct Matriz {int linha,coluna;float* v;};
float get(Matriz* m,int i,int j) {
 if(i<0||i>=m->linha||j<0||j>=m->coluna)throw runtime_error("out-of-range elite successor");
 return m->v[i*m->coluna+j];
}
'''
    main='''
int main(int argc,char** argv) {
 float values[50];for(int i=0;i<50;++i)values[i]=1000000.;
 values[0]=1;values[1]=2;values[2]=0;values[3]=10;values[4]=3;
 Matriz m={10,5,values};
 try {encode(&m,atoi(argv[1]),3);}catch(const exception& e){return 7;}
 return 0;
}
'''
    with tempfile.TemporaryDirectory(prefix='pcctp_elite_') as folder:
        folder=Path(folder);code=folder/'test.cpp';code.write_text(prefix+encoder+main)
        binary=folder/'test';subprocess.run([shutil.which('g++-15') or 'g++','-O2',str(code),'-o',str(binary)],check=True)
        original=subprocess.run([str(binary),'10'],cwd=folder,capture_output=True)
        fixed=subprocess.run([str(binary),'1'],cwd=folder,capture_output=True)
        assert original.returncode==7 and fixed.returncode==0
        assert (folder/'bd.txt').read_text().strip()=='1 5 6'
    result=dict(status='PASS',capacity_count_reads_unoccupied_rows=True,occupied_count_encodes_exact_cycle=True,
                scope='Actual encoder with checked matrix access and poisoned unoccupied rows; not a reconstruction of the full crashed execution',
                original_source_sha256=hashlib.sha256(source.encode()).hexdigest())
    (HERE/'elite_bounds_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
