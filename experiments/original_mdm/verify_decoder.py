"""Minimal decoder regression: a one-edge pattern needs component index zero.

The unmodified source declares cc_atual without initialization, then uses it
when a pattern has no adjacent pair. Test the isolated initialized variant;
do not silently replace the recorded original baseline.
"""
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
DRIVER='''#include "mineracao.h"
int decode(int,DMPadrao*&);
float get(Matriz*,int,int) { return 0; }
int main() {
  DMPadrao* patterns=0;
  int count=decode(7,patterns);
  if(count!=1 || patterns[0].qtd_cc!=1) return 1;
  auto &component=patterns[0].ccs[0];
  if(component.tamanho!=1 || component.aresta_id[0]!=0) return 2;
  if(patterns[0].arestas[0].vertice_origem!=0 || patterns[0].arestas[0].vertice_destino!=1) return 3;
  return 0;
}
'''

def main():
    build=HERE/'build_decoder_initialized'
    with tempfile.TemporaryDirectory(prefix='pcctp-single-edge-') as folder:
        folder=Path(folder); driver=folder/'driver.cpp'; driver.write_text(DRIVER)
        (folder/'padroes.txt').write_text('1;2;1\n')
        binary=folder/'verify'
        command=[shutil.which('g++-15') or 'g++','-std=gnu++11','-O2','-I',str(build),
                 str(driver),str(build/'mineracao.cpp'),'-o',str(binary)]
        subprocess.run(command,check=True,capture_output=True,text=True)
        process=subprocess.run([str(binary)],cwd=folder,capture_output=True,text=True,timeout=10)
        assert process.returncode==0,process
    result=dict(status='PASS',fixture='one directed edge 0->1, support 2',
                expected_components=1,expected_component_index=0,
                method='separate legacy variant with initialized decoder index')
    (HERE/'decoder_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__': main()
