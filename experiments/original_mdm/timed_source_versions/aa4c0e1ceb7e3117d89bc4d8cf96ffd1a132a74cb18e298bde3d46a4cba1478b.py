"""Check strict exclusion of late improvements and debit of parent startup time."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

HERE=Path(__file__).resolve().parent
SOURCE='''
#include "research_hook.hpp"
#include <thread>
#include <chrono>
int main(int argc,char** argv) {
 research_process_start();research_initialize("toy",4,42,2,5,argv[1]);
 int early[4]={1,2,0,-1},late[4]={1,3,-1,0};
 research_progress(10.,early);
 std::this_thread::sleep_for(std::chrono::milliseconds(40));
 research_progress(5.,late);research_export(argv[1],5.,late);
}
'''

def main():
    with tempfile.TemporaryDirectory(prefix='pcctp_deadline_') as folder:
        folder=Path(folder);(folder/'main.cpp').write_text(SOURCE);binary=folder/'test'
        subprocess.run([shutil.which('g++-15') or 'g++','-std=c++11','-O2','-DPCCTP_TIMED',
                        '-I',str(HERE),str(folder/'main.cpp'),str(HERE/'research_hook.cpp'),'-o',str(binary)],check=True)
        env=os.environ.copy();env['PCCTP_BUDGET']='0.02';env.pop('PCCTP_STARTED_AT',None)
        out=folder/'normal.json';subprocess.run([str(binary),str(out)],env=env,check=True)
        r=json.loads(out.read_text());assert r['cost']==5 and r['deadline_cost']==10
        assert len(r['anytime'])==1 and r['deadline_route']==[0,1,2]
        assert r['anytime'][0]['seconds']<=.02
        env['PCCTP_STARTED_AT']=format(time.time()-1.,'.17g');out=folder/'startup.json'
        subprocess.run([str(binary),str(out)],env=env,check=True)
        r=json.loads(out.read_text());assert not r['feasible_before_deadline'] and r['deadline_cost'] is None
        assert r['startup_debit_seconds']>=1.
    result=dict(status='PASS',late_improvement_excluded=True,parent_startup_debited=True,
                fallback_after_deadline_marked_unavailable=True,
                source_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                               ['verify_deadline.py','research_hook.cpp','research_hook.hpp']})
    (HERE/'deadline_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
