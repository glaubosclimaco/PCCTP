"""Check supplied FPMax against exhaustive maximal frequent sets on tiny data."""
from itertools import combinations
from pathlib import Path
import json
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent

def parse(path):
    result={}
    for line in path.read_text().splitlines():
        size,support,items=line.split(';',2)
        pattern=frozenset(map(int,items.split()))
        assert len(pattern)==int(size)
        if pattern: result[pattern]=int(support)
    return result

def main():
    transactions=[{1,2,3},{1,2,4},{1,2},{2,3,4},{1,3,4},{1,2,3,4}]
    with tempfile.TemporaryDirectory(prefix='pcctp-fpmax-') as folder:
        folder=Path(folder); database=folder/'transactions.txt'
        database.write_text(''.join(' '.join(map(str,sorted(t)))+'\n' for t in transactions))
        allsets=[frozenset(s) for k in range(1,5) for s in combinations(range(1,5),k)]
        checked=0
        for threshold in range(1,len(transactions)+1):
            frequent={s:sum(s<=t for t in transactions) for s in allsets
                      if sum(s<=t for t in transactions)>=threshold}
            expected={s:n for s,n in frequent.items() if not any(s<other for other in frequent)}
            output=folder/'patterns.txt'
            subprocess.run([str(HERE/'build/FPMax/fim_maximal'),str(database),str(threshold),str(output)],
                           check=True,capture_output=True,text=True)
            actual=parse(output)
            assert actual==expected,(threshold,actual,expected)
            checked+=1
    result=dict(status='PASS',thresholds=checked,transactions=len(transactions),
                scope='Fixed-support FPMax correctness on exhaustive tiny fixture; wrapper ranking not certified')
    (HERE/'fpmax_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__': main()
