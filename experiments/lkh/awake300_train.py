"""Repeat all81 frozen TRAIN runs after documented host-sleep interference."""
import datetime,fcntl,hashlib,json,os,random,subprocess,threading,time
from pathlib import Path
import ruin_total300_train as original
HERE=Path(__file__).resolve().parent;OUT=HERE/'ruin_total300_awake_train'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

class ClockGuard:
    def __init__(self):
        self.stop_event=threading.Event();self.max_gap=0.;self.clock_difference=0.
        self.started=time.time();self.thread=threading.Thread(target=self.watch,daemon=True)
    def watch(self):
        last=time.monotonic();wall=time.time()
        while not self.stop_event.wait(1.):
            now=time.monotonic();current=time.time()
            self.max_gap=max(self.max_gap,now-last,current-wall)
            self.clock_difference=max(self.clock_difference,abs((current-wall)-(now-last)))
            last,wall=now,current
    def __enter__(self):self.thread.start();return self
    def __exit__(self,*args):self.stop_event.set();self.thread.join()

def sleep_events_since(start):
    text=subprocess.run(['pmset','-g','log'],capture_output=True,text=True,check=True).stdout
    events=[]
    for line in text.splitlines():
        if 'Entering Sleep' not in line:continue
        try:when=datetime.datetime.strptime(line[:25],'%Y-%m-%d %H:%M:%S %z').timestamp()
        except ValueError:continue
        if when>=start-1:events.append(line)
    return events

if __name__=='__main__':
    lock=(HERE/'awake300_train.queue.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    lock.seek(0);lock.truncate();lock.write(str(os.getpid()));lock.flush()
    assert (HERE/'ruin_total300_train/environment_interruption.json').exists()
    OUT.mkdir(exist_ok=True);original.OUT=OUT
    base=json.loads((HERE/'ruin_total300_train/protocol.json').read_text())
    for name,hashed in base['source_hashes'].items():assert digest(HERE/name)==hashed
    for name,hashed in base['dependency_hashes'].items():assert digest(original.ROOT/name)==hashed
    for name,hashed in base['binary_hashes'].items():assert digest(Path(name))==hashed
    protocol=dict(base,repeat_reason='entire81-run repetition after documented host sleep; no outcome selection',
        original_interruption_sha256=digest(HERE/'ruin_total300_train/environment_interruption.json'),
        awake_runner_sha256=digest(Path(__file__)),max_monitor_gap_seconds=5.,
        repetition_protocol_sha256=digest(original.ROOT/'research_notes/awake300_repetition_protocol.md'),
        max_wall_monotonic_difference_seconds=3.,require_no_pmset_sleep_event=True)
    target=OUT/'protocol.json'
    if target.exists():assert json.loads(target.read_text())==protocol
    else:target.write_text(json.dumps(protocol,indent=2)+'\n')
    for c in base['cases']:
        for seed in (70,71,72):
            block=OUT/f"environment.{c['instance']}.seed{seed}.json"
            if block.exists():
                prior=json.loads(block.read_text());assert prior['status']=='PASS'
                for mode,hashed in prior['result_hashes'].items():
                    assert digest(OUT/f"{c['instance']}.seed{seed}.{mode}.json")==hashed
                continue
            # A partial, unguarded block must never silently become a valid resume.
            assert not any((OUT/f"{c['instance']}.seed{seed}.{m}.json").exists() for m in base['arms'])
            order=list(base['arms']);random.Random(c['instance']+str(seed)).shuffle(order)
            with ClockGuard() as guard:
                for mode in order:
                    r=original.run(c['instance'],seed,mode)
                    print(json.dumps({k:r[k] for k in ('instance','seed','mode','cost','wall_seconds')}),flush=True)
            events=sleep_events_since(guard.started)
            status='PASS' if guard.max_gap<=5 and guard.clock_difference<=3 and not events else 'ENVIRONMENT_FAILURE'
            record=dict(status=status,instance=c['instance'],seed=seed,order=order,
                started_unix=guard.started,finished_unix=time.time(),max_monitor_gap_seconds=guard.max_gap,
                max_clock_difference_seconds=guard.clock_difference,sleep_events=events,
                result_hashes={m:digest(OUT/f"{c['instance']}.seed{seed}.{m}.json") for m in base['arms']})
            block.write_text(json.dumps(record,indent=2)+'\n')
            assert status=='PASS','Host suspension/scheduling interruption: preserved block, no decision or next run'
    blocks=list(OUT.glob('environment.g_*.json'));assert len(blocks)==27
    (OUT/'environment_audit.json').write_text(json.dumps(dict(status='PASS',blocks=27,
        records={p.name:digest(p) for p in blocks},scope='monitor+pmset screening; not a guarantee of equal CPU/thermal conditions'),indent=2)+'\n')
    (OUT/'completion.json').write_text(json.dumps(dict(runs=81))+'\n')
