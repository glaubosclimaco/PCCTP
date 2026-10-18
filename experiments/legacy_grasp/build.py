"""Build a read-only snapshot of legacy GENIUS/RVND without Gurobi/LEMON.

The checkout contains Utils.h and utils.h, which collide on case-insensitive
filesystems. Restore the capitalized heuristic header from Git in an isolated
directory, never modifying the user's working copy. This is NOT an MDM baseline.
"""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]
DEST = Path(__file__).with_name('build')
DEST.mkdir(exist_ok=True)
sources = ['grasp.cpp', 'genius.cpp', 'Utils.cpp', 'random_provider.cpp', 'PR.cpp', 'mineracao.cpp']
headers = ['Graph.h', 'PR.h', 'mineracao.h', 'random_provider.h']
manifest = {'purpose': 'legacy GENIUS/RVND only; mining inactive', 'sha256': {}}
for name in sources:
    data = (ROOT/'src'/name).read_bytes()
    if not (DEST/name).exists() or (DEST/name).read_bytes()!=data:
        (DEST/name).write_bytes(data)
    manifest['sha256']['src/'+name] = hashlib.sha256(data).hexdigest()
for name in headers:
    data = (ROOT/'include'/name).read_bytes()
    if not (DEST/name).exists() or (DEST/name).read_bytes()!=data:
        (DEST/name).write_bytes(data)
    manifest['sha256']['include/'+name] = hashlib.sha256(data).hexdigest()
data = subprocess.check_output(['git', 'show', 'HEAD:include/Utils.h'], cwd=ROOT)
manifest['sha256']['HEAD:include/Utils.h'] = hashlib.sha256(data).hexdigest()
text = data.decode().replace('const char tipo[4]', 'char tipo[4]')
if not (DEST/'Utils.h').exists() or (DEST/'Utils.h').read_text()!=text:
    (DEST/'Utils.h').write_text(text)
manifest['header_adjustment'] = 'Remove const from tipo storage; layout unchanged; loader writes type labels.'
compiler = shutil.which('g++-15') or shutil.which('g++')
flags = [compiler, '-std=gnu++11', '-O2', '-fpermissive', '-Wno-narrowing',
           '-I', str(DEST), '-include', str(DEST/'Graph.h'),
           '-I', str(Path(__file__).parent)]
objects=[]
for name in sources:
    obj=DEST/(name+'.o'); src=DEST/name
    if not obj.exists() or obj.stat().st_mtime < max(p.stat().st_mtime for p in [src,*DEST.glob('*.h')]):
        subprocess.run(flags+['-c',str(src),'-o',str(obj)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    objects.append(str(obj))
command = flags + objects + [str(Path(__file__).with_name('driver.cpp')), '-o', str(DEST/'legacy_grasp')]
process = subprocess.run(command, capture_output=True, text=True)
(DEST/'compiler.log').write_text(process.stdout+process.stderr)
manifest['command'] = command
manifest['compiler'] = subprocess.check_output([compiler, '--version'], text=True).splitlines()[0]
(DEST/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
if process.returncode:
    print(process.stderr[-12000:])
    raise SystemExit(process.returncode)
print(DEST/'legacy_grasp')
