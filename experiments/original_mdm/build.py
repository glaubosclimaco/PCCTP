"""Compile the supplied original MDM + FPMax in an isolated instrumented copy."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'mdm-grasp-genius-vnd'
BUILD=Path(__file__).with_name('build')
BUILD.mkdir(exist_ok=True)
manifest={'origin':str(SOURCE),'sha256':{},'changes':[]}
for name in ['Utils.cpp','genius.cpp','grasp.cpp','instances.cpp','instancia_loader.cpp',
             'mineracao.cpp','principal.cpp','random_provider.cpp','random_provider.h','mineracao.h']:
    data=(SOURCE/name).read_bytes()
    manifest['sha256'][name]=hashlib.sha256(data).hexdigest()
    text=data.decode().replace('\r\n','\n')
    if name=='principal.cpp':
        text='#include "research_hook.hpp"\n'+text
        text=text.replace('if (argc != 5)', 'if (argc != 6)').replace('if (argc == 5)', 'if (argc == 6)')
        text=text.replace('char file_name[100]', 'char file_name[1024]')
        anchor='clock_t start = clock();'
        assert text.count(anchor)==1
        text=text.replace(anchor,'research_initialize(file_name, qtd_n, seed, sup_min, percentual_estabilidade, argv[5]);\n\t'+anchor)
        anchor='clock_t end = clock();'
        assert text.count(anchor)==1
        text=text.replace(anchor,anchor+'\n\tresearch_export(argv[5], custo_obtido, s_melhor_rota);')
        manifest['changes'].append('principal: sixth output argument, larger path buffer, observational initialization/export')
    if name=='grasp.cpp':
        text='#include "research_hook.hpp"\n'+text
        anchor='printf("Custo solucao busca local %f \\n", custo_vizinho );'
        assert text.count(anchor)==1, 'Snapshot anchor changed'
        text=text.replace(anchor,'research_snapshot(iteracao+1, melhor_custo, s_melhor_rota);\n\t\t'+anchor)
        anchor='qtd_de_padroes = mining(conjunto_elite, elite_size, qtd_n, padroes, sup_min, max_padroes);'
        assert text.count(anchor)==1
        text=text.replace(anchor,'research_mining();\n\t\t\t'+anchor)
        manifest['changes'].append('grasp: observational checkpoint and mining-call hooks; algorithm unchanged')
    target=BUILD/name
    if not target.exists() or target.read_text()!=text: target.write_text(text)
for name in ['research_hook.hpp','research_hook.cpp']:
    shutil.copyfile(Path(__file__).with_name(name),BUILD/name)
fp=BUILD/'FPMax'; fp.mkdir(exist_ok=True)
for source in (SOURCE/'Debug/FPMax').iterdir():
    if source.suffix in ('.cpp','.c','.h'):
        data=source.read_bytes()
        text=data.replace(b'#include <malloc.h>', b'#include <cstdlib>')
        (fp/source.name).write_bytes(text)
        manifest['sha256']['FPMax/'+source.name]=hashlib.sha256(source.read_bytes()).hexdigest()
        if text!=data: manifest['changes'].append('FPMax/'+source.name+': malloc.h -> cstdlib for macOS')
compiler=shutil.which('g++-15') or shutil.which('g++')
flags=[compiler,'-std=gnu++11','-O2','-fpermissive','-Wno-narrowing']
log=[]
def run(command):
    process=subprocess.run(command,capture_output=True,text=True)
    log.append(process.stdout+process.stderr)
    if process.returncode:
        (BUILD/'compiler.log').write_text('\n'.join(log))
        print(process.stderr[-10000:]); raise SystemExit(process.returncode)
run(flags+['-DMFI',*[str(fp/name) for name in ['fpmax.cpp','data.cpp','fp_tree.cpp','buffer.cpp','fsout.cpp','fp_node.cpp']],
           '-o',str(fp/'fim_maximal')])
run(flags+[str(fp/'fpmax_hnmp.cpp'),'-o',str(fp/'fpmax_hnmp')])
sources=['Utils.cpp','genius.cpp','grasp.cpp','instances.cpp','instancia_loader.cpp','mineracao.cpp',
         'principal.cpp','random_provider.cpp','research_hook.cpp']
command=flags+['-I',str(BUILD),*[str(BUILD/name) for name in sources],'-o',str(BUILD/'original_mdm')]
run(command)
manifest['command']=command
manifest['compiler']=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]
(BUILD/'compiler.log').write_text('\n'.join(log))
(BUILD/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(BUILD/'original_mdm')
