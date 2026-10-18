"""Compile the supplied original MDM + FPMax in an isolated instrumented copy."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'mdm-grasp-genius-vnd'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--initialize-decoder-index',action='store_true')
parser.add_argument('--profile-moves',action='store_true')
parser.add_argument('--scheduled',action='store_true')
args=parser.parse_args()
if args.profile_moves and not args.initialize_decoder_index:
    parser.error('Move profiling requires the explicitly corrected decoder baseline')
if args.scheduled and (not args.initialize_decoder_index or args.profile_moves):parser.error('Scheduled build requires corrected decoder and no profiling')
BUILD=Path(__file__).with_name('build_scheduled' if args.scheduled else 'build_move_profile' if args.profile_moves else 'build_decoder_initialized' if args.initialize_decoder_index else 'build')
BUILD.mkdir(exist_ok=True)
manifest={'origin':str(SOURCE),'sha256':{},'changes':[]}
for name in ['Utils.cpp','genius.cpp','grasp.cpp','instances.cpp','instancia_loader.cpp',
             'mineracao.cpp','principal.cpp','random_provider.cpp','random_provider.h','mineracao.h']:
    data=(SOURCE/name).read_bytes()
    manifest['sha256'][name]=hashlib.sha256(data).hexdigest()
    text=data.decode().replace('\r\n','\n')
    if name=='mineracao.cpp' and args.initialize_decoder_index:
        anchor='int cc_atual;'
        assert text.count(anchor)==1
        text=text.replace(anchor,'int cc_atual = 0; // initialized for a singleton pattern component')
        manifest['changes'].append('EXPLICIT BUG FIX: initialize decoder component index before singleton-component access; separate from unchanged legacy')
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
        if args.profile_moves or args.scheduled:
            declaration='void obter_vizinho(Vertice *vertices, Matriz* C, int mov, int *&rota,\n\t\tint *vertices_solucao, int tamanho_rota, int d, int prize, int qtd_n)'
            assert text.count(declaration)==1
            text=text.replace(declaration,declaration.replace('obter_vizinho','profiled_obter_vizinho'))
            wrapper='''
void obter_vizinho(Vertice *vertices, Matriz* C, int mov, int *&rota,
        int *vertices_solucao, int tamanho_rota, int d, int prize, int qtd_n) {
    std::vector<int> before(rota,rota+qtd_n);
    double before_cost=custo(C,rota);
    auto start=std::chrono::steady_clock::now();
    profiled_obter_vizinho(vertices,C,mov,rota,vertices_solucao,tamanho_rota,d,prize,qtd_n);
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    research_move(mov,before_cost,custo(C,rota),elapsed,before.data(),rota);
}
'''
            anchor='void busca_local(Vertice *vertices'
            assert text.count(anchor)==1
            text=text.replace(anchor,wrapper+'\n'+anchor)
            text='#include <vector>\n#include <chrono>\n'+text
            manifest['changes'].append('OBSERVATIONAL move profiling: copy successor array, time each operator, export first 2000 calls; no RNG/decision change')
            if args.scheduled:
                text=text.replace('    std::vector<int> before(rota,rota+qtd_n);\n','')
                text=text.replace('research_move(mov,before_cost,custo(C,rota),elapsed,before.data(),rota);',
                                  'schedule_update(mov,before_cost-custo(C,rota),elapsed);')
                # Include after Vertice is defined; retain the original RNG permutation.
                marker='std::string nome_movimento(int movimentoEnum)'
                text=text.replace(marker,'#include "scheduler.hpp"\n\n'+marker)
                marker='\t//BUSCA LOCAL: VND\n\twhile (ordem < qtd_vizinhanca)'
                assert text.count(marker)==1
                text=text.replace(marker,'\tschedule_order(vertices,C,s_melhor_vizinho,prize,qtd_n,vizinhancas);\n'+marker)
                before,after=text.split('void busca_local_v2(',1)
                marker='\t\t\tordem = 0;'
                assert before.count(marker)==1
                before=before.replace(marker,marker+'\n\t\t\tschedule_order(vertices,C,s_melhor_vizinho,prize,qtd_n,vizinhancas);')
                text=before+'void busca_local_v2('+after
                manifest['changes'][-1]='EXPERIMENTAL scheduling: common per-operator timing; reorder all 15 operators after state improvements; original RNG permutation retained'
    target=BUILD/name
    if not target.exists() or target.read_text()!=text: target.write_text(text)
for name in ['research_hook.hpp','research_hook.cpp']:
    shutil.copyfile(Path(__file__).with_name(name),BUILD/name)
if args.scheduled:
    home=Path(__file__).parent
    frozen=home/'operator_models/models.52b8e9ae6089.json'
    models=json.loads(frozen.read_text())['models']
    names=['operator_only','operator_size','operator_generic','operator_state']
    header='// Frozen TRAIN-only models; ABI modes 0/1/2/3.\n'
    for key in ['mean','scale','probability','timing']:
        for i,name in enumerate(names):
            values=','.join(format(v,'.17g') for v in models[name][key])
            header+=f'static const double {key}_{i}[]={{'+values+'};\n'
        header+=f'static const double* model_{key}[]={{'+','.join(f'{key}_{i}' for i in range(4))+'};\n'
    (BUILD/'schedule_models.hpp').write_text(header)
    manifest['model_sha256']=hashlib.sha256(frozen.read_bytes()).hexdigest()
    for name in ['scheduler.hpp','operator_state.cpp','operator_inference.cpp']:
        shutil.copyfile(home/name,BUILD/name)
        manifest['sha256'][name]=hashlib.sha256((home/name).read_bytes()).hexdigest()
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
if args.initialize_decoder_index: flags.append('-DPCCTP_DECODER_INITIALIZED')
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
if args.scheduled:sources+=['operator_state.cpp','operator_inference.cpp']
command=flags+['-I',str(BUILD),*[str(BUILD/name) for name in sources],'-o',str(BUILD/'original_mdm')]
run(command)
manifest['command']=command
manifest['compiler']=subprocess.check_output([compiler,'--version'],text=True).splitlines()[0]
(BUILD/'compiler.log').write_text('\n'.join(log))
(BUILD/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(BUILD/'original_mdm')
