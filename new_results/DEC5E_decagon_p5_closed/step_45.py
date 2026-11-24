import subprocess, time, os, sys, pickle
PY = '<project>/.venv/bin/python'
# launch big table builds in background
procs = {}
for dps, N in [(1150, 3900), (1550, 5250)]:
    log = open(f'build_{dps}.log','w')
    p = subprocess.Popen(['nice', PY, 'step_28_build_tables.py', str(dps), str(N), '10'],
                         stdout=log, stderr=subprocess.STDOUT, cwd='.')
    procs[dps] = p
print("launched:", {k: v.pid for k, v in procs.items()})
# meanwhile inspect DEC5D's exact block symbols
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib as D
bs = pickle.load(open('../DEC5D_decagon_order5_depth/dec_block_symbols.pkl','rb'))
print(type(bs), list(bs.keys()) if isinstance(bs, dict) else len(bs))
for k in bs:
    v = bs[k]
    print(k, type(v), (len(v) if hasattr(v,'__len__') else v))
