import pickle, os, subprocess
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
for s in ['pp','pm','mp','mm']:
    d = pickle.load(open(f'work/step01_{s}.pkl','rb'))
    print(s, 'rank', d['rank'], '|S|', len(d['S']), '|C|', len(d['C']), 'fam', d['fam_counts'], 'fed', d.get('fed'), 'skipped', d.get('skipped'), 'ncols', d['ncols'])
print(subprocess.run(['bash','-c','grep -E "step_01b mm" step_01b.out | tail -6'], capture_output=True, text=True).stdout)
with open('progress.log','a') as f:
    f.write('2026-10-29 04:18:06 step_01b ALL DONE: pp |S|=12156 (eng12120+fw5 36), pm 11101 (eng11065+fw5 36), mp 11041 (eng only), mm 11999 (eng only). D_rat & D_sq in span mod p1 in all 4 sectors.\n')