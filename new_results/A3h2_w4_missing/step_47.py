import sys, pickle, datetime
sys.path.insert(0,'<project>/research/08_theory/A3h2_w4_missing')
BASE='<project>'
WORK=BASE+'/research/08_theory/A3h2_w4_missing'
with open(WORK+'/progress.log','a') as f:
    f.write(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} step3 DONE: 29 missing relations FOUND (single-embedding LLL C=1e300, vals110 N=2000; gap 1e46 to accidents; 2x-verified 1e-679 at N=3800/prec2500; rank 29 mod p; heights 1,122..264; support 6,105..109; only N0 form-level (6-term), other 28 exotic fail k=3,5). Saved w4_missing_relations.json, cands29_reduced.pickle, vals110_k1_700.pickle\n")
print('logged')
# load rref pickles
for p in (1000000009,1000000033):
    d=pickle.load(open(BASE+f'/research/08_theory/A3h_weight4_target/rref_T3_{p}.pkl','rb'))
    print(p, d.keys(), 'rank',d['rank'],'quotient',d['quotient'])
    print('  free[:5]',d['free'][:5],'len',len(d['free']),'nf keys',d['nf'].keys())
    print('  r',d['r'])
