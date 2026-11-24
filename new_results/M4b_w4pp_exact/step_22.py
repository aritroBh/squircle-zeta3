import os, json, pickle, sys, time
ROOT = '<project>'
WORK = os.path.join(ROOT, 'research/08_theory/M4b_w4pp_exact')
os.makedirs(WORK, exist_ok=True)
os.chdir(WORK)
A3H5 = os.path.join(ROOT,'research/08_theory/A3h5_exact_certificate')

# check stage1_exact_z.json 'chosen'
z = json.load(open(os.path.join(A3H5,'stage1_exact_z.json')))
print('chosen:', z['chosen'])
print('n labels:', len(z['labels']))
print(z['labels'][:5])
# check rref pickles
for p in [1000000009, 1000000033, 1000000097, 1000000241, 1000000289, 1000000297, 1000000321, 1000000409]:
    for d in [A3H5, os.path.join(ROOT,'research/08_theory/A3h_weight4_target')]:
        f = os.path.join(d, f'rref_T0_{p}.pkl')
        if os.path.exists(f):
            R = pickle.load(open(f,'rb'))
            print(p, os.path.basename(d), 'rank', R['rank'], 'free', len(R['free']))
            break
# check stage2_rows_S
S = pickle.load(open(os.path.join(A3H5,'stage2_rows_S.pkl'),'rb'))
print('stage2 rows S:', len(S['S']), 'p', S['p'], 'rank', S['rank'])
