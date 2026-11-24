import sys; sys.path.insert(0,'.')
import importlib, symbol_engine_r as SE
importlib.reload(SE)
import random, itertools

def shuffle_words(u, v):
    out = []
    def rec(i, j, pre):
        if i==len(u) and j==len(v): out.append(pre); return
        if i<len(u): rec(i+1, j, pre+(u[i],))
        if j<len(v): rec(i, j+1, pre+(v[j],))
    rec(0,0,())
    return out

def tshuffle1(S1, S2p):
    out = {}
    for (i,),c1 in S1.items():
        for (j,k),c2 in S2p.items():
            for key in [(i,j,k),(j,i,k),(j,k,i)]:
                out[key] = out.get(key,0)+c1*c2
                if out[key]==0: del out[key]
    return out

# shuffle homomorphism on FULLY convergent distinct-letter words only
random.seed(2)
ok = True; cnt = 0
for s in [1,2,3]:
    L, idx = SE.letters(s)
    A = [('0',),('1',)] + [('p',i) for i in range(s)]
    conv1 = [a for a in A if a not in (('0',),('1',))]          # p_i only: G(a;1) convergent & != 1
    for trial in range(60):
        a = random.choice(conv1)
        b = random.choice(A); c = random.choice(A)
        ws = shuffle_words((a,), (b,c))
        if any(len(set(w))<3 for w in ws): continue
        if any(w[0]==('0',) or w[2]==('1',) for w in ws): continue  # keep all convergent
        lhs = tshuffle1({(k,):v for k,v in SE.s1_vec(a, idx, s).items()}, SE.S2_pair(b,c,idx,s))
        rhs = SE.tadd(*[SE.S3(*w, idx, s) for w in ws])
        cnt += 1
        if lhs != rhs:
            ok = False; print('FAIL', s, a, b, c); break
print(f'shuffle-homomorphism exact checks: {cnt} trials,', 'PASS' if ok else 'FAIL')
