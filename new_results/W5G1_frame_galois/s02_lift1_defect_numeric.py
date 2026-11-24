"""s02: numerical check that u=3,5 conjugates of lift1 rows built on weight-1 branch relations are NOT relations,
while u=1,7 are. Uses the same evaluator (eval8) as W5C1/planner."""
import sys, pickle
T='<project>/research/08_theory'
sys.path.insert(0,T+'/A3f6_level16'); sys.path.insert(0,T+'/A3c_claude_elementary/Q'); sys.path.insert(0,'<project>/paper/squircle_note/scripts')
import mpmath as mp
from flint import ctx, acb, arb
import eval8 as ev
from relations import words
ctx.prec=300; ctx.cap=700; ev.N=600
ev.roots=[acb(0,arb.pi()*j/4).exp() for j in range(8)]
ev.roots[0],ev.roots[2],ev.roots[4],ev.roots[6]=acb(1),acb(0,1),acb(-1),acb(0,-1)
mp.mp.dps=60
W5=words(5)
lift1=pickle.load(open(T+'/W5C1_certificate_scoping/rels_lift1.pkl','rb'))
used=pickle.load(open('used_rows.pkl','rb'))
act=lambda w,u: tuple(-1 if x==-1 else (u*x)%8 for x in w)
def val(row):
    return mp.fsum(mp.mpf(c)*mp.mpc(*(mp.mpf(p.mid().str(55,radius=False)) for p in (ev.h(w).real,ev.h(w).imag))) for w,c in row)
seen=0
for tag in ('pp','pm'):
    for f,i,lam,yr,ys in used[tag]:
        if f!='lift1' or lift1[i][0][1] not in (0,1): continue
        row=[(W5[w],c) for w,c in lift1[i][1]]
        print(tag,'lift1 idx',i,'label',lift1[i][0],' '.join(f'u={u}: {mp.nstr(abs(val([(act(w,u),c) for w,c in row])),3)}' for u in (1,7,3,5)),flush=True)
        seen+=1
        if seen>=4: break
