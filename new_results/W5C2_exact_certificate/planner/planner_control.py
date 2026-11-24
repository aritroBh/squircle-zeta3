# Control for planner_verify_w5c2 (B): a perturbed frame-w5 row (one coefficient +1) must NOT vanish numerically.
import sys, pickle
sys.path.insert(0,'../A3f6_level16'); sys.path.insert(0,'../A3c_claude_elementary/Q'); sys.path.insert(0,'<project>/paper/squircle_note/scripts')
import mpmath as mp
from flint import ctx, acb, arb
import eval8 as ev
ctx.prec=400; ctx.cap=900; ev.N=800
ev.roots=[acb(0,arb.pi()*j/4).exp() for j in range(8)]
ev.roots[0],ev.roots[2],ev.roots[4],ev.roots[6]=acb(1),acb(0,1),acb(-1),acb(0,-1)
mp.mp.dps=100
row=pickle.load(open('validated_w5.pkl','rb'))['rows'][5]
act=lambda w,u: tuple(-1 if x==-1 else (u*x)%8 for x in w)
for u in (1,3):
    r=[(act(w,u),c) for w,c in row]
    v=lambda rr: abs(mp.fsum(mp.mpf(c)*mp.mpc(*(mp.mpf(p.mid().str(95,radius=False)) for p in (ev.h(w).real,ev.h(w).imag))) for w,c in rr))
    print(u,'row',mp.nstr(v(r),3),'perturbed',mp.nstr(v(r[:1]+[(r[0][0],1)]+r[1:]),3))
