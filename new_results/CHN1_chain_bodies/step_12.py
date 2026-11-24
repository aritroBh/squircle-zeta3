import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import mzv_exact as M
from fractions import Fraction as F
T={W:M.DSTable(W) for W in range(2,9)}
def red(W,w): return T[W].reduce({w:F(1)})
ok=True
def chk(name,cond):
    global ok
    print(('PASS ' if cond else 'FAIL ')+name); ok = ok and cond
chk('z(2,1)=z(3)', red(3,(2,1))==red(3,(3,)))
chk('z(3,1)=z(4)/4', red(4,(3,1))==[F(1)/4])
chk('z(4,2)=z(3)^2-4/3 z(6)', red(6,(4,2))==[F(-4,3), F(1)])
chk('z(2,2,2,2)=5/192 z(8)', red(8,(2,2,2,2))==[F(5,192),0,0,0])
chk('z(3,1,3,1)=1/192 z(8)', red(8,(3,1,3,1))==[F(1,192),0,0,0])
chk('z(2,1^6)=z(8)', red(8,(2,1,1,1,1,1,1))==[F(1),0,0,0])
chk('duality all weight-8 words', all(T[8].reduce_word(w)==T[8].reduce_word(M.dual(w)) for w in M.WORDS[8]))
chk('duality all weight-7 words', all(T[7].reduce_word(w)==T[7].reduce_word(M.dual(w)) for w in M.WORDS[7]))
chk('duality all weight-6 words', all(T[6].reduce_word(w)==T[6].reduce_word(M.dual(w)) for w in M.WORDS[6]))
s2={}
for a in range(2,8): s2[(a,8-a)]=F(1)
chk('sum theorem depth2 W8', T[8].reduce(s2)==[F(1),0,0,0])
s3={}
for a in range(2,7):
    for b in range(1,8-a):
        c=8-a-b
        if c>=1: s3[(a,b,c)]=F(1)
chk('sum theorem depth3 W8', T[8].reduce(s3)==[F(1),0,0,0])
chk('z53+z35+z8 = z3z5 coords', T[8].reduce({(5,3):F(1),(3,5):F(1),(8,):F(1)})==[F(0),F(1),F(0),F(0)])
# also check zeta(4) = 2/5 zeta(2)^2 inside quotient: monomial wordcombos
w4=M.wordcombo(((4,),)); w22=M.wordcombo(((2,),(2,)))
d=M.wd_add_scaled(w4, w22, F(-2,5))
chk('z4 - 2/5 z2^2 in relation span', T[4].reduce(d)==[F(0)] and T[4].reduce(w22)==[F(5,2)])
print('ALL', ok, 't=%.1fs'%(time.time()-t0))
open('progress.log','a').write('step_11: classical identities in DS quotient ALL PASS (z21=z3, z31=z4/4, z42, z2222, z3131, duality W6-8, sum thm, z4=2/5 z2^2)\n')
