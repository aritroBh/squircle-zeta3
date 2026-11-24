import sys, mpmath as mp, pickle, time
sys.path.insert(0,'.')
mp.mp.dps=70
import engine14 as E
t0=time.time()
mons14 = pickle.load(open('mons14.pkl','rb'))
print("assembling over", len(mons14), "monomials")
INTC={}
def INTc(n,a,b,wc):
    k=(n,a,b,wc)
    if k not in INTC: INTC[k]=E.INT(n,a,b,wc)
    return INTC[k]
# distinct types
types = set((n,ab,wc) for (n,ab,wc,fk) in mons14)
print("distinct (n,a,b,wc) types:", len(types))
S={}
for (n,ab,wc,fk),coef in mons14.items():
    J = INTc(n,ab[0],ab[1],wc)
    # letter monomial
    lm={}
    if fk[0]: lm['l2']=fk[0]
    if fk[1]: lm['la']=fk[1]
    if fk[2]: lm['lu']=fk[2]
    lm=E.mkey(lm)
    for mon,c in J.items():
        m2=E.mdot(mon,lm)
        E.e_addinto(S, m2, E.fmul(c,coef))
print("S assembled: %d terms, %.1fs, distinct INT computed: %d"%(len(S),time.time()-t0,len(INTC)))
pickle.dump(S, open('S14.pkl','wb'))
# evaluate
NUM={k:mp.mpf(v) for k,v in pickle.load(open('atoms_num.pkl','rb')).items()}
def e_eval(P):
    tot=mp.mpf(0)
    for mon,c in P.items():
        v=E.fev(c,70)
        for (s,ex) in mon: v*=NUM[s]**ex
        tot+=v
    return tot
Snum = e_eval(S)
LL = mp.sqrt(E.fev(E._g['L'],70))
P3 = 84*LL*Snum
print("S_num =", mp.nstr(Snum,30))
print("84*L*S =", mp.nstr(P3,40))
print("P3 truth=", mp.nstr(mp.mpf('-5.508083698053478283972318931980510244285194094582'),40))
print("diff =", mp.nstr(abs(P3-mp.mpf('-5.508083698053478283972318931980510244285194094582')),5))
open('progress.log','a').write('PG6 s_B6: assembled S (%d terms); 84LS vs P3 truth diff = %s\n'%(len(S),mp.nstr(abs(P3-mp.mpf('-5.508083698053478283972318931980510244285194094582')),3)))
