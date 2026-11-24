"""r02: (i) own float evaluator (rc.evaluate, W5G2R's independent code) on the three new atom vectors, on every atom
vector of R_14 and R_5, on each of the 28 differing classical terms, and on Delta; (ii) exact rebuild of Delta from
my own parse of the two JSON forms with MY shuffle (rc.shv) for products: atom vectors of cross-ratio arguments come
from cmzv.atom (project converter, hypothesis), the three new atoms from delta.pkl; compare with stored Delta exactly;
(iii) own direct quadrature of the three pulled-back iterated integrals in the v-plane (no word conversion) vs mpmath."""
from rc import *
from rc import _gl
import json, re, mpmath as mp
sys.path.insert(0, T + '/A3h4_frame_ds'); sys.path.insert(0, T + '/A3c_claude_elementary/Q')
import cmzv as cm
from qsqrt2 import QS2 as Q
src=open(HERE+'/r01_table.py').read()
ns={'F':F,'re':re}; exec(src[src.index('def jarg'):src.index('J={jkey')], ns); jkey=ns['jkey']
LOG='r02.log'
mp.mp.dps=40; r2=mp.sqrt(2)
t0=time.time()
H=evaluate(); log(f'evaluated {len(H)} admissible words of weight <=5 ({time.time()-t0:.1f}s)',LOG)
def V(d): return sum(complex(float(c))*H[w] for w,c in d.items())
D=pickle.load(open(T+'/W5C14_fourteen_form/delta.pkl','rb'))
NEW=D['new_atoms']
def xval(pq): return mp.mpf(pq[0].numerator)/pq[0].denominator+mp.mpf(pq[1].numerator)/pq[1].denominator*r2
# (i) new atoms
for pq,vec in NEW.items():
    z=V(vec); ref=complex(mp.polylog(5,xval(pq)))
    log(f'new atom Li5({pq[0]}{pq[1]:+}*r2): {len(vec)} words, own value {z.real:.16f}{z.imag:+.1e}i, mpmath {ref.real:.16f}, |diff| {abs(z-ref):.2e}',LOG)
# atom vectors
LCHI={(-1,-1,-1,-1,1):F(1,4),(-1,-1,-1,-1,3):F(-1,4),(-1,-1,-1,-1,5):F(-1,4),(-1,-1,-1,-1,7):F(1,4)}
def atom(k):
    t=k[0]
    if t=='lam': return cm.atom(('log',Q(2)))
    if t=='L': return cm.atom(('log',cm.A))
    if t=='pi2': return cm.atom(('pi2',))
    if t=='z3': return cm.atom(('z3',))
    if t=='z5': return {(-1,-1,-1,-1,0):F(1)}
    if t=='Lchi': return LCHI          # represents L(chi8,5)/sqrt2
    _,n,pq=k
    if n==5 and pq in NEW: return NEW[pq]
    return cm.atom(('li',n,Q(pq[0],pq[1])))
def aval(k):
    t=k[0]
    if t=='lam': return mp.log(2)
    if t=='L': return mp.log(1+r2)
    if t=='pi2': return mp.pi**2
    if t=='z3': return mp.zeta(3)
    if t=='z5': return mp.zeta(5)
    if t=='Lchi': return (mp.zeta(5,mp.mpf(1)/8)+mp.zeta(5,mp.mpf(7)/8)-mp.zeta(5,mp.mpf(3)/8)-mp.zeta(5,mp.mpf(5)/8))/mp.mpf(8)**5/r2
    return mp.polylog(k[1],xval(k[2]))
R14={jkey(k):(F(p),F(q)) for k,(p,q) in json.load(open(T+'/M4f_w5_extras/R5_best.json'))['terms'].items()}
R5 ={jkey(k):(F(p),F(q)) for k,(p,q) in json.load(open(T+'/M4f_w5_extras/R5_min3.json'))['terms'].items()}
log(f'own parse: R14 {len(R14)} terms, R5 {len(R5)} terms',LOG)
atoms=sorted({a for k in list(R14)+list(R5) for a in k},key=str)
worst=0
for a in atoms:
    v=atom(a); e=abs(V(v)-complex(aval(a))); worst=max(worst,e/max(1,abs(complex(aval(a)))))
    if e>1e-12: log(f'  atom {a}: |diff| {e:.2e}  <-- LARGE',LOG)
log(f'all {len(atoms)} atom vectors (incl. new ones, L(chi8,5)/sqrt2): max relative |own value - mpmath| = {worst:.2e}',LOG)
# (ii) rebuild Delta with my own shuffle
cache={}
def term_vec(k):
    if k in cache: return cache[k]
    v=None
    for a in k:
        av={w:F(c) for w,c in atom(a).items()}
        v=av if v is None else shv(v,av)
    cache[k]=v; return v
def form_vec(form):
    rat,sq={},{}
    for k,(cp,cq) in form.items():
        v=term_vec(k)
        if k==(('Lchi',),): cp,cq=2*cq,cp       # (cp+cq r2)*r2*X
        for w,c in v.items():
            rat[w]=rat.get(w,0)+cp*c; sq[w]=sq.get(w,0)+cq*c
    return rat,sq
a14=form_vec(R14); a5=form_vec(R5)
mine=[]
for j in range(2):
    d=dict(a14[j])
    for w,c in a5[j].items(): d[w]=d.get(w,0)-c
    mine.append({w:c for w,c in d.items() if c})
for j,comp in enumerate(('rat','sq')):
    st=D[comp]; bad=sum(1 for w in set(st)|set(mine[j]) if st.get(w,0)!=mine[j].get(w,0))
    log(f'exact rebuild of Delta_{comp} (own parse, own shuffle): support {len(mine[j])}, mismatches with stored delta.pkl: {bad}',LOG)
# shuffle sample: cm.mul vs my shv on a product
k=[k for k in R14 if len(k)==3 and k[0][0]=='Li'][0]
x=atom(k[0]); y=atom(k[1])
log(f'shuffle sample {k[0]} x {k[1]}: cmzv.mul == own shv: {cm.mul(x,y)=={w:F(c) for w,c in shv({w:F(c) for w,c in x.items()},{w:F(c) for w,c in y.items()}).items()}}',LOG)
# differing classical terms: converted vector value vs mpmath value
diff={k:(R14.get(k,(0,0))[0]-R5.get(k,(0,0))[0],R14.get(k,(0,0))[1]-R5.get(k,(0,0))[1]) for k in set(R14)|set(R5)}
diff={k:v for k,v in diff.items() if v!=(0,0)}
log(f'differing classical terms: {len(diff)}; Li5 arguments among them: {sorted(set(a[2] for k in diff for a in k if a[0]=="Li" and a[1]==5 and len(k)==1))}',LOG)
tot=0; totref=mp.mpf(0); wt=0
for k,(cp,cq) in sorted(diff.items(),key=str):
    v=term_vec(k); val=V(v)
    ref=mp.mpf(1)
    for a in k: ref*=aval(a)
    c=xval((F(cp),F(cq)))*(r2 if k==(('Lchi',),) else 1)
    e=abs(val-complex(ref)); wt=max(wt,e/max(1,abs(complex(ref))))
    tot+=complex(c)*val; totref+=c*ref
log(f'each differing term: max relative |own value of converted vector - mpmath| = {wt:.2e}',LOG)
log(f'sum over differing terms: own evaluator {tot.real:.3e}{tot.imag:+.1e}i ; mpmath (R14-R5 classical) {mp.nstr(totref,5)}',LOG)
vr=V(D['rat']); vs=V(D['sq'])
sc=sum(abs(float(c))*abs(H[w]) for w,c in D['rat'].items())
log(f'Delta_rat.V = {vr:.3e}, Delta_sq.V = {vs:.3e}  (scale sum|c||H| = {sc:.3e}; rel {abs(vr)/sc:.1e})',LOG)
# (iii) direct quadrature in the v-plane
a_=1+2**.5; b_=2**.5-1
def itint(forms,c,e,q=30,K=50):
    """I_{c->e}[forms], forms outer->inner, each a callable of v; cumulative Gauss-Legendre on panels dyadic toward both ends"""
    x,wq,Sm=_gl(q)
    br=[0.0]+[2.0**-k for k in range(K,0,-1)]
    br=br+[1-2.0**-k for k in range(2,K+1)]+[1.0]
    br=sorted(set(br))
    n=len(forms); Fv=[0.0]*(n+1); Fv[n]=1.0   # Fv[k] = value of iterated integral of forms[k:] so far
    for p0,p1 in zip(br[:-1],br[1:]):
        h=(p1-p0)/2; tau=(p0+p1)/2+h*x; v=c+(e-c)*tau; dv=(e-c)
        Fn=[None]*(n+1); Fn[n]=np.ones(q)
        new=[0.0]*(n+1); new[n]=1.0
        for k in range(n-1,-1,-1):
            g=forms[k](v,tau)*dv*Fn[k+1]
            Fn[k]=Fv[k]+h*(Sm@g); new[k]=Fv[k]+h*(g@wq)
        Fv=new
    return Fv[0]
dl=lambda p: (lambda v: 1/(v-p))
# forms take (v, tau); near the start v=a of the third path, v-a is computed as (b-a)*tau to avoid cancellation
f_v=lambda v,t: 2/v; f_pm=lambda v,t: 1/(v-a_)+1/(v+a_)
f_pm3=lambda v,t: 1/((b_-a_)*t)+1/(v+a_)
cases=[('b^4',(F(17),F(-12)),[f_v]*4+[lambda v,t:-f_pm(v,t)],0.0,b_),
       ('(4-3r2)/8',(F(1,2),F(-3,8)),[lambda v,t:f_v(v,t)-f_pm(v,t)]*4+[f_pm],0.0,b_),
       ('1-b^4',(F(-16),F(12)),[f_pm3]*4+[lambda v,t:-f_v(v,t)],a_,b_)]
for nm,pq,fs,c,e in cases:
    z=itint(fs,c,e); ref=float(mp.polylog(5,xval(pq)))
    log(f'direct v-plane quadrature, Li5({nm}): {z:.15f}  mpmath {ref:.15f}  |diff| {abs(z-ref):.1e}',LOG)
log(f'done ({time.time()-t0:.0f}s)',LOG)
