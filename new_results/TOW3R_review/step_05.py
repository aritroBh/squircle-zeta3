import sys, os, time, pickle
from fractions import Fraction as Fr
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work+'/tow3_copy')
import sympy as sp
import tow3lib as T3
d, m = T3.d, T3.m
i = sp.I

H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(p,q,k,ddelta):
    return H[(p,q)][ddelta][k]
# QSI arithmetic: value = A + B sqrt2, A=x0+i x2, B=x1+i x3
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def qsi_div(x, y):
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    den = (gmul(C,C)[0]-2*gmul(D,D)[0], gmul(C,C)[1]-2*gmul(D,D)[1])
    numA = (gmul(A,C)[0]-2*gmul(B,D)[0], gmul(A,C)[1]-2*gmul(B,D)[1])
    numB = (gmul(B,C)[0]-gmul(A,D)[0], gmul(B,C)[1]-gmul(A,D)[1])
    return (gdiv(numA, den), gdiv(numB, den))

# 1) symbolic: tow3lib.arho(delta) (Bell from L_parsed) == displayed closed forms
N5 = ( -16*i*d**7 + 304*i*d**6 + (-80-2520*i)*d**5 + 160*d**4*m + (1040+11840*i)*d**4
      - 2080*d**3*m + (-5140-34069*i)*d**3 + (10280-240*i)*d**2*m + (11440+60081*i)*d**2
      + 240*i*d*m**2 + (-22880+960*i)*d*m + (-9600-59920*i)*d - 960*i*m**2 + 19200*m + 25200*i )
a5c = N5/(480*(2*d-5)*(2*d-3))
N6 = ( -64*d**9 + 1728*d**8 + (-20976+480*i)*d**7 - 960*i*d**6*m + (150432-10080*i)*d**6
      + 20160*i*d**5*m + (-702876+89040*i)*d**5 + (-2880-178080*i)*d**4*m + (2221452-423840*i)*d**4
      + 2880*d**3*m**2 + (34560+847680*i)*d**3*m + (-4760069+1146870*i)*d**3
      - 34560*d**2*m**2 + (-138960-2293260*i)*d**2*m + (6677868-1670730*i)*d**2
      + (138960-1440*i)*d*m**2 + (186480+3341460*i)*d*m + (-5524995+1019460*i)*d
      + 960*i*m**3 - 186480*m**2 - 2038920*i*m + 1984500 )
a6c = N6/(5760*(2*d-7)*(2*d-5)*(2*d-3))
a5 = T3.arho(5); a6 = T3.arho(6)
print('a5 symbolic Bell==displayed:', sp.simplify(a5-a5c)==0, flush=True)
print('a6 symbolic Bell==displayed:', sp.simplify(a6-a6c)==0, flush=True)

# 2) data check: u(delta)/u(0) == a_delta (Gaussian rational, sqrt2-part 0), ALL keys, delta=3,4,5,6
t0=time.time()
res = {}
for delta, aa in [(3,None),(4,None),(5,a5),(6,a6)]:
    tested=0; bad=[]
    for key in H:
        if key=='_meta': continue
        p,q=key; dd=p+q
        if 2*dd-4 < delta or delta not in H[key] or dd not in H[key][delta]: continue
        a = aa if aa is not None else T3.arho(delta)
        rv = sp.expand(a.subs({d:dd, m:p}))
        re_, im_ = sp.Rational(sp.re(rv)), sp.Rational(sp.im(rv))
        ratio = qsi_div(u(p,q,dd,delta), u(p,q,dd,0))
        tested += 1
        if ratio[1]!=(Fr(0),Fr(0)) or ratio[0]!=(re_,im_):
            bad.append(key)
    res[delta]=(tested,len(bad),bad[:5])
    print(f'rho({delta}): tested {tested} bad {len(bad)} {bad[:5]}  ({time.time()-t0:.0f}s)', flush=True)
