import mpmath as mp, re
D="<project>/research/08_theory/D13b_Mk_constant"
def load(fn="mk_deg9.tsv"):
    out={}
    for line in open(D+"/"+fn):
        p=line.rstrip("\n").split("\t")
        if not p[0].isdigit(): continue
        mid=re.match(r"\[([-0-9.eE+]+)",p[1]).group(1)
        out[int(p[0])]=mp.mpf(mid)
    return out
def B(z):  # B(z)=1+e^{-z}-sqrt(1+e^{-2z})
    return 1+mp.exp(-z)-mp.sqrt(1+mp.exp(-2*z))
def H(sig): # int_0^inf sqrt(1-e^{-2t}) e^{-sig t} dt = (1/2) Beta(sig/2,3/2)
    return mp.beta(sig/2, mp.mpf(3)/2)/2
def d_cut(j, m, alpha=mp.pi/4, allcuts=True):
    """cut representation of d_j^(m)=int_0^inf B(z) e^{-mz} L_j((m+1)z) dz, j>=1"""
    q=m+1; w=mp.expjpi(alpha/mp.pi)
    def f(r):
        s=r*w
        th=q*mp.pi*s/2
        fac = 1/(2*mp.cos(th)) if allcuts else mp.expj(th)
        return (s/(1+s))**j * H(q*s) * fac /(1+s) * w
    rs=mp.sqrt(2*(j+1)/(mp.pi*q))
    pts=[0,rs/8,rs/2,rs,2*rs,4*rs,10*rs,mp.inf]
    I=mp.quad(f,pts)
    return -(2/mp.pi)*I.imag
