import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 1200
cc.N = 800
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 50

def hmp(w):
    hv = cc.h(w)
    return mp.mpc(mp.mpf(hv.real.mid().str(45,radius=False,more=True)), mp.mpf(hv.imag.mid().str(45,radius=False,more=True)))

# Known classical: zeta(1,1,1)=zeta(2,1)=zeta(1,2)=zeta(3)
print('h((0,0,0)) vs zeta3:', mp.nstr(abs(hmp((0,0,0)) - mp.zeta(3)),4))
print('h((0,-1,0)) vs zeta3:', mp.nstr(abs(hmp((0,-1,0)) - mp.zeta(3)),4))
print('h((-1,0,0)) vs zeta3:', mp.nstr(abs(hmp((-1,0,0)) - mp.zeta(3)),4))

# triple sum general: zeta(1,1,1; a,b,c) = sum over n1<n2<n3 a^n1 b^n2 c^n3/(n1 n2 n3)
def d111(e1,e2,e3,M=40000):
    z1,z2,z3 = [mp.exp(mp.pi*mp.j/4*e) for e in (e1,e2,e3)]
    s = mp.mpc(0)
    # iterate n3; need S(n3)=sum_{n1<n2<n3} z1^n1 z2^n2/(n1 n2)
    # accumulate incrementally: when n2 goes from n3-1 to n3: add z2^n2/n2 * P(n2-1), P(m)=sum_{n1<=m} z1^n1/n1
    S = mp.mpc(0); P = mp.mpc(0)
    out = mp.mpc(0)
    for n3 in range(1, M+1):
        # update P to n3-1? careful: first add n2 = n3-1 term to S using P up to n3-2
        n2 = n3-1
        if n2 >= 2:
            S += z2**n2/n2 * P
        P += z1**n3/n3  # P now sum_{n1<=n3}
        if n3 >= 3:
            out += z3**n3/n3 * S
    return out

hv = hmp((2,2,2)); dv = d111(0,0,2)
print('h(2,2,2) vs zeta(1,1,1;1,1,i):', mp.nstr(abs(hv-dv),4), mp.nstr(hv,12), mp.nstr(dv,12))
hv = hmp((1,1,1)); dv = d111(0,0,1)
print('h(1,1,1) vs zeta(1,1,1;1,1,xi8):', mp.nstr(abs(hv-dv),4))
hv = hmp((5,5,5)); dv = d111(0,0,5)
print('h(5,5,5) vs zeta(1,1,1;1,1,-xi8):', mp.nstr(abs(hv-dv),4))
hv = hmp((0,0,2)); dv = d111(0,0,2)
print('h(0,0,2):', mp.nstr(abs(hv-dv),4))
hv = hmp((7,4,0)); dv = d111(4,7,0)
print('h(7,4,0) vs zeta(1,1,1; xi8^4?...):', mp.nstr(abs(hv-dv),4))
# h(7,4,0): standard tuple (0,4,7): eta1=1, eta2=xi8^4=-1, eta3=xi8^7: product telescopes: z1=eta1/eta2, z2=eta2/eta3, z3=eta3
# z1 = 1/-1 = -1 = xi8^4, z2 = -1/xi8^7 = xi8^{4-7}=xi8^-3=xi8^5, z3 = xi8^7
dv2 = d111(4,5,7)
print('h(7,4,0) vs zeta(1,1,1;-1,-xi8,xi8^7):', mp.nstr(abs(hv-dv2),4))
