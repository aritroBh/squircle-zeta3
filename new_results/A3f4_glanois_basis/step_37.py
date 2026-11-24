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
mp.mp.dps = 30

# direct quadrature of iterated integral for word (2,2,2): innermost-first order (2,2,2)
# I = int_{0<t1<t2<t3<1} om(t1) om(t2) om(t3), om = i dt/(1 - i t)
om = lambda t: 1j/(1-1j*t)
f = lambda t3: om(t3) * mp.quad(lambda t2: om(t2)*mp.quad(lambda t1: om(t1), [0, t2]), [0, t3])
I = mp.quad(f, [0, 1])
print('quad I(2,2,2):', mp.nstr(I, 15))

def hmp(w):
    hv = cc.h(w)
    return mp.mpc(mp.mpf(hv.real.mid().str(25,radius=False,more=True)), mp.mpf(hv.imag.mid().str(25,radius=False,more=True)))
print('h(2,2,2)      :', mp.nstr(hmp((2,2,2)), 15))

# test half() pieces for w=(6,6) kind 0: G_0((6,6);1/2): int_{0<t1<t2<1/2} dt1/(t1+i)dt2/(t2+i)
omg = lambda t: 1/(t+1j)
G = mp.quad(lambda t2: omg(t2)*mp.quad(lambda t1: omg(t1), [0,t2]), [0, mp.mpf('0.5')])
hf = cc.half((6,6),0)
print('G0((6,6);1/2) quad:', mp.nstr(G,15), ' half:', mp.nstr(mp.mpc(mp.mpf(hf.real.mid().str(20,radius=False,more=True)),mp.mpf(hf.imag.mid().str(20,radius=False,more=True))),15))
# kind 1: poles 1-(-i) = 1+i: G_1?? what is half((6,6),1) supposed to be?
omg1 = lambda t: 1/(t-(1+1j))
G1 = mp.quad(lambda t2: omg1(t2)*mp.quad(lambda t1: omg1(t1), [0,t2]), [0, mp.mpf('0.5')])
hf1 = cc.half((6,6),1)
print('kind1 G(1+i;1/2) quad:', mp.nstr(G1,15), ' half:', mp.nstr(mp.mpc(mp.mpf(hf1.real.mid().str(20,radius=False,more=True)),mp.mpf(hf1.imag.mid().str(20,radius=False,more=True))),15))
