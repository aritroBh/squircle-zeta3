import sys, os
WORK = 'research/08_theory/A3f4_glanois_basis'
os.chdir(WORK)
sys.path.insert(0, '<home>)/paper/squircle_note/scripts')
sys.path.insert(0, os.path.join(WORK,'Q'))
import check_conversion as cc
from flint import ctx, acb, arb
ctx.prec = 3000
cc.N = 2500
cc.roots = [acb(0,arb.pi()*j/4).exp() for j in range(8)]
cc.roots[0],cc.roots[2],cc.roots[4],cc.roots[6] = acb(1),acb(0,1),acb(-1),acb(0,-1)
import mpmath as mp
mp.mp.dps = 40

def hmp_acb(z): return mp.mpc(mp.mpf(z.real.mid().str(35,radius=False,more=True)), mp.mpf(z.imag.mid().str(35,radius=False,more=True)))

# full g_at_one((-1,-1,6)) vs quad:
# G = int_{0<t1<t2<t3<1} [dt1/(t1+i)] [dt2/t2] [dt3/t3]
G = mp.quad(lambda t3: (1/t3)*mp.quad(lambda t2: (1/t2)*mp.quad(lambda t1: 1/(t1+1j), [0,t2]), [0,t3]), [0,1])
print('quad G:', mp.nstr(G,20))
g = cc.g_at_one((-1,-1,6))
print('g_at_one:', mp.nstr(hmp_acb(g),20))
print('h(-1,-1,2) = -G:', mp.nstr(hmp_acb(cc.h((-1,-1,2))),20))
print('-Li3(i):', mp.nstr(-mp.nsum(lambda n:(1j)**n/n**3,[1,mp.inf]),20))
# individual halves:
# half((-1,-1,6),0) = G0((-1,-1,6);1/2): int_{0<...<1/2}
H = mp.quad(lambda t3: (1/t3)*mp.quad(lambda t2: (1/t2)*mp.quad(lambda t1: 1/(t1+1j), [0,t2]), [0,t3]), [0,mp.mpf('0.5')])
print('half((-1,-1,6),0) vs quad:', mp.nstr(hmp_acb(cc.half((-1,-1,6),0)),18), mp.nstr(H,18))
# half((6,),0) = int_0^{1/2} dt/(t+i)
H1 = mp.quad(lambda t: 1/(t+1j), [0,mp.mpf('0.5')])
print('half((6,),0) vs quad:', mp.nstr(hmp_acb(cc.half((6,),0)),18), mp.nstr(H1,18))
# half((8,),1): letter 8, kind1 pole=1: int_0^{1/2} dt/(t-1)
H2 = mp.quad(lambda t: 1/(t-1), [0,mp.mpf('0.5')])
print('half((8,),1) vs quad:', mp.nstr(hmp_acb(cc.half((8,),1)),18), mp.nstr(H2,18))
# half((8,8),1): int_{0<t1<t2<1/2} dt1/(t1-1) dt2/(t2-1)  [innermost w[-1]=8]
H3 = mp.quad(lambda t2: 1/(t2-1)*mp.quad(lambda t1: 1/(t1-1), [0,t2]), [0,mp.mpf('0.5')])
print('half((8,8),1) vs quad:', mp.nstr(hmp_acb(cc.half((8,8),1)),18), mp.nstr(H3,18))
# half((6,8,8),1): poles: letter6 kind1 = 1-roots[6] = 1+i; letter 8 pole 1:
# innermost w[-1]=8 pole 1, then 8, then 6 outermost (w[0]):
# int_{0<t1<t2<t3<1/2} dt1/(t1-1) dt2/(t2-1) dt3/(t3-(1+i))
H4 = mp.quad(lambda t3: 1/(t3-(1+1j))*mp.quad(lambda t2: 1/(t2-1)*mp.quad(lambda t1: 1/(t1-1), [0,t2]), [0,t3]), [0,mp.mpf('0.5')])
print('half((6,8,8),1) vs quad:', mp.nstr(hmp_acb(cc.half((6,8,8),1)),18), mp.nstr(H4,18))
