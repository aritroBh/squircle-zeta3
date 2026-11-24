# verify_pgen1r.py — assembled from executed session steps (not re-run end-to-end)
import sys, itertools, mpmath as mp
from fractions import Fraction as Fr
from flint import fmpz_mat
sys.path.insert(0, '.')
import symbol_engine_r as SE

def taskA():
    for s in (1,2,3,4):
        L, idx = SE.letters(s); n = len(L)
        A = [('0',),('1',)] + [('p',i) for i in range(s)]
        Pts = A + [('inf',)]
        Ws = [SE.S3(*w, idx, s) for w in SE.words_dc(A, s)]
        CR = []
        for p4 in itertools.combinations(Pts, 4): CR += SE.cross_ratios(list(p4), idx, s)
        CRset = {tuple(sorted(c.items())): (c,o) for c,o in CR}
        CRs = [SE.li3_symbol(c,o) for c,o in CRset.values()]
        B = [SE.bracket_tensor(w) for w in SE.lyndon3(n)]
        Mcr = fmpz_mat([SE.pairing_row(T,B) for T in CRs])
        Mw  = fmpz_mat([SE.pairing_row(T,B) for T in Ws])
        Mall= fmpz_mat([[Mcr[i,j] for j in range(Mcr.ncols())] for i in range(Mcr.nrows())]
                       + [[Mw[i,j] for j in range(Mw.ncols())] for i in range(Mw.nrows())])
        assert not any(SE.pairing_row(g,B) for g in SE.shuffle_D(n))
        expect = {1:(6,2,2),2:(30,10,7),3:(90,30,17),4:(210,70,34)}[s]
        assert (len(CRset), Mcr.rank(), Mw.rank()) == expect and Mall.rank()==Mcr.rank(), s
        print(f'A s={s}: #CR={len(CRset)} rankCR={Mcr.rank()} rankW={Mw.rank()} inclusion OK')

def taskCD():
    mp.mp.dps = 250
    phi = (1+mp.sqrt(5))/2; s5 = mp.sqrt(5); z3 = mp.zeta(3)
    Li3 = lambda x: mp.re(mp.polylog(3, mp.mpf(x)))
    Li2 = lambda x: mp.re(mp.polylog(2, mp.mpf(x)))
    L3  = lambda x: Li3(x) - mp.log(abs(x))*Li2(x) - mp.log(abs(x))**2*mp.log(1-x)/3
    L5  = (mp.zeta(3,mp.mpf(1)/5)-mp.zeta(3,mp.mpf(2)/5)-mp.zeta(3,mp.mpf(3)/5)+mp.zeta(3,mp.mpf(4)/5))/125
    D = (Li3(1/phi)-Li3(-1/phi))/2 - ((25*s5/48)*L5 - mp.pi**2*mp.log(phi)/12 + mp.log(phi)**3/4)
    assert abs(D) < mp.mpf('1e-200'); print('D residual:', mp.nstr(D,3))
    c1=(Fr(-75,32),Fr(81,64)); c2=(Fr(-35,32),Fr(53,64)); c3=(Fr(33,64),Fr(-19,80)); cl=(Fr(5,32),Fr(-1,8))
    ev = lambda pq,v: mp.mpf(pq[0])*v + s5*mp.mpf(pq[1])*v
    N3 = ev(c1,L3(1/phi)) + ev(c2,L3(-1/phi)) + ev((c3[0]-cl[0],c3[1]-cl[1]), z3)
    rel = mp.pslq([N3, z3, s5*z3, L5, s5*L5], tol=mp.mpf('1e-200'), maxcoeff=10**12, maxsteps=100000)
    assert rel == [3840,-60,-372,-4375,2500]; print('C: pair (1/64+31s5/320, 875/768-125s5/192) confirmed')

def taskA_repair():
    s = 2; L, idx = SE.letters(s)
    p1={idx[('p',0)]:1}; p2={idx[('p',1)]:1}; q1={idx[('q',0)]:1}; q2={idx[('q',1)]:1}; r12={idx[('r',0,1)]:1}
    c  = SE.vadd(p1, r12, -1); cX = SE.vadd(q1, r12, -1)
    sN = SE.tadd(SE.li3_symbol(c, SE.vadd(p2,r12,-1)), {k:-v for k,v in SE.li3_symbol(cX, SE.vadd(q2,r12,-1)).items()})
    A = [('0',),('1',),('p',0),('p',1),('inf',)]
    CR = []
    for p4 in itertools.combinations(A,4): CR += SE.cross_ratios(list(p4), idx, s)
    CRs = [SE.li3_symbol(cc,o) for cc,o in {tuple(sorted(cc.items())):(cc,o) for cc,o in CR}.values()]
    B = [SE.bracket_tensor(w) for w in SE.lyndon3(len(L))]
    Mcr = fmpz_mat([SE.pairing_row(T,B) for T in CRs])
    Mall= fmpz_mat([[Mcr[i,j] for j in range(Mcr.ncols())] for i in range(Mcr.nrows())]
                   + [SE.pairing_row({k:2*v for k,v in sN.items()}, B)])
    assert Mall.rank()==Mcr.rank()==10; print('A repair: repeated-word Nielsen symbol in span s(CR)+D')

taskA(); taskCD(); taskA_repair()
print('verify_pgen1r: all checks passed')
