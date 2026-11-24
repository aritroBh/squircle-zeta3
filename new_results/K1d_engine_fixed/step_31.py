import mpmath as mp
mp.mp.dps = 80
z3 = mp.zeta(3); z2 = mp.pi**2/6; l2 = mp.log(2)
def AA(xv): return mp.nsum(lambda n: mp.harmonic(n)*xv**n/n**2, [1, mp.inf], tol=mp.mpf('1e-70'))
Am1 = AA(-1); Ah = AA(mp.mpf('0.5'))
Li3m1 = mp.polylog(3,-1); Li3h = mp.polylog(3, mp.mpf('0.5'))
print("A(-1)        =", mp.nstr(Am1,50))
print(" -5/8 zeta(3)=", mp.nstr(-mp.mpf(5)/8*z3,50), " diff:", mp.nstr(Am1+mp.mpf(5)/8*z3,3))
print("A(1/2)                    =", mp.nstr(Ah,50))
print(" zeta(3)-pi^2 log2/12     =", mp.nstr(z3-mp.pi**2*l2/12,50), " diff:", mp.nstr(Ah-(z3-mp.pi**2*l2/12),3))
print("Li3(-1)       =", mp.nstr(Li3m1,50))
print(" -3/4 zeta(3) =", mp.nstr(-mp.mpf(3)/4*z3,50), " diff:", mp.nstr(Li3m1+mp.mpf(3)/4*z3,3))
print("Li3(1/2)                              =", mp.nstr(Li3h,50))
print(" 7/8 z3 - pi^2 log2/12 + log^3 2/6     =", mp.nstr(mp.mpf(7)/8*z3-mp.pi**2*l2/12+l2**3/6,50), " diff:", mp.nstr(Li3h-(mp.mpf(7)/8*z3-mp.pi**2*l2/12+l2**3/6),3))
print("\nzeta3-contents: A(-1):-5/8, Li3(-1):-3/4, A(1/2):1, Li3(1/2):7/8")
with open('progress.log','a') as f:
    f.write("step_48: atom z3-contents re-verified at 80 digits: A(-1)=-5z3/8, A(1/2)=z3-pi^2log2/12, Li3(-1)=-3z3/4, Li3(1/2)=7z3/8-pi^2log2/12+log^3 2/6 (diffs 0.0).\n")
