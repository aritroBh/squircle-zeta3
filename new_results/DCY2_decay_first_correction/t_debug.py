
from flint import ctx
ctx.dps = 40
import dcy2_cert as C
for N in [400, 2000, 10000, 10**5, 10**6, 10**9, 10**12]:
    r, q = C.ratio(N)
    print(f"N={N:9d} ratio={r.str(10)}  E_mom={q['E_mom'].str(8)}")
    if N == 100000:
        M4=q['M4']; print("   E_a =", (M4/24*q['Gsplit'](4,q['B3'])).str(8))
        print("   alpha4 =", (q['alpha']**4/24*q['Gsplit'](12,q['alpha'])).str(8))
        print("   B6piece =", (q['B6']*(q['B6']*q['Delta']**6).exp()*q['Gsplit'](6,q['alpha'])).str(8))
        L=q['L']; L3,L4,L5=L(3),L(4),L(5); Delta=q['Delta']
        beta2 = L3*L5 + L4**2/2 + L4*L5*Delta + L5**2*Delta**2/2
        beta3 = (3*L3**2*L4 + (3*L3**2*L5+3*L3*L4**2)*Delta + (6*L3*L4*L5+L4**3)*Delta**2
                 + (3*L3*L5**2+3*L4**2*L5)*Delta**3 + 3*L4*L5**2*Delta**4 + L5**3*Delta**5)/6
        print("   beta2*M8 =", (beta2*q['Mtilde'](8)).str(8), " beta3*M10 =", (beta3*q['Mtilde'](10)).str(8))
