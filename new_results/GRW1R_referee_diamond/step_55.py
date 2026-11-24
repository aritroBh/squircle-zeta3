import mpmath as mp, json
mp.mp.dps = 50
# GRW1's intermediate claim: Gamma(n-5/2)/Gamma(n+1) <= 1.63 n^{-7/2} for n >= 10 -- check n=10
n=10
lhs = mp.gamma(n-mp.mpf('2.5'))/mp.gamma(n+1)
rhs = mp.mpf('1.63')*mp.mpf(n)**mp.mpf('-3.5')
print(f"Gamma(7.5)/Gamma(11) = {mp.nstr(lhs,15)}")
print(f"1.63 * 10^-3.5       = {mp.nstr(rhs,15)}")
print(f"ratio Gamma/(n^-3.5) = {mp.nstr(lhs*n**mp.mpf('3.5'),10)}  > 1.63 : {lhs*n**mp.mpf('3.5')>mp.mpf('1.63')}")
# with the corrected factor, E-term at n=10:
Eterm = (4/mp.pi)*(mp.mpf(5)/24)*(mp.pi/2)**mp.mpf('2.5')*lhs
print("E-term bound at n=10 in units rho^-10 10^-3.5:", mp.nstr(Eterm*(mp.pi/2)**10*10**mp.mpf('3.5'),10))
# sign pattern of f_n vs law for n=2..60
f450=[mp.mpf(v) for v in json.load(open('../D11b2_stokes_jump/step2_num450_N180.json'))['f450']]
pat={0:-1,1:+1,2:+1,3:-1}
bad=[n for n in range(2,60) if (f450[n]>0)*2-1 != pat[n%4] if f450[n]!=0]
print("n in 2..59 where sign(f_n) != law pattern:", bad)
with open('progress.log','a') as f:
    f.write("step_33: 1.63-intermediate falsified at n=10; sign pattern holds for n>=6 (deviations listed)\n")
