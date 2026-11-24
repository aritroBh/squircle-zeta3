work='<project>/research/08_theory/K2c_kappa2/'
code = open(work+'step_30.py').read()
code += r'''

print()
print("=== isolate: m=0 ratio and m=1 relative contribution ===")
a = mp.pi/2; Z = 1j*a; q = mp.sqrt(Z); g1 = g1f(a)
for N in [2000, 5000, 12000]:
    sig = 1/mp.sqrt(mpf(N))
    J0_1 = Jm(N, Z, 0, 1); J0_2 = Jm(N, Z, 0, 2)
    lead = 2*(s2-1)*q/sig
    D = phi_p(sig,2,q,Z) - phi_p(sig,1,q,Z)
    S_exp = exp(-D)
    w1 = x_p(sig,1,q)/y_p(sig,1,q)**2
    w2 = x_p(sig,2,q)/y_p(sig,2,q)**2
    Pref = sqrt(phi2der(sig,1,q,Z)/phi2der(sig,2,q,Z))*(x_p(sig,1,q)/x_p(sig,2,q))*(w2/w1)**mp.mpf('1.5')
    F3 = B0(s2*sig,q)/B0(sig,q)
    R0_quad = J0_2/J0_1*exp(lead)
    R0_model = S_exp*Pref*F3*exp(lead)
    print('N=%d m=0 ratio*e^lead: quad=%s model=%s' % (N, mp.nstr(R0_quad,14), mp.nstr(R0_model,14)))
    # m=1 relative contribution within each Iplus
    I1_0 = Iplus_m(N,a,1,0); I1_1 = Iplus_m(N,a,1,1)
    I2_0 = Iplus_m(N,a,2,0); I2_1 = Iplus_m(N,a,2,1)
    Nf = sig**-2
    m1_model = mp.mpf('1.5')*g1*w1/Nf * B1(sig,q)/B0(sig,q)
    m2_model = mp.mpf('1.5')*g1*w2/Nf * B1(s2*sig,q)/B0(s2*sig,q)
    print('   r1(p=1): quad=%s model=%s' % (mp.nstr(I1_1/I1_0,14), mp.nstr(m1_model,14)))
    print('   r1(p=2): quad=%s model=%s' % (mp.nstr(I2_1/I2_0,14), mp.nstr(m2_model,14)))
'''
open(work+'step_31.py','w').write(code)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python', work+'step_31.py'],capture_output=True,text=True,timeout=3600)
print(r.stdout); print(r.stderr[-500:])
