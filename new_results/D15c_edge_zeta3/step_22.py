
allcoeff = '''
import sympy as sp
import sym_G
# rebuild D3 with all coefficients accessible
G21,G12,syms = sym_G.build()
# Reconstruct D3 by re-running build internals: easier to recompute D3 here
s0,s1,c0,c1 = syms['s0'],syms['s1'],syms['c0'],syms['c1']
L0,L1 = syms['L0'],syms['L1']
tau1,tau2,tau3 = syms['tau1'],syms['tau2'],syms['tau3']
gam1,gam2,gam3 = syms['gam1'],syms['gam2'],syms['gam3']
S0,Sc,lS0 = syms['S0'],syms['Sc'],syms['lS0']
T0=S0
T1=tau1+s0*L0+s1*L1; T2=tau2+s0*L0**2+s1*L1**2; T3=tau3+s0*L0**3+s1*L1**3
Tp0=Sc
Tp1=gam1+c0*(L0+1)+c1*(L1+1); Tp2=gam2+c0*(L0**2+2*L0)+c1*(L1**2+2*L1); Tp3=gam3+c0*(L0**3+3*L0**2)+c1*(L1**3+3*L1**2)
T=[T0,T1,T2,T3]; Tp=[Tp0,Tp1,Tp2,Tp3]; S0p=Tp0
U1=T[1]/S0; U2=T[2]/(2*S0)-T[1]**2/(2*S0**2); U3=T[3]/(6*S0)-T[1]*T[2]/(2*S0**2)+T[1]**3/(3*S0**3)
U1p=Tp[1]/S0-T[1]*S0p/S0**2
U2p=Tp[2]/(2*S0)-T[2]*S0p/(2*S0**2)-T[1]*Tp[1]/S0**2+T[1]**2*S0p/S0**3
U3p=(Tp[3]/(6*S0)-T[3]*S0p/(6*S0**2)-(Tp[1]*T[2]/(2*S0**2)+T[1]*Tp[2]/(2*S0**2)-T[1]*T[2]*S0p/S0**3)+(T[1]**2*Tp[1]/S0**3-T[1]**3*S0p/S0**4))
L=S0p/S0
R1=lS0-U1; R2=-lS0+U1-U2; R3=lS0-U1+U2-U3
Q0=-L; Q1=L-U1p; Q2=-L+U1p-U2p; Q3=L-U1p+U2p-U3p
rho1=R1; rho2=R2+R1**2/2; rho3=R3+R1*R2+R1**3/6
W1=2*Q0*Q1; W2=2*Q0*Q2+Q1**2; W3=2*Q0*Q3+2*Q1*Q2
B=1+Q0**2
E1=W1/(2*B); E2=W2/(2*B)-W1**2/(8*B**2); E3=W3/(2*B)-W1*W2/(4*B**2)+W1**3/(16*B**3)
D0=(1/S0)*sp.sqrt(B)
D3=sp.expand(D0*(rho3+rho2*E1+rho1*E2+E3))
# Extract all coefficients c_{ij} for i+j<=3
coeffs={}
for i in range(4):
    for j in range(4):
        if i+j<=3:
            coeffs[(i,j)] = sp.expand(D3).coeff(L0,i).coeff(L1,j)
# Check: does c30, c03, c20, c02, c11 depend on tau/gam?
logvars=[tau1,tau2,tau3,gam1,gam2,gam3,lS0]
for k in [(3,0),(0,3),(2,0),(0,2),(1,1),(1,0),(0,1),(0,0)]:
    cc = coeffs[k]
    P = sp.Poly(cc, *logvars)
    nonzero=[mm for mm in P.monoms() if any(e!=0 for e in mm)]
    print(f"c{k}: #terms={len(sp.expand(cc).as_ordered_terms())}, depends on tau/gam/lS0: {len(nonzero)>0}", nonzero[:6])
import pickle
pickle.dump({k:sp.srepr(v) for k,v in coeffs.items()}, open('allcoeffs.pkl','wb'))
print("saved allcoeffs.pkl")
'''
with open('step_11_allcoeff.py','w') as f: f.write(allcoeff)
import subprocess
r = subprocess.run(['./.venv/bin/python','step_11_allcoeff.py'],capture_output=True,text=True,timeout=300)
print(r.stdout); print(r.stderr[-2000:] if r.stderr else "")
