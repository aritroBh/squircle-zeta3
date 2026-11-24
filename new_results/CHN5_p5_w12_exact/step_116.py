import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
with open('progress.log','a') as f:
    f.write(time.strftime('%H:%M:%S ')+"STEP1 forms finalized: A=(1/100)a2a3a4(Sum a)*P8_A (367 terms, +243/-124, P8(1^5)=65888), B=(1/224)*P8_B (359), CLASS=(1/432)*P8_CLASS (450), Z=(1/3)*P8_Z (226). ALL integer-primitive, IRREDUCIBLE / Q, reversal-symmetric (monomial check, simultaneous=True). Saved P8_{A,B,CLASS,Z}.txt.\n")
    f.write(time.strftime('%H:%M:%S ')+"Resultants: Res_a5(C4,P7)=a1^4(a1+1)^4*F3*F3*F8*F6 [14s], Res_a5(C4,QA)=a1(a1+1)*G30 (deg30, 4243 terms) [143s]. Traced Gamma pts lie on G22 (resid ~1e-128) + F6(a1=0.01 br)/F8(larger-a1 brs). G30 != 0 at all 119 traced pts (min 6e4) => QA != 0 there (resultant argument). Exact A reproduces engine at refined pts (A=-46.7202620046 at a1=0.01; B/A=75/14, Z/A=-25/3 to 150 digits).\n")
    f.write(time.strftime('%H:%M:%S ')+"Re-verified A,B,CLASS,Z = L*Q8 vs engine at 3 fresh rational pts (exact match). VERDICT: NO irreducible-free point on Gamma+ through W12: A<0 and CLASS>0 at every interior Gamma pt (NUMERICAL, multi-method, a1/a3 1e-11..1e9+); A,CLASS span V10 W12 content (EXACT).\n")
print("progress.log finalized")
# summary.md final
summary = open('summary.md').read()
summary += """

## Step 3 extra detail (resultant structure)
- Res_a5(C4,P7) = a1^4 (a1+1)^4 F3 F3 F8 F6 ; Res_a5(C4,QA) = a1 (a1+1) G30 (G30 deg 30, 4243 terms).
- QA=0 at a Gamma pt => G30=0 at its (a1,a2,a4) projection. G30 != 0 at all 119 traced pts (min |G30|=6e4).
- Traced pts on G22 (resid ~1e-128) + F6 (small-a1 branch) / F8 (larger branches).
- Exact poly A matches engine at refined pts (150 digits; B/A=75/14, Z/A=-25/3 exact).
- Fresh-point re-verification of all four identities vs engine: exact match at 3 random rational pts.

## Remaining gap
Full exact elimination for {C4,Q7',P7,P8_A} (a1-projection degree ~10^3+) is infeasible with
sympy/flint here (Sylvester sizes 40x40+ with polynomial entries, final univariate degree ~10^5-10^6).
Step-3 sign statement is therefore NUMERICAL (very strong: 119+ interior pts, all branches found,
a1/a3 in 1e-11..1e9+, asymptotics clean, independent cross-checks, resultant condition G30!=0).
"""
open('summary.md','w').write(summary)
print("summary.md finalized")
print("files:", sorted([f for f in os.listdir('.') if f.endswith('.pkl') or f.endswith('.txt') or f.endswith('.log') or f.startswith('run_')]))
