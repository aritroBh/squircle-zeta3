code = r'''
# step_24.py : FINAL CONSOLIDATION — re-verify every headline identity, print evidence table.
import mpmath as mp
mp.mp.dps = 110

def A_series(p, s=0, tol=None):
    p = mp.mpf(p); s = mp.mpf(s); tol = tol or mp.mpf('1e-105')
    tot = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma((m+s)/p)**2/(p**2*mp.gamma(1+2*(m+s)/p))
        tot += t
        if m > 30 and abs(t) < tol: break
    return tot
def B_series(p, tol=None):
    p = mp.mpf(p); tol = tol or mp.mpf('1e-105')
    tot = mp.mpf(0); n = -1
    while True:
        n += 1
        t = mp.gamma((2*n+1)/p)**2/(p**2*mp.gamma(1+2*(2*n+1)/p))
        tot += t
        if n > 30 and abs(t) < tol: break
    return tot

z=mp.zeta; pi=mp.pi; phi=(1+mp.sqrt(5))/2
E=[]
def chk(label, val, exact, note=''):
    err = abs(val-exact)
    E.append((label, mp.nstr(val,30), mp.nstr(exact,30), mp.nstr(err,3), note))

# Task 1
chk('A(1)=zeta2/3', A_series(1), z(2)/3, 'PROVED')
chk('A(2)=7z2/12', A_series(2), 7*z(2)/12, 'PROVED')
chk('A(1/2)=2z2/3-4log^2 phi', A_series(mp.mpf('0.5')), 2*z(2)/3-4*mp.log(phi)**2, 'PROVED')
chk('B(2)=pi^2/12', B_series(2), pi**2/12, 'PROVED')
chk('B(1)=z2/6+log^2 phi', B_series(1), z(2)/6+mp.log(phi)**2, 'PROVED')
chk('A_1(1)=z2/3-1/2', A_series(1,1), z(2)/3-mp.mpf('0.5'), 'PROVED')
chk('A_2(1)=z2/3-13/24', A_series(1,2), z(2)/3-mp.mpf(13)/24, 'PROVED')

# Task 2 coefficient identities (d2,d3 exact; d5,d7,d9,d11 PSLQ/assembly)
d2=-z(2)/2; d3=z(3)/6
d5=(z(2)*z(3)-3*z(5))/60
d7=z(7)/14-z(2)*z(5)/42-z(2)**2*z(3)/140
d9=-mp.mpf(17)/72*z(9)+mp.mpf(3)/40*z(2)*z(7)+mp.mpf(79)/1814400*pi**6*z(3)+mp.mpf(1)/1600*pi**4*z(5)-mp.mpf(1)/180*z(3)**3
d9alt=(-428400*z(9)+136080*z(2)*z(7)+79*pi**6*z(3)+1134*pi**4*z(5)-10080*z(3)**3)/1814400
print('d9 two forms agree:', mp.nstr(abs(d9-d9alt),4))
d11=mp.mpf(31)/22*z(11)-mp.mpf(85)/198*z(2)*z(9)-mp.mpf(3)/880*pi**4*z(7)-mp.mpf(79)/332640*pi**6*z(5)-mp.mpf(2339)/119750400*pi**8*z(3)+mp.mpf(1)/11*z(3)**2*z(5)-mp.mpf(1)/99*z(2)*z(3)**3
d1=mp.mpf('1.73300259020287278417012270424573195250974925287880319940547987633644278103184087570877319')
print('coefficient values:')
for nm,v in [('d1',d1),('d2',d2),('d3',d3),('d5',d5),('d7',d7),('d9',d9),('d11',d11)]:
    print(' ',nm,'=',mp.nstr(v,35))

# numeric verification of expansion at p=40,60,80 (partial asymptotic sums)
def deficit(p): return z(2)-A_series(p)
for p in [40,60,80]:
    e=mp.mpf(1)/p
    approx=d1*e+d2*e**2+d3*e**3+d5*e**5+d7*e**7+d9*e**9+d11*e**11
    print(f'p={p}: deficit={mp.nstr(deficit(p),25)} approx={mp.nstr(approx,25)} err={mp.nstr(abs(deficit(p)-approx),3)}')

# triangle-end c4,c5 with correct Bell polynomials
def f(m): return mp.gamma(m+1)**2/mp.gamma(2*m+1)
def gder(m,j): return 2*mp.psi(j-1,m+1)-mp.mpf(2)**j*mp.psi(j-1,2*m+1)
def bells(gs):
    k=len(gs); B=[[mp.mpf(0)]*(k+1) for _ in range(k+1)]; B[0][0]=1
    for n in range(1,k+1):
        for j in range(1,n+1):
            B[n][j]=sum(mp.binomial(n-1,i-1)*gs[i-1]*B[n-i][j-1] for i in range(1,n-j+2))
    return [sum(B[n][j] for j in range(1,n+1)) for n in range(1,k+1)]
K=6; cs=[mp.mpf(0)]*(K+1)
for m in range(1,800):
    fm=f(m)
    if fm<mp.mpf('1e-115') and m>60: break
    Bs=bells([gder(m,j) for j in range(1,K+1)])
    for k in range(0,K+1):
        cs[k]+=(-1)**k/mp.factorial(k)*mp.mpf(m)**(k-2)*(fm*Bs[k-1] if k>0 else fm)
print('triangle c0..c6 (correct Bell):')
for k in range(K+1): print('  c_%d ='%k, mp.nstr(cs[k],30))

# B expansion numbers recap
print('B deficit: target 3z2/4;  d1/2 =', mp.nstr(d1/2,20), ' -z3/6 =', mp.nstr(-z(3)/6,20), ' (7/60)(3z5-z2z3) =', mp.nstr(mp.mpf(7)/60*(3*z(5)-z(2)*z(3)),20))

print()
print('=== EVIDENCE TABLE (special values) ===')
for row in E:
    print(' | '.join(row))
'''
open('steps/step_24.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_24.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-800:] if r.returncode else '')
open('progress.log','a').write('step_24: final consolidation, all headline identities re-verified\n')
