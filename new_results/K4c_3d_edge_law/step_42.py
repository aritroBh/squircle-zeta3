
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 50

def pslq_debug2(x, maxsteps=12):
    n = len(x)
    gam = mp.sqrt(mpf(4)/3)
    nx = mp.sqrt(sum(xi**2 for xi in x))
    y = [xi/nx for xi in x]
    A = [[mpf(1) if i==j else mpf(0) for j in range(n)] for i in range(n)]
    B = [[mpf(1) if i==j else mpf(0) for j in range(n)] for i in range(n)]
    s = [mpf(0)]*n
    for k in range(n-1, -1, -1):
        s[k] = mp.sqrt(sum(y[i]**2 for i in range(k, n)))
    H = [[mpf(0)]*(n-1) for _ in range(n)]
    for j in range(n-1):
        for i in range(n):
            if i < j: H[i][j] = mpf(0)
            elif i == j: H[i][j] = s[j+1]/s[j]
            else: H[i][j] = -y[i]*y[j]/(s[j]*s[j+1])
    def checkYH(tag):
        for j in range(n-1):
            d = sum(y[i]*H[i][j] for i in range(n))
            if abs(d) > mpf('1e-40'):
                print(f"   !! y^T H col{j} = {mp.nstr(d,3)} at {tag}")
    for step in range(maxsteps):
        best = -1; r = 0
        for i in range(n-1):
            val = gam**i * abs(H[i][i])
            if val > best: best = val; r = i
        y[r], y[r+1] = y[r+1], y[r]
        H[r], H[r+1] = H[r+1], H[r]
        A[r], A[r+1] = A[r+1], A[r]
        for k in range(n):
            B[k][r], B[k][r+1] = B[k][r+1], B[k][r]
        if r < n-2:
            t0 = mp.sqrt(H[r][r]**2 + H[r][r+1]**2)
            if t0 == 0:
                print(f"step {step}: t0=0; y={[mp.nstr(t,4) for t in y]}"); break
            t1 = H[r][r]/t0; t2 = H[r][r+1]/t0
            for i in range(r, n):
                t3, t4 = H[i][r], H[i][r+1]
                H[i][r] = t1*t3 + t2*t4
                H[i][r+1] = -t2*t3 + t1*t4
        checkYH(f"after corner step{step}")
        for i in range(r+1, n):
            for j in range(min(i, n-2), -1, -1):
                if H[j][j] == 0: continue
                t = mp.nint(H[i][j]/H[j][j])
                if t != 0:
                    y[j] += t*y[i]
                    for k in range(j+1):
                        H[i][k] -= t*H[j][k]
                    for k in range(n):
                        A[i][k] -= t*A[j][k]
                        B[k][j] += t*B[k][i]
        checkYH(f"after reduce step{step}")
        print(f"step {step}: r={r} y={[mp.nstr(t,4) for t in y]}  Hdiag={[mp.nstr(H[i][i],4) for i in range(n-1)]}")
        for i in range(n):
            dotx = sum(B[jj][i]*x[jj] for jj in range(n))
            if abs(dotx) < mpf('1e-40'):
                print(f"  relation col {i}: {[int(mp.nint(B[jj][i])) for jj in range(n)]} dot={mp.nstr(dotx,3)}")

q0 = mpf(1)/2 + 2*mp.log(2)
pslq_debug2([q0, mpf(1), mp.log(2)], maxsteps=12)
