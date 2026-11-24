
"""Numeric MZV values: exact partial sums with Richardson extrapolation in 1/N.
Only for words with all indices >=2 (no log terms)."""
from mpmath import mp, mpf
def mzv_partial(word, N):
    prev=None
    for i in range(len(word)-1,-1,-1):
        s=word[i]; T=[mpf(0)]*(N+1)
        if prev is None:
            for m in range(1,N+1): T[m]=T[m-1]+mpf(1)/m**s
        else:
            for m in range(1,N+1): T[m]=T[m-1]+prev[m-1]/m**s
        prev=T
    return prev[N]
def neville(xs, ys, x0=mpf(0)):
    n=len(xs); Q=list(ys)
    for j in range(1,n):
        for i in range(n-1,j-1,-1):
            Q[i]=((x0-xs[i-j])*Q[i]-(x0-xs[i])*Q[i-1])/(xs[i]-xs[i-j])
    return Q[n-1]
def mzv_num(word, npts=16, N0=400, step=140):
    assert all(s>=2 for s in word), 'Richardson only for all indices >=2'
    xs=[]; ys=[]
    for i in range(npts):
        N=N0+i*step
        xs.append(mpf(1)/N); ys.append(mzv_partial(word,N))
    return neville(xs,ys)
