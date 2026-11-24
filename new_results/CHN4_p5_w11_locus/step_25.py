import random
def rank_mod(rows, p):
    M=[[x%p for x in r] for r in rows]; nr=len(M); nc=len(M[0]); r=0
    for c in range(nc):
        piv=next((i for i in range(r,nr) if M[i][c]%p), None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        inv=pow(M[r][c]%p,p-2,p); M[r]=[(v*inv)%p for v in M[r]]
        for i in range(nr):
            if i!=r and M[i][c]%p:
                f=M[i][c]%p; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(nc)]
        r+=1
        if r==nr: break
    return r
p=(1<<61)-1
# identity 170
I=[[1 if i==j else 0 for j in range(170)] for i in range(170)]
print("rank identity 170:", rank_mod(I,p))
# random integer matrix 170x170
random.seed(5)
R=[[random.randint(0,10**6) for _ in range(170)] for _ in range(170)]
print("rank random 170x170:", rank_mod(R,p))
# rank-deficient: duplicate a row
R2=[row[:] for row in R]; R2[5]=R2[3][:]
print("rank with duplicated row (should be 169):", rank_mod(R2,p))
