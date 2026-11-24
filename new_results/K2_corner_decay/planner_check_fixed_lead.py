import mpmath as mp, sys, re
mp.mp.dps=50
angles={'pi_3':mp.pi/3,'pi_2':mp.pi/2,'2pi_3':2*mp.pi/3,'3pi_4':3*mp.pi/4,'5pi_6':5*mp.pi/6}
def load(tag):
    M={}
    for line in open(f'out/mk_{tag}.tsv'):
        if line.startswith('#') or not line.strip(): continue
        k,v=line.split('\t'); t=v.strip().strip('()').replace(' ','')[:-1]; i=max(j for j in range(1,len(t)) if t[j] in '+-' and t[j-1] not in 'eE'); M[int(k)]=mp.mpf(t[:i])
    return M
for tag,a in angles.items():
    M=load(tag)
    for fixed in (False, True):
        rows=[];rhs=[]
        for k in range(150,301):
            N=k+1; th=mp.sqrt(2*a*N)
            y=M[k]*N*mp.e**th/mp.sqrt(2*a*mp.sin(a))
            if fixed: y=y+mp.cos(th)          # impose leading term -cos(th)
            row=[] if fixed else [mp.cos(th),mp.sin(th)]
            for r in range(1,7):
                row+=[N**(-mp.mpf(r)/2)*mp.cos(th), N**(-mp.mpf(r)/2)*mp.sin(th)]
            rows.append(row); rhs.append(y)
        A=mp.matrix(rows); b=mp.matrix(rhs); x=mp.lu_solve(A.T*A,A.T*b)
        res=A*x-b; rms=mp.sqrt(sum(v**2 for v in res)/len(res))
        lead = '' if fixed else f' lead=({mp.nstr(x[0],6)},{mp.nstr(x[1],6)})'
        print(f'{tag:6s} {"A=1,phi=pi fixed" if fixed else "free":17s} rms={mp.nstr(rms,3)}{lead}')
