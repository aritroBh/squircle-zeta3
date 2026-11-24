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
    # model: M_k N e^{th}/sqrt(2 a sin a) = sum_r N^{-r/2}(P_r cos th + Q_r sin th), th = sqrt(2aN), r=0..R
    for c_mult,label,K0 in [(2,'c=2a',120),(2,'c=2a',180),(None,'c=pi',120)]:
        rows=[];rhs=[]
        R=4
        for k in range(K0,301):
            N=k+1; c=(2*a if c_mult else mp.pi); th=mp.sqrt(c*N)
            y=M[k]*N*mp.e**th/mp.sqrt(2*a*mp.sin(a))
            row=[]
            for r in range(R+1):
                row+= [N**(-mp.mpf(r)/2)*mp.cos(th), N**(-mp.mpf(r)/2)*mp.sin(th)]
            rows.append(row); rhs.append(y)
        A=mp.matrix(rows); b=mp.matrix(rhs)
        x=mp.lu_solve(A.T*A, A.T*b)
        res=A*x-b; rms=mp.sqrt(sum(v**2 for v in res)/len(res))
        amp=mp.sqrt(x[0]**2+x[1]**2); phi=mp.atan2(-x[1],x[0])
        print(f'{tag:6s} {label:5s} lead amp={mp.nstr(amp,8)} phase={mp.nstr(phi,8)} (phase/pi={mp.nstr(phi/mp.pi,6)}, alpha/2={mp.nstr(a/2,6)}) rms={mp.nstr(rms,3)}')
