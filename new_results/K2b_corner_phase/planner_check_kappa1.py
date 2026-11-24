import mpmath as mp
mp.mp.dps=40
angles={'pi_3':mp.pi/3,'pi_2':mp.pi/2,'2pi_3':2*mp.pi/3,'3pi_4':3*mp.pi/4,'5pi_6':5*mp.pi/6}
def load(tag):
    M={}
    for line in open(f'../K2_corner_decay/out/mk_{tag}.tsv'):
        if line.startswith('#') or not line.strip(): continue
        k,v=line.split('\t'); t=v.strip().strip('()').replace(' ','')[:-1]
        i=max(j for j in range(1,len(t)) if t[j] in '+-' and t[j-1] not in 'eE'); M[int(k)]=mp.mpf(t[:i])
    return M
def kappa1(a):
    C=1/(2*mp.sqrt(a))+a**1.5/12+3*mp.sqrt(a)/8*mp.cot(a)
    return (C+mp.sqrt(a)/2)/mp.sqrt(2)+1j*(C-mp.sqrt(a)/2)/mp.sqrt(2)
d13b=((24+12*mp.pi+mp.pi**2)+1j*(24-12*mp.pi+mp.pi**2))/(48*mp.sqrt(mp.pi))
print('kappa1(pi/2) formula', mp.nstr(kappa1(mp.pi/2),15), ' D13b', mp.nstr(d13b,15), ' conj(D13b)', mp.nstr(mp.conj(d13b),15))
for tag,a in angles.items():
    M=load(tag); A=mp.sqrt(2*a*mp.sin(a))
    out=[]
    for lab,k1 in [('lead only',0),('kappa1',kappa1(a)),('conj kappa1',mp.conj(kappa1(a)))]:
        mx=0
        for k in range(150,301):
            N=k+1; th=mp.sqrt(2*a*N)
            pred=-(A/N)*mp.e**(-th)*mp.re(mp.e**(1j*th)*(1+k1/mp.sqrt(N)))
            rho=(M[k]-pred)*N*mp.e**th/A
            mx=max(mx,abs(rho)*(N if k1!=0 else mp.sqrt(N)))
        out.append(f'{lab}: max|rho|*{"N" if lab!="lead only" else "sqrtN"} = {mp.nstr(mx,4)}')
    print(f'{tag:6s} '+' | '.join(out))
