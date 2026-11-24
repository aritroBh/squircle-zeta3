import mpmath as mp, json
from mpmath import mpf
mp.mp.dps = 70
WORK='research/08_theory/K4d_identify'
def sermul(a,b,K):
    o=[mpf(0)]*(K+1)
    for i in range(K+1):
        for j in range(i+1): o[i]+=a[j]*b[i-j]
    return o
def serinv(a,K):
    o=[mpf(0)]*(K+1); o[0]=1/a[0]
    for n in range(1,K+1): o[n]=-sum(a[j]*o[n-j] for j in range(1,n+1))/a[0]
    return o
L = (mpf(-6)/5)*mp.log(2) + (mpf(13)/5)*mp.log(3)
u=[mpf(0),1,-1,1,-1]
Lu=[L*x for x in u]; e=[mpf(0)]*5; e[0]=1; pw=[mpf(0)]*5; pw[0]=1
for k in range(1,5):
    pw=sermul(pw,Lu,4)
    for i in range(4,-1,-1): e[i]+=pw[i]/mp.factorial(k)
finv=serinv(e,4)
def Tseries(SA):
    return sermul(list(SA),finv,4)
T64 = Tseries([mpf(s) for s in json.load(open(WORK+'/modelA_h64.json'))['Sdk']])
T96 = Tseries([mpf(s) for s in json.load(open(WORK+'/modelA_h96.json'))['Sdk']])
S0 = T96[0]
for k in [2,3]:
    a=T64[k]/S0; b=T96[k]/S0
    print("T%d/S0  h64=%s"%(k,mp.nstr(a,50)))
    print("        h96=%s"%mp.nstr(b,50))
    print("        |h96-h64|=%s  => ~%d digits agree\n"%(mp.nstr(abs(a-b),2), -mp.log10(abs(a-b))))
