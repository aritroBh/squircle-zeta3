import mpmath as mp
mp.mp.dps=60
# G(r,s;1) = int_0^1 log(1-t/s)/(t-r) dt. Candidate: Li2((1-r)/(s-r)) - Li2(-r/(s-r)) + logs?
def G2num(r,s):
    return mp.quad(lambda t: mp.log(1-t/s)/(t-r), [0,1])
r,s = mp.mpf('2.4'), mp.mpf('-0.4')
v = G2num(r,s)
c1 = Li21 = mp.polylog(2,(1-r)/(s-r))
c2 = mp.polylog(2,(-r)/(s-r))
print("G:", mp.nstr(v,18))
print("Li2((1-r)/(s-r))-Li2(-r/(s-r)):", mp.nstr(c1-c2,18))
# test the shuffle-derived formula: G(r,s;1)+G(s,r;1)=log(1-1/r)log(1-1/w)...
v2 = G2num(s,r)
print("shuffle check:", mp.nstr(v+v2,18), mp.nstr(mp.log(1-1/r)*mp.log(1-1/s),18))
# G(r,s;1) = -Li2(s/(s-r)) + ... ? try integral formula: int log(1-t/s)/(t-r) 
# substitute u=(t-r)/(s-r): t = r+u(s-r), 1-t/s = (s-r)(1-u)/(s)... 
# t-r = u(s-r), dt = (s-r)du: G = int_0^{(1-r)/(s-r)} log((s-r)(1-u)/s)/u du
#   = int log((s-r)/s)/u + log(1-u)/u du  (limits 0 -> (1-r)/(s-r))
#   = log((s-r)/s)*log((1-r)/(s-r)) - Li2((1-r)/(s-r)) + ... wait int log(1-u)/u du = -Li2(u)
print("candidate:", mp.nstr(mp.log((s-r)/s)*mp.log((1-r)/(s-r)) - mp.polylog(2,(1-r)/(s-r)),18))
