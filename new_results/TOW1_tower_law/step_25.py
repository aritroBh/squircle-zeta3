from towlib import *
for p in (2,3,4):
    maxd = min(4*p-4, 2*(2*p-1)-4)  # parent max depth / child max depth
    # also limited by data: P needs levels 2p+1+floor(d'/2) <= 14
    maxd = min(maxd, 2*(14-2*p-1))
    print(f'--- p={p} ---')
    Deltas = {}
    for dp in range(0, maxd+1):
        D = rho((p,p),dp) - polyval_diag(p,dp,2*p)
        Deltas[dp] = D
        kind = 'real' if D.is_real() else ('pure-i' if (D.p0==0 and D.p1==0) else 'COMPLEX')
        print(f'  Delta({dp}) = {D}  [{kind}]')
    for delta in range(0, maxd+1):
        feed = feed_diag_into_child(p, 2*p, delta, lambda dp: Deltas.get(dp, ZERO))
        print(f'  feed into child at delta={delta}: zero={feed.is_zero()}  {feed if not feed.is_zero() else ""}')
