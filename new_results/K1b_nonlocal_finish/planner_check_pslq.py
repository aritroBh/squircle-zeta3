import json, mpmath as mp
mp.mp.dps=150
d=json.load(open('../10_frontier/E1_edge_law/verify/out/t3_b0_150.json'))
for a in d['arcs']: print({k:(str(v)[:40]) for k,v in a.items()})
z3=mp.zeta(3); L=mp.log(2); pi=mp.pi
basis=lambda X:[X,z3,pi**2*L,L**3,pi**2,L**2,L,mp.mpf(1)]
r1=[500,-231,-111,-52,80,1440,-2400,0]
r3=[128,15,-4,-1024,4,3072,-1536,0]
for a in d['arcs']:
    P3=mp.mpf(a['P3_arc']); P0=mp.mpf(a['P0'])
    for X,name in [(P3,'P3'),(P3/P0,'P3/P0')]:
        for r,rn in [(r1,'r1'),(r3,'r3')]:
            v=mp.fsum(c*b for c,b in zip(r,basis(X)))
            if abs(v)<mp.mpf(10)**-100: print('edge',a['edge'],name,rn,'residual',mp.nstr(v,5))
tot=mp.mpf(d['P3_total']); s=mp.fsum(mp.mpf(a['P3_arc']) for a in d['arcs']); print('sum arcs - total', mp.nstr(s-tot,5), 'n arcs', len(d['arcs']))
print('target coeff', mp.nstr(mp.mpf(231)/125*mp.sqrt(5)-mp.mpf(15)/64*mp.sqrt(2),20), ' N_e(2)', mp.nstr(mp.mpf(159)/500*mp.sqrt(5),30))
