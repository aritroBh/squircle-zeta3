import mpmath, json
from mpmath import mp, mpf, pi, nstr
C = json.load(open('ckpt_constants.json'))
paper = json.load(open('ckpt_paper_Mk.json'))
def gv(k): return mpf(C[k])
mp.dps = 245
sq2v, lamv, Lv, Li2bv, piv = gv('sq2'), gv('lam'), gv('L'), gv('Li2_b'), pi
z3v = gv('zeta3')
d3 = +mpf(paper['d3'])
grp = [('zeta3',z3v),('sq2*zeta3',sq2v*z3v),('pi^2*lam',piv**2*lamv),('pi^2*L',piv**2*Lv),
 ('lam^3',lamv**3),('lam^2*L',lamv**2*Lv),('lam*L^2',lamv*Lv**2),('L^3',Lv**3),
 ('lam*Li2(b)',lamv*Li2bv),('L*Li2(b)',Lv*Li2bv),('sq2*pi^2*lam',sq2v*piv**2*lamv),
 ('sq2*pi^2*L',sq2v*piv**2*Lv),('sq2*lam^3',sq2v*lamv**3),('sq2*lam^2*L',sq2v*lamv**2*Lv),
 ('sq2*lam*L^2',sq2v*lamv*Lv**2),('sq2*L^3',sq2v*Lv**3),('sq2*lam*Li2(b)',sq2v*lamv*Li2bv),
 ('sq2*L*Li2(b)',sq2v*Lv*Li2bv),('pi^2',piv**2),('lam^2',lamv**2),('L^2',Lv**2),('lam*L',lamv*Lv),
 ('Li2(b)',Li2bv),('sq2*pi^2',sq2v*piv**2),('sq2*lam^2',sq2v*lamv**2),('sq2*L^2',sq2v*Lv**2),
 ('sq2*lam*L',sq2v*lamv*Lv),('sq2*Li2(b)',sq2v*Li2bv),('lam',lamv),('L',Lv),('1',mpf(1)),
 ('sq2*lam',sq2v*lamv),('sq2*L',sq2v*Lv),('sq2*1',sq2v),('d3',d3)]
cand230 = {'zeta3':411493,'sq2*zeta3':207550,'pi^2*lam':211903,'pi^2*L':-191244,'lam^3':-326364,
 'lam^2*L':-179897,'lam*L^2':-482537,'L^3':-159755,'lam*Li2(b)':-279643,'L*Li2(b)':-310614,
 'sq2*pi^2*lam':-104004,'sq2*pi^2*L':-331892,'sq2*lam^3':-83395,'sq2*lam^2*L':104384,
 'sq2*lam*L^2':-144252,'sq2*L^3':142066,'sq2*lam*Li2(b)':-153017,'sq2*L*Li2(b)':135276,
 'pi^2':294585,'lam^2':-306796,'L^2':6153,'lam*L':182264,'Li2(b)':101962,'sq2*pi^2':318680,
 'sq2*lam^2':-164443,'sq2*L^2':110650,'sq2*lam*L':-94743,'sq2*Li2(b)':464285,'lam':404365,
 'L':-381088,'1':-279728,'sq2*lam':-211111,'sq2*L':-149603,'sq2*1':-53896,'d3':188907}
res = sum(mpf(a)*v for nm, v in grp for a in [cand230.get(nm,0)])
print("C=1e230 candidate residual @245dps:", nstr(res, 6))
cand180 = {'zeta3':-22240,'sq2*zeta3':29844,'pi^2*lam':-41041,'pi^2*L':-6306,'lam^3':-64729,
 'lam^2*L':44509,'lam*L^2':42638,'L^3':20256,'lam*Li2(b)':-35776,'L*Li2(b)':-158,
 'sq2*pi^2*lam':60636,'sq2*pi^2*L':10270,'sq2*lam^3':-39008,'sq2*lam^2*L':-65405,
 'sq2*lam*L^2':-52062,'sq2*L^3':-54820,'sq2*lam*Li2(b)':-18903,'sq2*L*Li2(b)':-31836,
 'pi^2':701,'lam^2':20261,'L^2':19155,'lam*L':-7766,'Li2(b)':-14806,'sq2*pi^2':24571,
 'sq2*lam^2':20925,'sq2*L^2':-4235,'sq2*lam*L':50906,'sq2*Li2(b)':-23245,'lam':34785,
 'L':1719,'1':-21888,'sq2*lam':72989,'sq2*L':42956,'sq2*1':-63638,'d3':79366}
res2 = sum(mpf(a)*v for nm, v in grp for a in [cand180.get(nm,0)])
print("C=1e180 candidate residual @245dps:", nstr(res2, 6))
with open('progress.log','a') as f:
    f.write("step_52: no-trilog candidates confirmed spurious (residuals ~1e-176, ~1e-225, not exact): d3 needs trilogarithms (NUMERICAL)\n")
