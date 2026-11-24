import json
d = json.load(open('ckpt_pslq_d3.json'))
got = d['d3_rel']
expected = {'zeta3':-228,'Lchi8_3':-7680,'Li3(b^2)':-528,'sq2*zeta3':-567,'sq2*Li3(1/sq2)':5184,
 'pi^2*lam':8,'pi^2*L':112,'lam^3':-32,'lam*L^2':192,'L^3':-800,'L*Li2(b)':-2304,
 'sq2*pi^2*lam':260,'sq2*lam^3':-108,'sq2*lam^2*L':-264,'sq2*lam*L^2':336,'sq2*L^3':-32,
 'sq2*lam*Li2(b)':672,'pi^2':-24,'lam^2':144,'L^2':-384,'sq2*pi^2':152,'sq2*L^2':-960,
 'sq2*lam*L':768,'sq2*Li2(b)':-1920,'lam':-192,'sq2*L':192,'d3':24}
for k in sorted(set(got)|set(expected)):
    g0, e0 = got.get(k,0), expected.get(k,0)
    if g0 != e0: print("DIFF", k, "got", g0, "expected", e0)
print("done")
