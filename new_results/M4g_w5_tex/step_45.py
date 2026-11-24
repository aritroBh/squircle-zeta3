import mpmath as mp
import json
cert=json.load(open('../M4e_w5_lchi/Lchi8_5.json'))['value']
for dps in (450, 650):
    mp.dps=dps
    h=[mp.zeta(5,mp.mpf(k)/8) for k in (1,3,5,7)]
    Lchi=(h[0]+h[3]-h[1]-h[2])/mp.mpf(8)**5
    print('dps',dps,'LCHI-cert =',mp.nstr(Lchi-mp.mpf(cert),3))
    print('dps',dps,'zeta3 =',mp.nstr(mp.zeta(3),30))
    print('dps',dps,'zeta5 =',mp.nstr(mp.zeta(5),30))
    print('dps',dps,'log2  =',mp.nstr(mp.log(2),30))
