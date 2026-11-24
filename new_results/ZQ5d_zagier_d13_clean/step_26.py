import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
for f in ['lll_capd17_P1500','lll_capd17_P1800','lll_capfull_P1500','lll_capfull_P1800']:
    J = json.load(open(c+f'/out/d13/{f}.json'))
    print(f, 'P=',J['P'],'dt=',round(J['dt'],1),'nrel=',J['nrel'], 'pool=', J['pool'])
# also d29 timings
for f in ['lll_S7h2_P1500','lll_S7h2full_P1500','lll_timing']:
    try:
        J = json.load(open(c+f'/out/d29/{f}.json'))
        print('d29', f, {k:J[k] for k in J if k in ('P','dt','nrel','dim')})
    except Exception as e:
        print('d29', f, 'err', e)
