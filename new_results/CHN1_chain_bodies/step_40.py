import os, time, json
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
print(open('blind_check.out').read()[-400:] if os.path.exists('blind_check.out') else 'no out yet')
if os.path.exists('blind_result.json'):
    res=json.load(open('blind_result.json'))
    pred=open('blind_prediction.txt').read().strip()
    print('predicted:', pred[:60])
    print('computed :', res['8'][0][:60])
    print('imag     :', res['8'][1])
    from mpmath import mp, mpf, nstr
    mp.dps=90
    d=mpf(res['8'][0])-mpf(pred)
    print('|diff| =', nstr(abs(d),3))
    print('c0:', res['0'])
