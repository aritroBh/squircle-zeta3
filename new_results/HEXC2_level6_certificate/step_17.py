import json, re
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle'
def load_fit(name):
    D = json.load(open(G1+'/'+name))
    out=[]
    for key,val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if key=="1": e,w = 0,()
        elif key=="c": e,w = 1,()
        else:
            e=int(m.group(1)); w=tuple(int(z) for z in m.group(2).split(',')) if m.group(2).strip() else ()
        out.append((e,w,val))
    return out
PSI = load_fit('hpl_fit_k3_(1,0,0).json')
for e,w,v in PSI: print('   c^%d H%s  %s'%(e,list(w),v))
base = '<project>/research/08_theory'
print()
print(open(base+'/HEXC1_hexagon_certificate/word_mpl_conv.json').read()[:600])
print()
print(open(base+'/HEXC1_hexagon_certificate/w4_F4_hpl_coords.json').read())
print(open(base+'/HEXC1_hexagon_certificate/w4_F4_closed_coords.json').read())
