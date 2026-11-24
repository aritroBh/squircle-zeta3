import os, json, mpmath as mp
WORK = 'research/08_theory/OCT4_octagon_letters'
BASE = 'research/08_theory/'
def load():
    hidata = json.load(open(BASE+'OCT3_octagon_basis/hidata.json'))
    basis34 = json.load(open(BASE+'OCT3_octagon_basis/basis34.json'))
    consts62 = json.load(open(BASE+'OCT3_octagon_basis/consts62.json'))
    vals = json.load(open(BASE+'D15b_polygon_law/vals_m2_8_200.json'))
    return hidata, basis34, consts62, vals
def log(msg):
    with open(WORK+'/progress.log','a') as f:
        f.write(msg.rstrip()+"\n")
