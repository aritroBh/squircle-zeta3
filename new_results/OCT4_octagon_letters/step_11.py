import os, json, mpmath as mp
work = 'research/08_theory/OCT4_octagon_letters'
os.chdir(work)
base = 'research/08_theory/'

# init progress log
with open('progress.log','w') as f:
    f.write("OCT4 progress log\n")

def log(msg):
    with open('progress.log','a') as f:
        f.write(msg.rstrip()+"\n")
    print(msg)

# Load data
hidata = json.load(open(base+'OCT3_octagon_basis/hidata.json'))
basis34 = json.load(open(base+'OCT3_octagon_basis/basis34.json'))
consts62 = json.load(open(base+'OCT3_octagon_basis/consts62.json'))
vals = json.load(open(base+'D15b_polygon_law/vals_m2_8_200.json'))

print("hidata keys:", list(hidata.keys()))
print("basis34 keys:", list(basis34.keys())[:10], "..." if len(basis34)>10 else "")
print("consts62 keys:", list(consts62.keys())[:10])
print("vals keys:", list(vals.keys()))
print("vals['4'] entries:", len(vals['4']))
for e in vals['4']:
    print("   ", e[:60])
log("step_03: loaded hidata, basis34, consts62, vals. vals['4'] has %d entries" % len(vals['4']))
