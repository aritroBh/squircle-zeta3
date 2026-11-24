
import sys, time, pickle
sys.dont_write_bytecode = True
import rul1, flint

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

data = pickle.load(open('rels8_dsdual.pkl','rb'))
cols8 = sorted(rul1.conv_words_of_weight(8))
colpos = {c:i for i,c in enumerate(cols8)}
P = 2**62 - 57   # prime? check
import sympy as sp
P = int(sp.prevprime(2**62))
print("prime", P)

t = time.time()
A = flint.nmod_mat(len(data), len(cols8), P)
for i, (fam, tag, r) in enumerate(data):
    for k, (n, d) in r.items():
        A[i, colpos[k]] = (n * pow(d, -1, P)) % P
print("build matrix %.0fs"%(time.time()-t), flush=True)
t = time.time()
r = A.rank()
print("flint rank: %d in %.0fs" % (r, time.time()-t), flush=True)
log("flint nmod_mat rank w8: %d (codim %d, DG 34)" % (r, len(cols8)-r))
