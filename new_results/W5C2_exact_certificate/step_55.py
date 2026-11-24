import os, signal, subprocess, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
pid = int(open('step_01.pid').read())
try:
    os.killpg(pid, signal.SIGKILL)  # kill the whole session group of the bash driver
    print('killed pgid', pid)
except Exception as e:
    print('kill err', e)
    try:
        os.kill(pid, signal.SIGKILL); print('killed pid')
    except Exception as e2: print(e2)
time.sleep(2)
# add DenseEchelon to the library
lib = open('w5c2_lib.py').read()
dense = '''

class DenseEchelon:
    """Row echelon over GF(p) with DENSE pivot rows (fast at high density).
    Invariant: pivot row created at time t has zeros in all pivot columns of time < t
    (rem contains no pivot columns at creation). Reduction sweeps pivot columns in
    increasing creation time; eliminating pivot t only introduces nonzeros in pivot
    columns of time > t, so a rescan after each elimination is correct."""
    def __init__(self, ncols, p, maxpiv=None):
        self.NC = ncols; self.p = p
        self.maxpiv = maxpiv or ncols
        self.piv = np.zeros((self.maxpiv, ncols), np.int64)
        self.pivcol = np.full(self.maxpiv, -1, np.int64)
        self.pivrow = []
        self.r = 0
        self.elims = 0
    def reduce(self, cols, vals):
        p = self.p; NC = self.NC
        buf = np.zeros(NC, np.int64)
        buf[cols] = np.asarray(vals, np.int64) % p
        r = self.r
        if r:
            pc = self.pivcol[:r]
            elims = 0
            while True:
                nz = pc[buf[pc] != 0]
                if nz.size == 0: break
                c = int(nz[0])           # smallest creation-time pivot column present
                # find its pivot index: pivot cols are unique, creation order = index order
                t = int(np.searchsorted(pc, c))
                a = int(buf[c])
                buf = (buf - a*self.piv[t]) % p
                elims += 1
            self.elims += elims
        return buf
    def add_pivot(self, buf, src=-1):
        nz = np.nonzero(buf)[0]
        assert nz.size
        c = int(nz[0])
        lead = int(buf[c]) % self.p
        inv = pow(lead, -1, self.p)
        self.piv[self.r] = (buf*inv) % self.p
        self.pivcol[self.r] = c
        self.pivrow.append(src)
        self.r += 1
        return c
    def feed(self, rows, label='', t0=None, logevery=4000, stopcheck=None):
        for ri in range(len(rows)):
            cols, vals = rows[ri]
            buf = self.reduce(cols, vals)
            if buf.any(): self.add_pivot(buf, src=ri)
            if (ri+1) % logevery == 0:
                el = time.monotonic()-t0 if t0 else 0
                log(f'{label} row {ri+1}/{len(rows)} pivots={self.r} elims={self.elims/1e6:.1f}e6 el={el:.0f}s')
            if stopcheck and (ri+1) % stopcheck[0] == 0:
                if stopcheck[1](self):
                    log(f'{label} early stop at row {ri+1} pivots={self.r}')
                    return ri+1
        return None
'''
open('w5c2_lib.py','w').write(lib + dense)
print('DenseEchelon added')
print(subprocess.run(['bash','-c','ls work/'], capture_output=True, text=True).stdout)