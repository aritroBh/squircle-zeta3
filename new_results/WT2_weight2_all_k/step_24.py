src = open('fastred.py').read()
old = """            Drem2 = pdeflate(Drem, r)
            if peval(Drem2, r).is_zero():
                raise ValueError(f"higher real pole at {r}")
            res = peval(N, r) / peval(Drem2, r)"""
new = """            Drem2 = pdeflate(Drem, r)
            if peval(Drem2, r).is_zero():
                raise ValueError(f"higher real pole at {r}")
            res = peval(N, r) / peval(pdeflate(list(D), r), r)"""
assert old in src
src = src.replace(old, new)
open('fastred.py','w').write(src)
print("patched")
