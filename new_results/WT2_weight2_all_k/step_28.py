src = open('fastred.py').read()
old = """    def integ(coeffs):
        out = []
        for j, e in enumerate(coeffs, start=1):
            if j == 1:
                assert e.is_zero(), "log term in higher antiderivative"
                out.append(ZERO)
            else:
                out.append(e.scale(Fr(-1, j-1)))
        while out and out[-1].is_zero(): out.pop()
        return out"""
new = """    def integ(coeffs):
        out = []
        for j, e in enumerate(coeffs, start=1):
            if j == 1:
                assert e.is_zero(), "log term in higher antiderivative"
                continue
            out.append(e.scale(Fr(-1, j-1)))
        while out and out[-1].is_zero(): out.pop()
        return out"""
assert old in src
open('fastred.py','w').write(src.replace(old, new))
print("patched")
