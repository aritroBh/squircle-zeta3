import time, json, os
T0 = time.monotonic()
def log(s): print(f"[{time.monotonic()-T0:6.1f}s] {s}", flush=True)
exec(open('step_13.py').read().split("results = {}")[0])  # reuse MyCyclo, prime_for, map_count defs

base = '<project>/research/08_theory'
stored = json.load(open(os.path.join(base,'COR13C_third_corner_certificate','map_counts.json')))
for N in (18, 24, 30, 36):
    strict, extra, exps, p = map_count(N)
    st = stored[str(N)]
    match = (strict == st['strict'] and extra == st['relaxed_extra'] and sorted(exps) == sorted(tuple(e) for e in st['exps']))
    log(f"N={N}: prime {p} (k={(p-1)//N}); exact strict={strict} (2N={2*N}) extra={extra}; matches map_counts.json: {match}")
    assert strict == 2*N and extra == 0 and match
with open('progress.log','a') as f:
    f.write("referee step(map counts hi): exact strict counts N=18,24,30,36 = 36,48,60,72 = 2N; match map_counts.json\n")
print("STEP PASS")
