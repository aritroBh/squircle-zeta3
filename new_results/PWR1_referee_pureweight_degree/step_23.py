import pwr1_engine as E, sympy as sp, json

tot2, pure2, B2 = E.total_log(2)
P = E.ser_exp(pure2).scale(4)   # P/sqrt2 = 4 exp(log(P/(4sqrt2)))
print("P/sqrt2 through tau^8 (formal Z):")
for k in range(0, 9):
    print(f"  tau^{k}: {P.coeff(k)}")

# Compare with D16's tau_series_2d.json
js = json.load(open('../D16_pure_weight/tau_series_2d.json'))
print("\njson keys:", list(js.keys())[:20])
