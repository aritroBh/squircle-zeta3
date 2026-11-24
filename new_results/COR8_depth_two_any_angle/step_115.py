
import pickle, time
import mpmath as mp
mp.mp.dps = 90
DS = {eval(k): mp.mpf(v) for k,v in pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/series_D_v3.pkl','rb'))['DS'].items()}
B = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/invisible_series.pkl','rb'))
print(type(B), type(B[0]) if isinstance(B, list) else B.keys())
