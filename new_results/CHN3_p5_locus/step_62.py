import sys; sys.path.insert(0,'.')
import point_engine as PE
import cProfile, pstats, io
pr = cProfile.Profile()
pr.enable()
r = PE.reduced_at((1,2,3,4,5), 10)   # W=10 for speed
pr.disable()
s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
ps.print_stats(12)
print(s.getvalue()[:2500])
