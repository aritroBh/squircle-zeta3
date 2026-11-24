import json
j3 = json.load(open('../../10_frontier/H1_octahedron_end/routeA_n3_K12.json'))
print("eta list first 4:", j3['3']['eta'][:4])
print("d list first 4:", j3['3']['d'][:4])
src = open('../../10_frontier/H1_octahedron_end/01_series_routeA.py').read()
print(src[:2500])
