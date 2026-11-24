import os
for d in ["../OCT6_octagon_closed", "../OCT7_generic_octagon"]:
    print("="*30, d)
    try:
        print(os.listdir(d))
    except Exception as e:
        print("ERR", e)
