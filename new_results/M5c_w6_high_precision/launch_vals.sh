#!/bin/bash
cd "<project>/research/08_theory/M5c_w6_high_precision"
nohup nice -n 10 "<project>/.venv/bin/python" -B w6_values.py 3500 > w6_vals_3500.out 2>&1 &
echo $! > w6_vals_3500.pid
