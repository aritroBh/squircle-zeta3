#!/bin/bash
cd "<project>/research/08_theory/M5c_w6_high_precision"
nohup nice -n 10 "<project>/.venv/bin/python" -B w6c_lll.py 3300 w6x14_15_9 > w6c_lll_x14_15_9_3300.out 2>&1 &
echo $! > w6c_lll_x14_15_9_3300.pid
