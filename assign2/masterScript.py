#!/usr/bin/env python
import subprocess

ITERATIONS = 5
nodes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
delta = [1, 2, 3, 4, 5, 6, 7, 8, 10]

for N in nodes:
  for d in delta:
    subprocess.call(['benchmark.py', '-n', str(N), '-d', str(d), '-i', str(ITERATIONS)])
