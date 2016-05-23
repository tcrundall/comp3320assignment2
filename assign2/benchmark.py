#!/usr/bin/env python
import argparse
import subprocess
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--iter',  dest='iter',  help='number of trials')
parser.add_argument('-n', '--nodes', dest='nodes', help='nodes per dimension')
parser.add_argument('-d', '--delta', dest='delta', help='range of bonds')

args = parser.parse_args()

seperation = 10./int(args.nodes)

eval_per_call = []
loop_per_call = []
cumulative    = []

for i in range(int(args.iter)):

  subprocess.call(["kernel_main", '-n', args.nodes, '-d', args.delta, '-s', str(seperation)])

  with open("prof_output", 'w') as prof_output:
    subprocess.call(["gprof", "kernel_main", "gmon.out"], stdout=prof_output)


  with open("prof_output", 'r') as results:
    for j in range(4):
      results.readline()

    line = results.readline().strip().split()
    units = line[4]

    line = results.readline().strip().split()
    eval_per_call.append(float(line[4]))
  
    line = results.readline().strip().split()
    loop_per_call.append(float(line[4]))
  
    line = results.readline().strip().split()
    cumulative.append(float(line[1]))

savefile = "results/N" + args.nodes + "d" + args.delta

units = units[:2]

with open(savefile, 'w') as f:
  eval_time = sum(eval_per_call)/len(eval_per_call)
  loop_time = sum(loop_per_call)/len(loop_per_call)
  total_time = sum(cumulative)/len(cumulative)

  f.write(args.nodes + " nodes and " + args.delta + " delta\n")
  f.write(str(eval_time) + ' ' + units + '\n')
  f.write(str(loop_time) + ' ' + units + '\n')
  f.write(str(total_time) + ' s\n')

