import matplotlib.pyplot as plt

nodes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
deltas = [1, 2, 3, 4, 5, 6, 7, 8, 10]

data = [ [0.0 for delta in deltas] for node in nodes]

units = [ ["" for delta in deltas] for node in nodes]

node = 10
for i, node in enumerate(nodes):
  for j, delta in enumerate(deltas):
    with open("SSE_RESULTS/N"+str(node)+"d"+str(delta),'r') as f:
      f.readline()
      line = f.readline().strip().split()
      units[i][j] = line[1]
      eval_per_call = float(line[0])
      loop_per_call = float(f.readline().strip().split()[0])
      total         = float(f.readline().strip().split()[0])
      if units[i][j] == "us":
        eval_per_call = eval_per_call/1000
        loop_per_call = loop_per_call/1000

      data[i][j] = (eval_per_call, loop_per_call, total)

#plotting eval_per_call over nodes
xlist = nodes
delta = 2

ylist = [data[i][delta][0] for i, node in enumerate(nodes)]

print ylist
