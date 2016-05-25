import matplotlib.pyplot as plt
import math

codes = ["ORIG", "SEQ", "SSE"]

nodes = [10, 20, 30, 40, 50, 60]
deltas = [1, 2, 3, 4, 5, 6, 7, 8, 10]

data =[ [ [0.0 for delta in deltas] for node in nodes] for i in range(3)]

units =[ [ ["" for delta in deltas] for node in nodes] for i in range(3)]

node = 10

for k, code in enumerate(codes):
  for i, node in enumerate(nodes):
    for j, delta in enumerate(deltas):
      with open(str(code)+"_RESULTS/N"+str(node)+"d"+str(delta),'r') as f:
        f.readline()
        line = f.readline().strip().split()
        units[k][i][j] = line[1]
        eval_per_call = float(line[0])
        loop_per_call = float(f.readline().strip().split()[0])
        total         = float(f.readline().strip().split()[0])
        if units[k][i][j] == "us":
          eval_per_call = eval_per_call/1000
          loop_per_call = loop_per_call/1000

        data[k][i][j] = (eval_per_call, loop_per_call, total)

#plotting eval_per_call over nodes

xlist = nodes
delta = 3
d_index = deltas.index(delta)
def main():
    for j, code in enumerate(codes):
        ylist = [math.sqrt(data[j][i][d_index][2]) for i, node in enumerate(nodes)]

        plt.xlabel("Nodes per dimension")
        plt.ylabel("sqrt(time) (s)")
        plt.plot(xlist,ylist,'--ro')
        string = code + " sqrt total time vs Nodes, delta=" + str(delta)
        plt.title(string)
        plt.savefig(string)
        plt.clf()

def main2():
    for j, code in enumerate(codes):
        node_index = 4
        xlist = deltas
        ylist = [(data[j][node_index][i][0]) for i, delta in enumerate(deltas)]

        plt.xlabel("delta")
        plt.ylabel("time (ms)")
        plt.plot(xlist,ylist,'--ro')
        string = code + " eval per call vs deltas, nodes=" + str(nodes[node_index])
        plt.title(string)
        plt.savefig(string)
        plt.clf()

def main4():
    node_index = 0
    xlist = deltas
    ylist = [math.sqrt(data[0][node_index][i][0]) for i, delta in enumerate(deltas)]
    zlist = [math.sqrt(data[1][node_index][i][0]) for i, delta in enumerate(deltas)]
    qlist = [math.sqrt(data[2][node_index][i][0]) for i, delta in enumerate(deltas)]

    plt.xlabel("delta")
    plt.ylabel("time (ms)")
    plt.plot(xlist,ylist,'--ro', label='original')
    plt.plot(xlist,zlist,'--bo', label='improved')
    plt.plot(xlist,qlist,'--go', label='SSE')
    plt.legend(loc=2)
    string = "sqrt eval per call vs deltas, nodes=" + str(nodes[node_index])
    plt.title(string)
    plt.savefig(string)
    plt.clf()

def main3():
    ylist = [math.sqrt(data[0][i][d_index][0]) for i, node in enumerate(nodes)]
    zlist = [math.sqrt(data[1][i][d_index][0]) for i, node in enumerate(nodes)]
    qlist = [math.sqrt(data[2][i][d_index][0]) for i, node in enumerate(nodes)]
    plt.xlabel("Nodes per dimension")
    plt.ylabel("sqrt(time) (ms^1/2)")
    plt.plot(xlist,ylist,'--ro', label='original')
    plt.plot(xlist,zlist,'--bo', label='improved')
    plt.plot(xlist,qlist,'--go', label='SSE')
    plt.legend(loc=2)
    string = "Improvement of sqrt eval per call vs Nodes, delta=" + str(delta)
    plt.title(string)
    plt.savefig(string)
    plt.clf()
main4()
