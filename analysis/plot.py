#!/bin/python3

import matplotlib.pyplot as plt
import load

def runningAverage(data):
    Av  = len(data) * [0]
    Sum = 0

    for idx in range(0, len(data)):
        Sum += data[idx]
        Av[idx] = Sum / (idx + 1)

    return Av

t1 = load.Col("energy_single_0.002.log", 1); t1 = [0.002 * x for x in t1]
E1 = load.Col("energy_single_0.002.log", 8)

t2 = load.Col("energy_double_0.002.log", 1); t2 = [0.002 * x for x in t2]
E2 = load.Col("energy_double_0.002.log", 8)

t3 = load.Col("energy_single_0.00002.log", 1); t3 = [0.00002 * x for x in t3]
E3 = load.Col("energy_single_0.00002.log", 8)

t4 = load.Col("energy_single_0.0000002.log", 1); t4 = [0.0000002 * x for x in t4]
E4 = load.Col("energy_single_0.0000002.log", 8)

t5 = load.Col("energy_double_0.0000002.log", 1); t5 = [0.0000002 * x for x in t5]
E5 = load.Col("energy_double_0.0000002.log", 8)

plt.figure(figsize=(12, 5))

# plt.plot(t1, E1, color='red', label="single, dt = 2 fs")
# plt.plot(t1, runningAverage(E1), color='red')

# plt.plot(t3, E3, color='green', label="single, dt = 0.02 fs")

# plt.plot(t2, E2, color='blue', label="double, dt = 2 fs")
# plt.plot(t2, runningAverage(E2), color='blue')

plt.plot(t4, E4, label="single, dt = 0.0002 fs", color='red')

plt.plot(t5, E5, label="double, dt = 0.0002 fs", color='blue')

plt.xlabel("Time (ps)")
plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 3))
plt.ylabel("Conserved Energy (kJ/mol)")
plt.legend()
plt.tight_layout()
plt.savefig("plot2.pdf")
plt.show()
