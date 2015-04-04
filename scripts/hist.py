#!/usr/bin/env python

import sys
import math

import numpy as np
import matplotlib.pyplot as plt

plus = []
minus = []

for line in sys.stdin:
    line = line.strip()
    if not line or line.startswith("ORF"):
        continue
    tokens = line.split()
    if len(tokens) < 3:
        continue

    status, e_val = tokens[2], float(tokens[4])
    if status == "+":
        plus.append(e_val)
    else:
        minus.append(e_val)

plus = map(lambda p: -math.log(p), plus)
minus = map(lambda p: -math.log(p), minus)

bins = np.linspace(0, 140, 30)

plt.hist(plus, bins, alpha=0.5, label="matched")
plt.hist(minus, bins, alpha=0.5, label="not matched")
plt.legend(loc="upper right")
plt.show()
