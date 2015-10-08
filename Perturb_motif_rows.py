import random
import numpy as np
import sys

f = open('perturbations.txt', 'w')
with open("pfm_all.txt") as handle:
 for m in motifs.parse(handle, "jaspar"):
    counts = m.counts
    values = list()
    ncol = len(counts[1,:])
    for x in range(0,ncol): 
      for y in range(0,4):
        values.append(counts[y,x])
    new_counts = np.reshape(np.matrix(values), (4,ncol), order="F")
    for x in range(0,ncol):
      for y in range(0,20):
        a = random.randint(0,3)
        b = random.randint(0,3)
        old_a = new_counts[a,x]
        old_b = new_counts[b,x]
        new_counts[a,x] = old_b
        new_counts[b,x] = old_a
    f.write(">%s %s\n"%(m.matrix_id,m.name))
    for x in range(0,4):
      for y in range(0, ncol):
        f.write(str(int(new_counts[x,y])))
        f.write("\t")
      f.write("\n")

f.close() 
