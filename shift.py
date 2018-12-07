import numpy as np

coords = []
masses = {"O":15.9990,"H":1.0079,"M":0.0000,"D":0.0e0}

newx = []
newy = []
newz = []
ids  = []

xcom = 0.0e0
ycom = 0.0e0
zcom = 0.0e0
mtot = 0.0e0

with open("restart.xyz","r") as f:
  for i,line in enumerate(f):
    if i > 1: 
      cl = line.strip().split()
      mid = cl[0]
      mx  = float(cl[1])
      my  = float(cl[2])
      mz  = float(cl[3])
      xcom += mx*masses[mid]
      ycom += my*masses[mid]
      zcom += mz*masses[mid]
      mtot += masses[mid]
      ids.append(mid)
      newx.append(mx)
      newy.append(my)
      newz.append(mz)
xcom /= mtot
ycom /= mtot
zcom /= mtot
newx = (np.array(newx))-xcom + 7.00
newy = (np.array(newy))-ycom + 7.00
newz = (np.array(newz))-zcom + 7.00

with open("restart.xyz","w") as f:
  n = len(ids)
  f.write("{}\n\n".format(n))
  for i in range(n):
    f.write("{} {:<8.4f} {:<8.4f} {:<8.4f}\n".format(ids[i],newx[i],newy[i],newz[i]))

