
import sys
import pandas as pd
import numpy as np
import math
import mysymnmf as sy


k = (int)(sys.argv[1])
goal = sys.argv[2]
txtFile = sys.argv[3]
EPS = 0.0001
iter = 300
df1 = pd.read_csv(txtFile,sep = ",",header=None,float_precision='high')
points = df1.values.tolist()
if (k >= len(points)):
    print("Invalid number of clusters!")
    exit()
if(goal == "sym"):
    result = sy.sym(points)
if(goal == "ddg"):
    result = sy.ddg(points)
if(goal == "norm"):
    result = sy.norm(points)
if (goal == "symnmf"):
    w = sy.norm(points)
    m = np.mean(w)
    n = len(points)
    np.random.seed(0)
    H = [[np.random.uniform(0,2*math.sqrt(m/k)) for _ in range(k)] for _ in range(n)]
    result = sy.symnmf(H,w,iter,EPS)


for i in range(len(result)):
    for j in range(len(result[i])):
        sys.stdout.write(f"{result[i][j]:.4f},")
    print()

