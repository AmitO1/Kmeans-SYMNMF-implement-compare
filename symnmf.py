
import sys
import pandas as pd
import numpy as np
import math
import mysymnmf as sy

def print_result(result):
    for i in range(len(result)):
        for j in range(len(result[i])):
            if (j != len(result[i]) -1):
                sys.stdout.write(f"{result[i][j]:.4f},")
            else:
                sys.stdout.write(f"{result[i][j]:.4f}")
        print()

num_args = len(sys.argv) - 1
if(num_args == 3):
    k = (int)(sys.argv[1])
    goal = sys.argv[2]
    txtFile = sys.argv[3]
elif(num_args == 2):
    goal = sys.argv[1]
    txtFile = sys.argv[2] 
else:
    print("An error has occured!")
    exit()

EPS = 0.0001
iter = 300
df1 = pd.read_csv(txtFile,sep = ",",header=None,float_precision='high')
points = df1.values.tolist()

if(goal == "sym"):
    result = sy.sym(points)
    print_result(result)
elif(goal == "ddg"):
    result = sy.ddg(points)
    print_result(result)
elif(goal == "norm"):
    result = sy.norm(points)
    print_result(result)
elif (goal == "symnmf"):
    if (k >= len(points)):
        print("Invalid number of clusters!")
        exit()
    w = sy.norm(points)
    m = np.mean(w)
    n = len(points)
    np.random.seed(0)
    H = [[np.random.uniform(0,2*math.sqrt(m/k)) for _ in range(k)] for _ in range(n)]
    result = sy.symnmf(H,w,iter,EPS)
    print_result(result)

else:
    print("An error has occured")
    exit()



