import math
import sys
import os
import pandas as pd
import mysymnmf as sy
import numpy as np
from sklearn.metrics import silhouette_score

#function to assign the correct label for each cluster after kmeans algorithm
def assign_labels(X, cluster_centers):
    X = np.array(X)
    cluster_centers = np.array(cluster_centers)
    n_samples = X.shape[0]
    n_clusters = cluster_centers.shape[0]
    labels = np.empty(n_samples, dtype=int)

    for i in range(n_samples):
        distances = np.linalg.norm(X[i] - cluster_centers, axis=1)
        labels[i] = np.argmin(distances)

    return labels

def print_result(result):
    for i in range(len(result)):
        for j in range(len(result[i])):
            if (j != len(result[i]) -1):
                sys.stdout.write(f"{result[i][j]:.4f},")
            else:
                sys.stdout.write(f"{result[i][j]:.4f}")
        print()

#used in Kmeans
def isNatural(number):
        if (number - math.floor(number) == 0):
                return True
        return False
#used in Kmeans
def euclideanDistance(vector1,vector2):
        sum =0
        distance =0
        for i in range(len(vector1)):
                distance = (vector1[i]-vector2[i])
                sum += math.pow(distance,2)
        return math.sqrt(sum)
#used in Kmeans
def closestCluster (kCluster,vector):
        index = -1
        minDis = 0
        temp =0
        for i in range(len(kCluster)):
                if (index == -1):
                        index = i
                        minDis = euclideanDistance(kCluster[i],vector)
                else:   
                        temp = euclideanDistance(kCluster[i],vector)
                        if (temp < minDis):
                                minDis = temp
                                index = i
        return index
                
num_args = len(sys.argv) - 1
EPS = 0.0001
iterNum = 300

if(num_args == 2):
        K = (int)(sys.argv[1])
        txtFile = sys.argv[2]
        check = os.path.isfile(txtFile)
        if not (check):
                print("An Error Has Occurred")
                exit()
else:
        print("An Error Has Occurred")
        exit()

df1 = pd.read_csv(txtFile,sep = ",",header=None,float_precision='high')
vectorsArray = df1.values.tolist()

d = len(vectorsArray[0])
N = len(vectorsArray)

if (K >= N):
       print("An error has occured")
       exit()   
#run the Kmeans algorithm from assignment 1 
distanceSmallerThenEPS = 0
kClusters = [[0] * d for _ in range(K)]
num_rows = len(kClusters)
num_cols = len(kClusters[0])
tempClusters = [[0] * num_cols for _ in range(num_rows)]
totalPoints = [0] * K
w = sy.norm(vectorsArray)
m = np.mean(w)
np.random.seed(0)
H = [[np.random.uniform(0,2*math.sqrt(m/K)) for _ in range(K)] for _ in range(N)]
H_result = sy.symnmf(H,w,iterNum,EPS)
for i in range(K):
        for j in range(d):
                kClusters[i][j] = vectorsArray[i][j]
                
while(iterNum > 0 and distanceSmallerThenEPS == 0):
        
        for i in range(N):
                clusterIndex = closestCluster(kClusters,vectorsArray[i])
                totalPoints[clusterIndex]+=1
                for j in range(d):
                        tempClusters[clusterIndex][j] += vectorsArray[i][j]
        for i in range(K):
                for j in range(d):
                        tempClusters[i][j] = (tempClusters[i][j] / totalPoints[i])
        
        distanceSmallerThenEPS = 1
        for i in range(K):
                
                if (euclideanDistance(kClusters[i],tempClusters[i]) >= EPS):
                        distanceSmallerThenEPS = 0
                if (totalPoints[i] != 0):
                        for j in range(d):
                                kClusters[i][j] =0
                                kClusters[i][j] = tempClusters[i][j]
                                tempClusters[i][j] = 0
                totalPoints[i] = 0
        iterNum -=1

#find label for Kmeans and Symnmf and calculate the silhouette score     
cluster_labels_kmean = assign_labels(vectorsArray,kClusters)
kmeans_silhouette_score = silhouette_score(vectorsArray, cluster_labels_kmean)
kmeans_silhouette_score = f"{kmeans_silhouette_score:.4f}"
cluster_labels_H = np.argmax(H_result, axis=1)
symnmf_silhouette_score = silhouette_score(vectorsArray, cluster_labels_H)
symnmf_silhouette_score = f"{symnmf_silhouette_score:.4f}"
print("nmf:", symnmf_silhouette_score)
print("kmeans:", kmeans_silhouette_score)
        

