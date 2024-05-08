#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double euclidean_distance(double *vector1, double* vector2,int length){
    int i;
    double dist = 0;
    for (i = 0;i<length;i++){
        dist += pow((vector1[i] - vector2[i]),2);
    }
    return dist;
}
double frobenius_norm_pow(double **matrix1,double **matrix2,int n,int k){
    int i,j;
    double sum=0.0;
    double res;
    for (i = 0;i<n;i++){
        for (j = 0;j<k;j++){
            res = matrix1[i][j] - matrix2[i][j];
            sum +=pow(res,2);
        }
    }
    return sum;
}
void symc (double **sym_matrix, double **points,int rows,int columns){
    int i,j;
    double data;
    for (i = 0;i<rows;i++){
        for (j = 0;j < rows;j++){
            if (i != j){
                data = (euclidean_distance(points[i],points[j],columns))/2;
                sym_matrix[i][j] = exp(-data);
            }
            else{
                sym_matrix[i][j] = 0;
            }
        }
    }
}
void ddgc (double **sym_matrix, double **diag_matrix,int rows){
    int i,j;
    double sum;
    for (i = 0;i<rows;i++){
        for (j = 0;j<rows;j++){
            diag_matrix[i][j] = 0;
            sum+= sym_matrix[i][j];
        }
        diag_matrix[i][i] = sum;
        sum =0;
    }
}
void normc (double **diag_matrix,double **sym_matrix,double **norm_matrix,int rows){
    int i,j,k;
    double sum;
    double *norm_repP, **norm_rep;
    norm_repP = calloc(rows*rows,sizeof(double));
    if (norm_repP == NULL){
        printf("An error has occured!\n");
    }
    norm_rep = calloc(rows,sizeof(double*));
    if (norm_rep == NULL){
        free(norm_repP);
        printf("An error has occured!\n");
    }
    for (i = 0 ;i<rows;i++){
        norm_rep[i] = i*rows + norm_repP;
    }
    for (i = 0;i<rows;i++){
        diag_matrix[i][i] = 1 / sqrt(diag_matrix[i][i]);
    }
    for (i = 0;i<rows;i++){
        for (j = 0;j<rows;j++){
            for (k =0;k<rows;k++){
                sum += diag_matrix[i][k] * sym_matrix[k][j];
            }
            norm_rep[i][j] = sum;
            sum = 0;
        }
    }
    for (i = 0;i<rows;i++){
        for (j = 0;j<rows;j++){
            for (k =0;k<rows;k++){
                sum += norm_rep[i][k] * diag_matrix[k][j];
            }
            norm_matrix[i][j] = sum;
            sum = 0;
        }
    }
    free(norm_repP);
    free(norm_rep);
}
void transpose (double **matrix, double **matrix_t,int n,int k){
    int i,j;
    for (i=0;i<n;i++){
        for (j =0;j<k;j++){
            matrix_t[j][i] = matrix[i][j];
        }
    }
}
void matrix_mult(double **mn_matrix,double **np_matrix,double **result,int m,int n,int p){
    int i,j,k;
    for(i=0;i<m;i++){
        for(j =0;j<p;j++){
            result[i][j] = 0;
            for (k=0;k<n;k++){
               result[i][j] += mn_matrix[i][k] * np_matrix[k][j]; 
            }
        }
    }
}
void symnmfc(double **h_matrix, double **norm_matrix, int n,int k,double EPS,int iter){
    double beta = 0.5,from_norm,num;
    double *h_copyP,*ht_matrixP,*wh_matrixP,*hht_matrixP,*hhth_matrixP;
    double **h_copy,**ht_matrix,**wh_matrix,**hht_matrix,**hhth_matrix;
    int i,j;
    h_copyP = calloc(n*k,sizeof(double));
    if (h_copyP == NULL){
        printf("An error has occured!\n");
    }
    h_copy = calloc(n,sizeof(double*));
    if (h_copy == NULL){
        free(h_copyP);
        printf("An error has occured\n");
    }
    wh_matrixP = calloc(n*k,sizeof(double));
    if (wh_matrixP == NULL){
        free(h_copyP);
        free(h_copy);
        printf("An error has occured!\n");
    }
    wh_matrix = calloc(n,sizeof(double*));
    if (wh_matrix == NULL){
        printf("An error has occured!\n");
        free(h_copyP);
        free(h_copy);
        free(wh_matrixP);
    }
    hhth_matrixP = calloc(n*k,sizeof(double));
    if (hhth_matrixP == NULL){
        printf("An error has occured!\n");
        free(h_copyP);
        free(h_copy);
        free(wh_matrixP);
        free(wh_matrix);
    }
    hhth_matrix = calloc(n,sizeof(double*));
    if (hhth_matrix == NULL){
        printf("An error has occured!\n");
        free(h_copyP);
        free(h_copy);
        free(wh_matrixP);
        free(wh_matrix);
        free(hhth_matrixP);
    }
    for (i = 0; i<n;i++){
        h_copy[i] = h_copyP + i*k;
        wh_matrix[i] = wh_matrixP + i*k;
        hhth_matrix[i] = hhth_matrixP + i*k;
    }
    hht_matrixP = calloc(n*n,sizeof(double));
    if (hht_matrixP == NULL){
        printf("An error has occured!\n");
        free(h_copyP);
        free(h_copy);
        free(wh_matrixP);
        free(wh_matrix);
        free(hhth_matrixP);
        free(hhth_matrix);
    }
    hht_matrix = calloc(n,sizeof(double*));
    if (hht_matrix == NULL){
        printf("An error has occured!\n");
        free(h_copyP);
        free(h_copy);
        free(wh_matrixP);
        free(wh_matrix);
        free(hhth_matrixP);
        free(hhth_matrix);
        free(hht_matrixP);
    }
    for (i = 0;i<n;i++){
        hht_matrix[i] = hht_matrixP + i*n;
    }
    ht_matrixP = calloc(k*n,sizeof(double));
    if (ht_matrixP == NULL){
        printf("An error has occured!\n");
        free(h_copyP);
        free(h_copy);
        free(wh_matrixP);
        free(wh_matrix);
        free(hhth_matrixP);
        free(hhth_matrix);
        free(hht_matrixP);
        free(hht_matrix);
    }
    ht_matrix = calloc(k,sizeof(double*));
    if (ht_matrixP == NULL){
        printf("An error has occured!\n");
        free(h_copyP);
        free(h_copy);
        free(wh_matrixP);
        free(wh_matrix);
        free(hhth_matrixP);
        free(hhth_matrix);
        free(hht_matrixP);
        free(hht_matrix);
        free(ht_matrixP);
    }
    for (i =0;i<k;i++){
        ht_matrix[i] = ht_matrixP + i*n;
    }
    while(iter > 0){
        for (i=0;i<n;i++){
            for (j =0;j<k;j++){
                h_copy[i][j] = h_matrix[i][j];
            }
        }
        matrix_mult(norm_matrix,h_matrix,wh_matrix,n,n,k);
        transpose(h_matrix,ht_matrix,n,k);
        matrix_mult(h_matrix,ht_matrix,hht_matrix,n,k,n);
        matrix_mult(hht_matrix,h_matrix,hhth_matrix,n,n,k);
        for (i=0;i<n;i++){
            for (j =0;j<k;j++){
                num = wh_matrix[i][j] / hhth_matrix[i][j];
                h_matrix[i][j] = h_copy[i][j] * (1-beta + (beta*num));
            }
        }
        from_norm = frobenius_norm_pow(h_matrix,h_copy,n,k);
        if (from_norm < EPS){
            break;
        }
        iter--;
    }
    free(h_copyP);
    free(h_copy);
    free(wh_matrixP);
    free(wh_matrix);
    free(hhth_matrixP);
    free(hhth_matrix);
    free(hht_matrixP);
    free(hht_matrix);
    free(ht_matrixP);
    free(ht_matrix);
}


int main(){
    /*
    double num1,num2;
    double *arrp,*diagp,*normp,*hp,*pointsP;
    double **arr,**diag,**normm,**h,**points;
    int i,j;
    arrp = calloc(10000,sizeof(double));
    arr = calloc(100,sizeof(double*));
    diagp = calloc(36,sizeof(double));
    diag = calloc(6,sizeof(double*));
    normp = calloc(36,sizeof(double));
    normm = calloc(6,sizeof(double*));
    pointsP = calloc(36,sizeof(double));
    points = calloc(6,sizeof(double*));
    hp = calloc(200,sizeof(double));
    h = calloc(100,sizeof(double*));
    for (i = 0;i<6;i++){
        diag[i] = diagp + 6*i;
        points[i] = pointsP + 6*i;
        normm[i] = normp + 6*i;
    }
    FILE *simFile;
    simFile = fopen("input1.txt","r");
    if (simFile == NULL){
        printf("error reading file\n");
    }
    for (i = 0;i<6;i++){
        for (j = 0;j<6;j++){
            fscanf(simFile,"%lf,",&num1);
            points[i][j] = num1;
        }
    }
    fclose(simFile); 
    ddgc(points,diag,6);
    normc(diag,points,normm,6);
    for (i = 0;i<6;i++){
        for (j=0;j<6;j++){
            printf("%.4f,",normm[i][j]); 
        }
        printf("\n");
    }   */
}