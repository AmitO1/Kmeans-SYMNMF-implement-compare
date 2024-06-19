# define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "symnmf.h"
static PyObject* norm(PyObject *self, PyObject *args){
    PyObject *points_list, *norma_matrix,*norma_row;
    PyObject *item1;
    int n_check,k_check;
    double num;
    double *sym_matrixP,*diag_matrixP,*norm_matrixP,*points_matrixP;
    double **sym_matrix,**diag_matrix,**norm_matrix,**points_matrix;
    int i,j;
    if (!PyArg_ParseTuple(args,"O",&points_list)){
        return NULL;
    }
    /*find vectors amount and dimension for memory allocating*/
    n_check = PyObject_Length(points_list);
    if (n_check < 0){
        return NULL;
    }
    k_check = PyObject_Length(PyList_GetItem(points_list,0));
    /*allocate memory for matrixed needed for send to C function*/
    sym_matrixP = calloc(n_check*n_check,sizeof(double));
    if (sym_matrixP == NULL){
        printf("An error has occured!\n");
        return NULL;
    }
    sym_matrix = calloc(n_check,sizeof(double*));
    if (sym_matrix == NULL){
        free(sym_matrixP);
        printf("An error has occured!\n");
        return NULL;
    }
    diag_matrixP = calloc(n_check*n_check,sizeof(double));
    if (diag_matrixP == NULL){
        free(sym_matrixP);
        free(sym_matrix);
        printf("An error has occured!\n");
        return NULL;
    }
    diag_matrix = calloc(n_check,sizeof(double*));
    if (diag_matrix == NULL){
        free(sym_matrixP);
        free(sym_matrix);
        free(diag_matrixP);
        printf("An error has occured!\n");
        return NULL;       
    }
    norm_matrixP = calloc(n_check*n_check,sizeof(double));
    if (norm_matrixP == NULL){
        free(sym_matrixP);
        free(sym_matrix);
        free(diag_matrixP);
        free(diag_matrix);
        printf("An error has occured!\n");
        return NULL;
    }
    norm_matrix = calloc(n_check,sizeof(double*));
    if (norm_matrix == NULL){
        free(sym_matrixP);
        free(sym_matrix);
        free(diag_matrixP);
        free(diag_matrix);
        free(norm_matrixP);
        printf("An error has occured!\n");
        return NULL;        
    }
    for (i = 0;i<n_check;i++){
        sym_matrix[i] = sym_matrixP + i*n_check;
        diag_matrix[i] = diag_matrixP + i*n_check;
        norm_matrix[i] = norm_matrixP + i*n_check;
    }
    points_matrixP = calloc(n_check*k_check,sizeof(double));
    if(points_matrixP == NULL){
        free(sym_matrixP);
        free(sym_matrix);
        free(diag_matrixP);
        free(diag_matrix);
        free(norm_matrixP);
        free(norm_matrix);
        printf("An error has occured!\n"); 
        return NULL;         
    }
    points_matrix = calloc(n_check,sizeof(double*));
    if (points_matrix == NULL){
        free(sym_matrixP);
        free(sym_matrix);
        free(diag_matrixP);
        free(diag_matrix);
        free(norm_matrixP);
        free(norm_matrix);
        free(points_matrixP);
        printf("An error has occured!\n");   
        return NULL;       
    }
    for (i =0;i<n_check;i++){
        points_matrix[i] = points_matrixP + i*k_check;
    }
    /*get vectors inside C matrix from vectors recieved from python*/
    for(i=0;i<n_check;i++){
        for(j =0;j<k_check;j++){
            item1 = PyList_GetItem(PyList_GetItem(points_list,i),j);
            num = PyFloat_AsDouble(item1);
            points_matrix[i][j] = num;
        }
    }
    /*run needed c function in order to get the matrixes needed to calculate norm matrix*/
    symc(sym_matrix,points_matrix,n_check,k_check);
    ddgc(sym_matrix,diag_matrix,n_check);
    normc(diag_matrix,sym_matrix,norm_matrix,n_check);
    /*transfer c matrix to python list and return norm matrix*/
    norma_matrix = PyList_New(n_check);
    for (i=0;i<n_check;i++){
        norma_row = PyList_New(n_check);
        for (j=0;j<n_check;j++){
            item1 = PyFloat_FromDouble(norm_matrix[i][j]);
            PyList_SetItem(norma_row,j,item1);
        }
        PyList_SetItem(norma_matrix,i,norma_row);
    }
    free(sym_matrixP);
    free(sym_matrix);
    free(diag_matrixP);
    free(diag_matrix);
    free(norm_matrixP);
    free(norm_matrix);
    free(points_matrixP); 
    free(points_matrix);
    return norma_matrix;   
}
static PyObject* sym(PyObject *self,PyObject *args){
    PyObject *points_list, *sy_matrix,*sy_row;
    PyObject *item1;
    int n_check,k_check;
    double num;
    double *sym_matrixP,*points_matrixP;
    double **sym_matrix,**points_matrix;
    int i,j;
    if (!PyArg_ParseTuple(args,"O",&points_list)){
        return NULL;
    }
    /*finding the needed matrix size according to the vectors recieved from python*/
    n_check = PyObject_Length(points_list);
    if (n_check < 0){
        return NULL;
    }
    k_check = PyObject_Length(PyList_GetItem(points_list,0));
    /*allocate memory for matrixes needed for the similarity matrix*/
    sym_matrixP = calloc(n_check*n_check,sizeof(double));
    if(sym_matrixP == NULL){
        printf("An error has occured!\n");
        return NULL;
    }
    sym_matrix = calloc(n_check,sizeof(double*));
    if (sym_matrix == NULL){
        free(sym_matrixP);
        printf("An error has occured!\n");
        return NULL;
    }
    for (i = 0;i<n_check;i++){
        sym_matrix[i] = sym_matrixP + i*n_check;
    }
    points_matrixP = calloc(n_check*k_check,sizeof(double));
    if (points_matrixP == NULL){
        free(sym_matrixP);
        free(sym_matrix);
        printf("An error has occured!\n");
        return NULL;
    }
    points_matrix = calloc(n_check,sizeof(double*));
    if (points_matrix == NULL){
        free(sym_matrixP);
        free(sym_matrix);
        free(points_matrixP);
        printf("An error has occured!\n");
        return NULL;
    }
    for (i =0;i<n_check;i++){
        points_matrix[i] = points_matrixP + i*k_check;
    }
    /*update the points matrix according to the vecotrs list recived from python*/
    for(i=0;i<n_check;i++){
        for(j =0;j<k_check;j++){
            item1 = PyList_GetItem(PyList_GetItem(points_list,i),j);
            num = PyFloat_AsDouble(item1);
            points_matrix[i][j] = num;
        }
    }
    /*run the symc function from c and then update the sy_matrix in order to return for python*/
    symc(sym_matrix,points_matrix,n_check,k_check);
    sy_matrix = PyList_New(n_check);
    for (i=0;i<n_check;i++){
        sy_row = PyList_New(n_check);
        for (j=0;j<n_check;j++){
            item1 = PyFloat_FromDouble(sym_matrix[i][j]);
            PyList_SetItem(sy_row,j,item1);
        }
        PyList_SetItem(sy_matrix,i,sy_row);
    }
    free(sym_matrixP);
    free(sym_matrix);
    free(points_matrixP);
    free(points_matrix);
    return sy_matrix;

}
static PyObject* ddg(PyObject *self, PyObject *args){
    PyObject *points_list;
    PyObject *w_matrix,*w_row;
    PyObject *item1;
    int n_check,k_check,i,j;
    double num;
    double *diag_matrixP,*sym_matrixP,*points_matrixP;
    double **diag_matrix,**sym_matrix,**points_matrix;
    if (!PyArg_ParseTuple(args,"O",&points_list)){
        return NULL;
    }
    /*find vectors dimension and amount*/
    n_check = PyObject_Length(points_list);
    if (n_check < 0){
        return NULL;
    }
    k_check = PyObject_Length(PyList_GetItem(points_list,0));
    /*allocate memory needed for each matrix*/
    diag_matrixP = calloc(n_check*n_check,sizeof(double));
    if (diag_matrixP == NULL){
        printf("An error has occured!\n");
        return NULL;
    }
    diag_matrix = calloc(n_check,sizeof(double*));
    if (diag_matrix == NULL){
        free(diag_matrixP);
        printf("An error has occured!\n");
        return NULL;
    }
    sym_matrixP = calloc(n_check*n_check,sizeof(double));
    if (sym_matrixP == NULL){
        free(diag_matrixP);
        free(diag_matrix);
        printf("An error has occured!\n");
        return NULL;
    }
    sym_matrix = calloc(n_check,sizeof(double*));
    if (sym_matrix == NULL){
        free(diag_matrixP);
        free(diag_matrix);
        free(sym_matrixP);
        printf("An error has occured!\n");
        return NULL;
    }
    for(i=0;i<n_check;i++){
        diag_matrix[i] = diag_matrixP + i*n_check;
        sym_matrix[i] = sym_matrixP + i*n_check;
    }
    points_matrixP = calloc(n_check*k_check,sizeof(double));
    if (points_matrixP == NULL){
        free(diag_matrixP);
        free(diag_matrix);
        free(sym_matrixP);
        free(sym_matrix);
        printf("An error has occured!\n");
        return NULL;
    }
    points_matrix = calloc(n_check,sizeof(double*));
    if (points_matrix == NULL){
        free(diag_matrixP);
        free(diag_matrix);
        free(sym_matrixP);
        free(sym_matrix);
        free(points_matrixP);
        printf("An error has occured!\n");
        return NULL;
    }
    for(i=0;i<n_check;i++){
        points_matrix[i] = points_matrixP + i*k_check;
    }
    /*fill points matrix from points recieved from python list*/
    for(i=0;i<n_check;i++){
        for(j =0;j<k_check;j++){
            item1 = PyList_GetItem(PyList_GetItem(points_list,i),j);
            num = PyFloat_AsDouble(item1);
            points_matrix[i][j] = num;
        }
    }
    /*calculate the similarity matrix and diagonal matrix*/
    symc(sym_matrix,points_matrix,n_check,k_check);
    ddgc(sym_matrix,diag_matrix,n_check);
    /*transfer the diagonal matrix from C matrix to python list and return it*/
    w_matrix = PyList_New(n_check);
    for (i=0;i<n_check;i++){
        w_row = PyList_New(n_check);
        for (j=0;j<n_check;j++){
            item1 = PyFloat_FromDouble(diag_matrix[i][j]);
            PyList_SetItem(w_row,j,item1);
        }
        PyList_SetItem(w_matrix,i,w_row);
    }
    free(diag_matrixP);
    free(diag_matrix);
    free(sym_matrixP);
    free(sym_matrix);
    free(points_matrixP);
    free(points_matrix);
    return w_matrix;
}
static PyObject* symnmf(PyObject *self, PyObject *args){
    PyObject *h_matrix, *w_matrix;
    PyObject *item1;
    double EPS,num;
    int iter,i,j,n_check,k_check,check2;
    double *h_matP,*w_matP;
    double **h_mat,**w_mat;
    if (!PyArg_ParseTuple(args,"OOid",&h_matrix,&w_matrix,&iter,&EPS)){
        return NULL;
    }
    /*find vectors amount and dimension*/
    n_check = PyObject_Length(h_matrix);
    check2 = PyObject_Length(w_matrix);
    if (n_check < 0 || check2 < 0){
        return NULL;
    }
    k_check = PyObject_Length(PyList_GetItem(h_matrix,0));
    /*allocate memory for needed matrixes for the calculation*/
    h_matP = calloc(n_check*k_check,sizeof(double));
    if (h_matP == NULL){
        printf("An error has occured!\n");
        return NULL;
    }
    h_mat = calloc(n_check,sizeof(double*));
    if (h_mat == NULL){
        free(h_matP);
        printf("An error has occured!\n");
        return NULL;
    }
    for (i =0;i<n_check;i++){
        h_mat[i] = h_matP +i*k_check;
    }
    w_matP = calloc(n_check*n_check,sizeof(double*));
    if (w_matP == NULL){
        free(h_matP);
        free(h_mat);
        printf("An error has occured!\n");
        return NULL;
    }
    w_mat = calloc(n_check,sizeof(double*));
    if (w_mat == NULL){
        free(h_matP);
        free(h_mat);
        free(w_matP);
        printf("An error has occured!\n");
        return NULL;
    }
    for (i =0;i<n_check;i++){
        w_mat[i] = w_matP + i*n_check;
    }
    /*tansfer initial H and W values from python list to C matrix and find the final H*/
    for (i =0;i<n_check;i++){
        for (j = 0;j<n_check;j++){
            item1 = PyList_GetItem(PyList_GetItem(w_matrix,i),j);
            num = PyFloat_AsDouble(item1);
            w_mat[i][j] = num;
        }
    }
    for (i =0;i<n_check;i++){
        for (j = 0;j<k_check;j++){
            item1 = PyList_GetItem(PyList_GetItem(h_matrix,i),j);
            num = PyFloat_AsDouble(item1);
            h_mat[i][j] = num;
        }
    }
    symnmfc(h_mat,w_mat,n_check,k_check,EPS,iter);
    /*transfer final H from C matrix to python list*/
    for (i =0;i<n_check;i++){
        for (j = 0;j<k_check;j++){
            item1 = PyFloat_FromDouble(h_mat[i][j]);
            PyList_SetItem(PyList_GetItem(h_matrix,i),j,item1);
        }
    }
    free(h_matP);
    free(h_mat);
    free(w_matP);
    free(w_mat);
    return h_matrix;        
}
static PyMethodDef mysymnmfMethods[] = {
    {
        "ddg",
        ddg,
        METH_VARARGS,
        "the function calculate the diagonle matrix"
    }, {
        "sym",
        sym,
        METH_VARARGS,
        "calculate the symatrical matrix"
    }, {
        "norm",
        norm,
        METH_VARARGS,
        "calculate the norm matrix"
    }, {
        "symnmf",
        symnmf,
        METH_VARARGS,
        "the complete symnmf algorithm"
    }
    ,{NULL,NULL,0,NULL}
};
static struct PyModuleDef mysymnmfModule = {
    PyModuleDef_HEAD_INIT,
    "mysymnmf",
    "need to check what to do",
    -1,
    mysymnmfMethods
};
PyMODINIT_FUNC PyInit_mysymnmf(void){
    return PyModule_Create(&mysymnmfModule);
}