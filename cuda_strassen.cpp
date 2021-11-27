%%cu
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <bits/stdc++.h>

using namespace std;

void print(int n, int** mat)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int* allocateMatrix(int n)
{
    int* data = (int*)malloc(n * n * sizeof(int));
    return data;
}

int** allocateMatrix2D(int n)
{
    int* data = (int*)malloc(n * n * sizeof(int));
    int** array = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++)
    {
        array[i] = &(data[n * i]);
    }
    return array;
}

void fillMatrix(int n, int*& mat)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat[i * n + j] = rand() % 5;
        }
    }
}

void fillMatrix2D(int n, int** &mat)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat[i][j] = rand() % 5;
        }
    }
}

int** getSlice(int n, int** mat, int offseti, int offsetj)
{
    int m = n / 2;
    int** slice = allocateMatrix2D(m);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            slice[i][j] = mat[offseti + i][offsetj + j];
        }
    }
    return slice;
}

int** addMatrices(int n, int** mat1, int** mat2, bool add)
{
    int** result = allocateMatrix2D(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (add)
                result[i][j] = mat1[i][j] + mat2[i][j];
            else
                result[i][j] = mat1[i][j] - mat2[i][j];
        }
    }

    return result;
}

int** combineMatrices(int m, int** c11, int** c12, int** c21, int** c22)
{
    int n = 2 * m;
    int** result = allocateMatrix2D(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i < m && j < m)
                result[i][j] = c11[i][j];
            else if (i < m)
                result[i][j] = c12[i][j - m];
            else if (j < m)
                result[i][j] = c21[i - m][j];
            else
                result[i][j] = c22[i - m][j - m];
        }
    }

    return result;
}

void freeMatrix(int n, int* mat)
{
    free(mat);
}

void freeMatrix2D(int n, int** mat)
{
    free(mat[0]);
    free(mat);
}

int** naive(int n, int** mat1, int** mat2)
{
    int** prod = allocateMatrix2D(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            prod[i][j] = 0;
            for (int k = 0; k < n; k++)
            {
                prod[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return prod;
}

__global__ void multiply(int* mat1, int* mat2, int* product, int n)
{
    int prod = blockIdx.x * blockDim.x + threadIdx.x;
    int i = prod / n;
    int j = prod % n;
    for (int k = 0; k < n; k++) {
        product[i * n + j] += mat1[i * n + k] * mat2[k * n + j];
    }
}

int** cudaNaive(int n, int** mat1, int** mat2)
{
    int* h_mat1 = allocateMatrix(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            h_mat1[i*n + j] = mat1[i][j];
        }
    }

    int* h_mat2 = allocateMatrix(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            h_mat2[i*n + j] = mat2[i][j];
        }
    }

    int* h_product = allocateMatrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            h_product[i * n + j] = 0;
        }
    }

    size_t bytes = n * n * sizeof(int);

    int *d_mat1, *d_mat2, *d_product;

    cudaMalloc(&d_mat1, bytes);
    cudaMalloc(&d_mat2, bytes);
    cudaMalloc(&d_product, bytes);

    cudaMemcpy(d_mat1, h_mat1, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_mat2, h_mat2, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_product, h_product, bytes, cudaMemcpyHostToDevice);

    int threads = min(1024, n);
    int blocks = (n * n) / threads;
    dim3 gridSize(blocks, 1, 1);
    dim3 blockSize(threads, 1, 1);

    multiply<<<gridSize, blockSize>>>(d_mat1, d_mat2, d_product, n);
    cudaDeviceSynchronize();

    cudaMemcpy(h_product, d_product, bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_mat1);
    cudaFree(d_mat2);
    cudaFree(d_product);

    freeMatrix(n, h_mat1);
    freeMatrix(n, h_mat2);

    int** product = allocateMatrix2D(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            product[i][j] = h_product[i*n + j];
        }
    }
    return product;
}


int** strassen(int n, int** mat1, int** mat2)
{

    if (n <= 32)
    {
        return naive(n, mat1, mat2);
    }

    int m = n / 2;

    int** a = getSlice(n, mat1, 0, 0);
    int** b = getSlice(n, mat1, 0, m);
    int** c = getSlice(n, mat1, m, 0);
    int** d = getSlice(n, mat1, m, m);
    int** e = getSlice(n, mat2, 0, 0);
    int** f = getSlice(n, mat2, 0, m);
    int** g = getSlice(n, mat2, m, 0);
    int** h = getSlice(n, mat2, m, m);

    int** bds = addMatrices(m, b, d, false);
    int** gha = addMatrices(m, g, h, true);
    int** s1 = cudaNaive(m, bds, gha);
    freeMatrix2D(m, bds);
    freeMatrix2D(m, gha);

    int** ada = addMatrices(m, a, d, true);
    int** eha = addMatrices(m, e, h, true);
    int** s2 = cudaNaive(m, ada, eha);
    freeMatrix2D(m, ada);
    freeMatrix2D(m, eha);

    int** acs = addMatrices(m, a, c, false);
    int** efa = addMatrices(m, e, f, true);
    int** s3 = cudaNaive(m, acs, efa);
    freeMatrix2D(m, acs);
    freeMatrix2D(m, efa);

    int** aba = addMatrices(m, a, b, true);
    int** s4 = cudaNaive(m, aba, h);
    freeMatrix2D(m, aba);
    freeMatrix2D(m, b);

    int** fhs = addMatrices(m, f, h, false);
    int** s5 = cudaNaive(m, a, fhs);
    freeMatrix2D(m, fhs);
    freeMatrix2D(m, a);
    freeMatrix2D(m, f);
    freeMatrix2D(m, h);

    int** ges = addMatrices(m, g, e, false);
    int** s6 = cudaNaive(m, d, ges);
    freeMatrix2D(m, ges);
    freeMatrix2D(m, g);

    int** cda = addMatrices(m, c, d, true);
    int** s7 = cudaNaive(m, cda, e);
    freeMatrix2D(m, cda);
    freeMatrix2D(m, c);
    freeMatrix2D(m, d);
    freeMatrix2D(m, e);

    int** s1s2a = addMatrices(m, s1, s2, true);
    int** s6s4s = addMatrices(m, s6, s4, false);
    int** c11 = addMatrices(m, s1s2a, s6s4s, true);
    freeMatrix2D(m, s1s2a);
    freeMatrix2D(m, s6s4s);
    freeMatrix2D(m, s1);

    int** c12 = addMatrices(m, s4, s5, true);
    freeMatrix2D(m, s4);

    int** c21 = addMatrices(m, s6, s7, true);
    freeMatrix2D(m, s6);

    int** s2s3s = addMatrices(m, s2, s3, false);
    int** s5s7s = addMatrices(m, s5, s7, false);
    int** c22 = addMatrices(m, s2s3s, s5s7s, true);
    freeMatrix2D(m, s2s3s);
    freeMatrix2D(m, s5s7s);
    freeMatrix2D(m, s2);
    freeMatrix2D(m, s3);
    freeMatrix2D(m, s5);
    freeMatrix2D(m, s7);

    int** prod = combineMatrices(m, c11, c12, c21, c22);

    freeMatrix2D(m, c11);
    freeMatrix2D(m, c12);
    freeMatrix2D(m, c21);
    freeMatrix2D(m, c22);

    return prod;
}

int main()
{
    int n;
    n = 2;

    int** mat1 = allocateMatrix2D(n);
    fillMatrix2D(n, mat1);

    int** mat2 = allocateMatrix2D(n);
    fillMatrix2D(n, mat2);

    clock_t start, end;
    start = clock();

    int** prod = strassen(n, mat1, mat2);

    end = clock();
    double time = double(end - start) / double(CLOCKS_PER_SEC);
    cout<<"Parallel Strassen Runtime (CUDA): "<<time<<" seconds\n";

    return 0;
}