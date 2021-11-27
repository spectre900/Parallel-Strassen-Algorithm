%%cu
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <bits/stdc++.h>

using namespace std;

void print(int n, int* mat)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << mat[i * n + j] << " ";
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

void fillMatrix(int n, int*& mat)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mat[i * n + j] = rand() % 5;
        }
    }
}

void freeMatrix(int n, int* mat)
{
    free(mat);
}

__global__ void matrixMultiplication(int* mat1, int* mat2, int* product, int n)
{
    int prod = blockIdx.x * blockDim.x + threadIdx.x;
    int i = prod / n;
    int j = prod % n;
    for (int k = 0; k < n; k++) {
        product[i * n + j] += mat1[i * n + k] * mat2[k * n + j];
    }
}

int main()
{
    int n;
    cin >> n;
    n = 1024;

    int* h_mat1 = allocateMatrix(n);
    fillMatrix(n, h_mat1);

    int* h_mat2 = allocateMatrix(n);
    fillMatrix(n, h_mat2);

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

    clock_t start, end;
    start = clock();

    matrixMultiplication<<<gridSize, blockSize>>>(d_mat1, d_mat2, d_product, n);
    cudaDeviceSynchronize();

    end = clock();
    double time = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Parallel Naive Runtime (CUDA): " << time << " seconds\n";

    cudaMemcpy(h_product, d_product, bytes, cudaMemcpyDeviceToHost);

    cudaFree(d_mat1);
    cudaFree(d_mat2);
    cudaFree(d_product);

    freeMatrix(n, h_mat1);
    freeMatrix(n, h_mat2);
    freeMatrix(n, h_product);

    return 0;
}