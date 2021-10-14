#include <omp.h>
#include <bits/stdc++.h>

using namespace std;

void print(int n, int** mat)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << mat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int** allocateMatrix(int n)
{
    int** mat = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++)
    {
        mat[i] = (int*)malloc(n * sizeof(int));
        for (int j = 0; j < n; j++)
        {
            mat[i][j] = rand() % 5;
        }
    }
    return mat;
}

void freeMatrix(int n, int** mat)
{
    for (int i = 0; i < n; i++)
        delete[] mat[i];
    delete[] mat;
}

int** naive(int n, int** mat1, int** mat2)
{
    int** prod = allocateMatrix(n);

    int i, j;

    #pragma omp parallel for collapse(2)
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
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

int** getSlice(int n, int** mat, int offseti, int offsetj)
{
    int m = n / 2;
    int** slice = allocateMatrix(m);
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
    int** result = allocateMatrix(n);
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
    int** result = allocateMatrix(n);

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

int** strassen(int n, int** mat1, int** mat2)
{

    if (n <= 512)
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

    int** s1;
    #pragma omp task shared(s1)
    {
        int** bds = addMatrices(m, b, d, false);
        int** gha = addMatrices(m, g, h, true);
        s1 = strassen(m, bds, gha);
        freeMatrix(m, bds);
        freeMatrix(m, gha);
    }

    int** s2;
    #pragma omp task shared(s2)
    {
        int** ada = addMatrices(m, a, d, true);
        int** eha = addMatrices(m, e, h, true);
        s2 = strassen(m, ada, eha);
        freeMatrix(m, ada);
        freeMatrix(m, eha);
    }

    int** s3;
    #pragma omp task shared(s3)
    {
        int** acs = addMatrices(m, a, c, false);
        int** efa = addMatrices(m, e, f, true);
        s3 = strassen(m, acs, efa);
        freeMatrix(m, acs);
        freeMatrix(m, efa);
    }

    int** s4;
    #pragma omp task shared(s4)
    {
        int** aba = addMatrices(m, a, b, true);
        s4 = strassen(m, aba, h);
        freeMatrix(m, aba);
    }

    int** s5;
    #pragma omp task shared(s5)
    {
        int** fhs = addMatrices(m, f, h, false);
        s5 = strassen(m, a, fhs);
        freeMatrix(m, fhs);
    }

    int** s6;
    #pragma omp task shared(s6)
    {
        int** ges = addMatrices(m, g, e, false);
        s6 = strassen(m, d, ges);
        freeMatrix(m, ges);
    }

    int** s7;
    #pragma omp task shared(s7)
    {
        int** cda = addMatrices(m, c, d, true);
        s7 = strassen(m, cda, e);
        freeMatrix(m, cda);
    }

    #pragma omp taskwait

    freeMatrix(m, a);
    freeMatrix(m, b);
    freeMatrix(m, c);
    freeMatrix(m, d);
    freeMatrix(m, e);
    freeMatrix(m, f);
    freeMatrix(m, g);
    freeMatrix(m, h);

    int** c11;
    #pragma omp task shared(c11)
    {
        int** s1s2a = addMatrices(m, s1, s2, true);
        int** s6s4s = addMatrices(m, s6, s4, false);
        c11 = addMatrices(m, s1s2a, s6s4s, true);
        freeMatrix(m, s1s2a);
        freeMatrix(m, s6s4s);
    }

    int** c12;
    #pragma omp task shared(c12)
    {
        c12 = addMatrices(m, s4, s5, true);
    }

    int** c21;
    #pragma omp task shared(c21)
    {
        c21 = addMatrices(m, s6, s7, true);
    }

    int** c22;
    #pragma omp task shared(c22)
    {
        int** s2s3s = addMatrices(m, s2, s3, false);
        int** s5s7s = addMatrices(m, s5, s7, false);
        c22 = addMatrices(m, s2s3s, s5s7s, true);
        freeMatrix(m, s2s3s);
        freeMatrix(m, s5s7s);
    }

    #pragma omp taskwait

    freeMatrix(m, s1);
    freeMatrix(m, s2);
    freeMatrix(m, s3);
    freeMatrix(m, s4);
    freeMatrix(m, s5);
    freeMatrix(m, s6);
    freeMatrix(m, s7);

    int** prod = combineMatrices(m, c11, c12, c21, c22);

    freeMatrix(m, c11);
    freeMatrix(m, c12);
    freeMatrix(m, c21);
    freeMatrix(m, c22);

    return prod;
}

bool check(int n, int** prod1, int** prod2)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (prod1[i][j] != prod2[i][j])
                return false;
        }
    }
    return true;
}

int main()
{
    int n;
    cout << "\nEnter matrix Dimension: ";
    cin >> n;

    int** mat1 = allocateMatrix(n);
    int** mat2 = allocateMatrix(n);

    double startParStrassen = omp_get_wtime();
    int** prod;
    #pragma omp parallel
    {
        #pragma omp single
        {
            prod = strassen(n, mat1, mat2);
        }
    }
    double endParStrassen = omp_get_wtime();
    cout << "\nParallel Strassen Runtime: ";
    cout << setprecision(5) << endParStrassen - startParStrassen << endl;

    cout << endl;

    return 0;
}