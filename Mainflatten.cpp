
#include <iostream>
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <math.h>
#include "Functionsflatten.cpp"
#include <time.h>

using namespace std;
 
const int n = 16;
const int m = 50;
const double d1 = 1e-2;
const double tol = 1e-4;

/*****************************************************
* INPUTS:                                            *
* d1  : precision of Rutishauser                     *
* tol : precision of Inverse Iteration               *
* m   : maximum number of iterations (integer)       *
* n   : order of matrix A (integer)                  *
* A   : input matrix to study (of MAT type)          *
* -------------------------------------------------- *
* OUTPUTS:                                           *
* R   : table of eigenvalues (of VEC type)           *
* VX  : table of eigenvectors (of MAT type)          *
*****************************************************/

// Driver code
int main()
{    
    //allocate vector and matrices
    double   *Mat = new double [n * n];                    // given square symmetric matrix
    double   *Eigenvalues = new double [n];                // eigenvalues solution vector
    double   *Eigenvectors = new double [n * n];           // eigenvectors solution matrix   
    double   *lambda = new double [n]; 

    double A[][16] = {
    {3.7450,    0.0839,    0.8675,    0.8287,    2.0374,    1.7158,   -0.8208,    2.3104,         0,    0.6090,    0.8836,   -2.2513,    1.1669,    0.1956,    0.6991,   -1.3373},
    {0.0839,    3.0862,    0.2318,    0.5619,    0.5974,    1.7734,    1.5466,   -0.7386,   -0.6090,         0,    0.6258,    0.8013,   -1.6870,    0.8037,    0.3392,    0.6059},
    {0.8675,    0.2318,    3.3715,    0.0973,    0.8360,    0.8869,    1.9295,    1.5585,   -0.8836,   -0.6258,         0,    0.3838,    0.7301,   -1.9622,    1.0107,   -0.0538},
    {0.8287,    0.5619,    0.0973,    3.5665,   -0.0063,    0.6810,    0.6569,    2.0806,    2.2513,   -0.8013,   -0.3838,         0,    0.8428,    0.8335,    -1.8880,    1.2303},
    {2.0374,    0.5974,    0.8360,   -0.0063,    3.1479,    0.2583,    0.6119,    0.7910,   -1.1669,    1.6870,   -0.7301,   -0.8428,         0,    0.3315,    0.7888,   -2.0354},
    {1.7158,    1.7734,    0.8869,    0.6810,    0.2583,    3.2225,    0.1897,    0.8974,   -0.1956,   -0.8037,    1.9622,   -0.8335,   -0.3315,         0,    0.7787,    0.5938},
   {-0.8208,    1.5466,    1.9295,    0.6569,    0.6119,    0.1897,    3.3317,    0.0532,   -0.6991,   -0.3392,   -1.0107,    1.8880,   -0.7888,   -0.7787,         0,    0.5647},    
    {2.3104,   -0.7386,    1.5585,    2.0806,    0.7910,    0.8974,    0.0532,    3.5365,    1.3373,   -0.6059,    0.0538,   -1.2303,    2.0354,   -0.5938,    -0.5647,         0},
    {     0,   -0.6090,   -0.8836,    2.2513,   -1.1669,   -0.1956,   -0.6991,    1.3373,    3.7450,    0.0839,    0.8675,    0.8287,    2.0374,    1.7158,    -0.8208,    2.3104},
    {0.6090,         0,   -0.6258,   -0.8013,    1.6870,   -0.8037,   -0.3392,   -0.6059,    0.0839,    3.0862,    0.2318,    0.5619,    0.5974,    1.7734,    1.5466,   -0.7386},
    {0.8836,    0.6258,         0,   -0.3838,   -0.7301,    1.9622,   -1.0107,    0.0538,    0.8675,    0.2318,    3.3715,    0.0973,    0.8360,    0.8869,    1.9295,    1.5585},
   {-2.2513,    0.8013,    0.3838,         0,   -0.8428,   -0.8335,    1.8880,   -1.2303,    0.8287,    0.5619,    0.0973,    3.5665,   -0.0063,    0.6810,    0.6569,    2.0806},
    {1.1669,   -1.6870,    0.7301,    0.8428,         0,   -0.3315,   -0.7888,    2.0354,    2.0374,    0.5974,    0.8360,   -0.0063,    3.1479,    0.2583,    0.6119,    0.7910},
    {0.1956,    0.8037,   -1.9622,    0.8335,    0.3315,         0,   -0.7787,   -0.5938,    1.7158,    1.7734,    0.8869,    0.6810,    0.2583,    3.2225,    0.1897,    0.8974},
    {0.6991,    0.3392,    1.0107,   -1.8880,    0.7888,    0.7787,         0,   -0.5647,   -0.8208,    1.5466,    1.9295,    0.6569,    0.6119,    0.1897,    3.3317,    0.0532},
   {-1.3373,    0.6059,   -0.0538,    1.2303,   -2.0354,    0.5938,    0.5647,         0,    2.3104,   -0.7386,    1.5585,    2.0806,    0.7910,    0.8974,    0.0532,    3.5365}
    };
    
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++)
        {
            Mat[i * n + j] = A[i][j];
        }
    }
    for (int i = 0; i < n; i++) 
    {
        lambda[i] = 0;
    }

    //Displaying the orignal matrix
    cout << "---------------------------------------------------\n" << endl;
    cout <<  "Symmetric matrix is:" << endl;
    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++)
        {
            cout << Mat[i * n + j] << "\t";
        }
        cout << "\t" << endl;
    }
    cout << "\t" << endl;  
    
    clock_t t;
    t=clock();
    CAL(Mat, Eigenvalues, Eigenvectors, n, d1, m, tol, lambda);                //InverseIteration.cpp
    t = clock() - t;
    printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    
    //Displaying the results 
    cout << "Eigenvalues(Rutishauser):" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << Eigenvalues[i] << "\t";
    }
    cout << "\t" << endl;
    cout << "\t" << endl;

    cout << "Eigenvalues (Inverse Iteration):" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << lambda[i] << "\t";
    }
    cout << "\t" << endl;
    cout << "\t" << endl;

    cout << "Eigenvectors (in lines):" << endl;
    for (int j = 0; j < n; j++) 
    {
        for (int i = 0; i < n; i++)
        {
            cout << Eigenvectors[i * n + j] << "\t";
        }
	    cout << "\t" << endl;
    }
    cout << "\t" << endl;

    delete [] Mat;
    delete [] Eigenvalues;
    delete [] Eigenvectors;
    delete [] lambda;

    return 0;
}

