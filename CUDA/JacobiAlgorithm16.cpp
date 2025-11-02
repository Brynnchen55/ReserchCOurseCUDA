#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>     /*  for  fabs, sqrt, pow, exp, sin, cos, log,   */
                      /*       atan, acos                             */
#include <string.h>
#include <algorithm>
//#include <map>
#include <time.h>

using namespace std;
const int n = 16;


void Jacobi(double **A, int N, double *D, double **V, int *NROT) 
{
    double  c, g, h, s, sm, t, tau, theta, tresh;
    int     i, j, ip, iq;
    

    //allocate vectors B, Z
    double   *B = new double [n];
    double   *Z = new double [n];

    for (ip = 0; ip < N; ip++)         //initialize V to identity matrix
    {
        for (iq = 0; iq < N; iq++)
        {
            V[ip][iq] = 0.0; 
        }
        V[ip][ip] = 1.0;
    } 

    for (ip = 0; ip < N; ip++) 
    {
        B[ip] = A[ip][ip];
        D[ip] = B[ip];
        Z[ip] = 0.0;
    }

    *NROT = 0;

    for (i = 0; i < 50; i++) 
    {
        //i = 50;
        sm = 0;
        for (ip = 0; ip < N - 1; ip++)    //sum off-diagonal elements
        {
            for (iq = ip + 1; iq < N; iq++)
            {
                sm = sm + fabs(A[ip][iq]);
            }
        }

        if (sm == 0)
        {
            delete [] B;
            delete [] Z; 
            return;       //normal return
        }

        if (i < 4)
        {
            tresh = 0.2 * sm * sm;
            //tresh = 0.2 * sm / (n * n);
        }
        else
        {
            tresh = 0.0;
        }

        for (ip = 0; ip < N - 1; ip++)        
        {
            for (iq = ip + 1; iq < N; iq++)            
            {
                g = 100 * fabs(A[ip][iq]);            // after 4 sweeps, skip the rotation if the off-diagonal element is small
                if ((i > 4) && (fabs(D[ip]) + g == fabs(D[ip])) && (fabs(D[iq]) + g == fabs(D[iq])))
                {
                    A[ip][iq] = 0.0;
                }
                else if (fabs(A[ip][iq]) > tresh)
                {
                    h = D[iq] - D[ip];
                    if (fabs(h) + g == fabs(h))
                    {
                        t = (A[ip][iq]) / h;
                    }
                    else
                    {
                        theta = 0.5 * h / A[ip][iq];
                        t = 1 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0)
                        {
                            t = -t;
                        }
                    }

                    c = 1.0 / sqrt(1.0 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * A[ip][iq];
                    Z[ip] -= h;
                    Z[iq] += h;
                    D[ip] -= h;
                    D[iq] += h;
                    A[ip][iq] = 0.0;
//clock_t t;
    //t=clock();
                    for (j = 0; j < ip; j++)
                    {
                        g = A[j][ip];
                        h = A[j][iq];
                        A[j][ip] = g - s * (h + g * tau);
                        A[j][iq] = h + s * (g - h * tau);
                    }
                    for (j = ip + 1; j < iq; j++)
                    {
                        g = A[ip][j];
                        h = A[j][iq];
                        A[ip][j] = g - s * (h + g * tau);
                        A[j][iq] = h + s * (g - h * tau);
                    }
                    for (j = iq + 1; j < N; j++)
                    {
                        g = A[ip][j];
                        h = A[iq][j];
                        A[ip][j] = g - s * (h + g * tau);
                        A[iq][j] = h + s * (g - h * tau);
                    }
                    for (j = 0; j < N; j++)
                    {
                        g = V[j][ip];
                        h = V[j][iq];
                        V[j][ip] = g - s * (h + g * tau);
                        V[j][iq] = h + s * (g - h * tau);
                    }
//t = clock() - t;
    //printf ("It took me for 4 %d clicks (%f seconds).\n",t,((double)t)/CLOCKS_PER_SEC);
                    *NROT = *NROT + 1;
                }                 //end ((i.gt.4)...else if
            }                     // main iq loop
        }                         // main ip loop

        for (ip = 0; ip < N; ip++)
        {
            B[ip] += Z[ip];
            D[ip] = B[ip];
            Z[ip] = 0.0;
        }
        //cout << "It took me " << endl;
    //cout << i << "\t";
    }                             //end of main i loop
    
    
    

    delete [] B;
    delete [] Z; 
    return;  //too many iterations
}    

int cmp(double a,double b)
{
    return a>b;
}

void Sorting(double **A, int N, double *D, double **V) 
{
      
    double   *D1 = new double [N];              
    int     *nu = new int [N];
    int     i, j;
    double   s;
    double **V1 = (double**)malloc(sizeof(double)*N);
    if(V1)
    {
        for (i = 0; i < n; i++)
        {
            V1[i] = (double*)malloc(sizeof(double)*N);
        }
    }         
    s = D[0];
    
    for (i = 0; i < N; i++)
    {
        D1[i] = D[i];
    }

    sort(D, D + N, cmp);
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (D[i] == D1[j])
            {
                nu[i] = j;
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            V1[j][i] = V[j][nu[i]];
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            V[i][j] = V1[i][j];
        }
    }

    /*int i, j;
    //Sorting eigenvalues and eigenvectors
    map<double, int> mapEigen;
    for (i = 0; i < N; i++) 
    {
        //Eigenvalues[i] = Mat[i][i];
        mapEigen.insert(make_pair(D[i], i));
    }
    double *pdbTmpVec = new double[N * N];
    map<double, int>::reverse_iterator iter = mapEigen.rbegin();
    for (j = 0; iter != mapEigen.rend(), j < N; ++iter, ++j)
    {
        for (i = 0; i < N; i++)
        {
            pdbTmpVec[i*N + j] = V[i][iter->second];
        }
        D[j] = iter->first;
    }
    for (i = 0; i < N; i++)
    {
        double dSumVec = 0;
        for (j = 0; j < N; j++)
        {
            dSumVec += pdbTmpVec[j * N + i];
        }
        if (dSumVec < 0)
        {
            for (j = 0; j < N; j++)
            {
                pdbTmpVec[j * N + i] *= -1;
            }
        }
    }
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            V[i][j] = pdbTmpVec[j * N + i];
        }
    }

    delete [] pdbTmpVec;*/
    delete [] D1;
    delete [] V1;
    return;

}

int main()  {
    
    int nrot;                        // Number of matrix rotations
    int i,j;                         // loop variables
    double vmax;                     // normalization value
    

    //allocate vector and matrices
    double   **Mat;                    // given square symmetric matrix
    double   *Eigenvalues;                // eigenvalues solution vector
    double   **Eigenvectors;           // eigenvectors solution matrix

    Mat = (double**)malloc(sizeof(double*)*n);
    Eigenvalues = new double [n];
    Eigenvectors = (double**)malloc(sizeof(double)*n);
    if(Mat)
    {
        for (i = 0; i < n; i++)
        {
            Mat[i] = (double*)malloc(sizeof(double)*n);
        }
    }
    if(Eigenvectors)
    {
        for (i = 0; i < n; i++)
        {
            Eigenvectors[i] = (double*)malloc(sizeof(double)*n);
        }
    }

    //double Rxx[][4] = {{1.0000,    0.5000,    0.3333,    0.2500},
    //                 {0.5000,    1.0000,    0.6667,    0.5000},
    //                {0.3333,    0.6667,    1.0000,    0.7500},
    //                {0.2500,    0.5000,    0.7500,    1.0000}};
    double Rxx[][16] = {
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

    for (i = 0; i < n; i++) 
    {
        for (j = 0; j < n; j++)
        {
            Mat[i][j] = Rxx[i][j];
        }
    }

    //Displaying the orignal matrix
    cout << "---------------------------------------------------\n" << endl;
    cout <<  "Symmetric matrix is:" << endl;
    for (i = 0; i < n; i++) 
    {
        for (j = 0; j < n; j++)
        {
            cout << Mat[i][j] << "\t";
        }
        cout << "\t" << endl;
    }
    cout << "\t" << endl;

    clock_t t;
    t=clock();
    Jacobi(Mat,n,Eigenvalues,Eigenvectors,&nrot);
    
    
    t = clock() - t;
    printf ("It took me %d clicks (%f seconds).\n",t,((double)t)/CLOCKS_PER_SEC);
    Sorting(Mat,n,Eigenvalues,Eigenvectors); 
    //Displaying the results 
    cout << "Eigenvalues:" << endl;
    for (i = 0; i < n; i++)
    {
        cout << Eigenvalues[i] << "\t";
    }
    cout << "\t" << endl;
    cout << "\t" << endl;

    cout << "Eigenvectors (in lines):" << endl;
    for (j = 0; j < n; j++) 
    {
        vmax = Eigenvectors[1][j];
	    for (i = 0; i < n; i++)
	    if(fabs(Eigenvectors[i][j]) > fabs(vmax))
        {
            vmax = Eigenvectors[i][j];
        }  
        for (i = 0; i < n; i++)
        {
            cout << Eigenvectors[i][j] / vmax << "\t";
        }
	    cout << "\t" << endl;
    }
    cout << "\t" << endl;
  
    cout << "Number of rotations:" << nrot << endl;
    cout << "\t" << endl;
  
    delete [] Mat;
    delete [] Eigenvalues;
    delete [] Eigenvectors; 

}
