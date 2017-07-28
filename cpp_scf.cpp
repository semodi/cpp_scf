#include <iostream>
#include <fstream>
#include <lawrap/blas.h>
#include <vector>

using namespace std;

int main()
{
  // Read in input from Python scripts:
  ifstream F_input("F.data");
  vector<double> F(49);
  int i = 0;
  for (i = 0; i <= 48; ++i)
    F_input >> F[i];

  ifstream S_input("S.data");
  vector<double> S(49);
  for (i = 0; i <= 48; ++i)
    S_input >> S[i];

  ifstream H_input("H.data");
  vector<double> H(49);
  for (i = 0; i <= 48; ++i)
    H_input >> H[i];

  ifstream C_input("C.data");
  vector<double> C(49);
  for (i = 0; i <= 48; ++i)
    C_input >> C[i];

  // Calculate D = 2*C*C^T
  // LAWrap::gemm(char transa, char transb, int m, int n, int k,
  //              double alpha, double* A, int lda, double* B,
  //              int ldb, double beta, double* C, int ldc)
  vector<double> D (49);
  LAWrap::gemm('T', 'N', 7, 7, 5, 2.0, C.data(), 7, C.data(), 7, 0.0, D.data(), 7);

  // Energy = Trace[(F+H)*D]/2
  // LAWrap::axpy(int n, double alpha, double* x, int incx,
  //              double* y, int incy)
  LAWrap::axpy(7, 1.0, F.data(), 1, H.data(), 1); // stores F+H in H

  // calculates the inner part which will be given to Trace
  LAWrap::gemm('N', 'N', 7, 7, 7, 1.0, H.data(), 7, D.data(), 7, 0.0, H.data(), 7);


  double trace = 0.0;
  // calculate the Trace
  for (i = 0; i <= 48; i+=8)
    trace += H[i];

  double energy = 0.0;
  energy = trace / 2.0;

  cout << "The energy is: " << energy << '\n';

  // sanity check - Trace[D*S]/2 == ndocc = 5
  LAWrap::gemm('N', 'N', 7, 7, 7, 1.0, D.data(), 7, S.data(), 7, 0.0, S.data(), 7);

  double check = 0.0;
  for (i = 0; i <= 48; i += 8)
    check += S[i];

  check /= 2;

  cout << "Sanity check: Trace[D*S]/2 == ndocc = 5\nTrace[D*S]\2 = " << check << '\n';
  return 0;
}
