#include <iostream>
#include <fstream>
#include <lawrap/blas.h>
#include <vector>

using namespace std;

int nel = 7;

vector<double> F (nel*nel);
vector<double> S (nel*nel);
vector<double> C (nel*nel);
vector<double> H (nel*nel);
vector<double> D (nel*nel);
vector<double> E (nel*nel);
vector<double> DS (nel*nel);

int main()
{
    ifstream F_stream("F.data");
    ifstream C_stream("C.data");
    ifstream S_stream("S.data");
    ifstream H_stream("H.data");
    for(int i = 0; i < nel*nel; i++) 
    {
        F_stream >> F[i];
        C_stream >> C[i];
        S_stream >> S[i];
        H_stream >> H[i];
    }

    // Compute the Density (C*C_T)
    // x.data() gives the pointer to vector x
    LAWrap::gemm('T', 'N', 7, 7, 5, 2.0, C.data(), 7, C.data(), 7, 0.0, D.data(), 7);
    
    // Compute Tr(D*S) to check the density matrix (should be == nocc)
    LAWrap::gemm('N', 'N', 7, 7, 7, 1.0, D.data(), 7, S.data(), 7, 0.0, DS.data(), 7);
    double trace = 0;
    for(int i=0; i< nel*nel; i+=8) trace += DS[i];
    cout << "Trace of DS = " << trace << endl; 

    // Compute F + H
    LAWrap::axpy(nel*nel, 1.0, F.data(), 1, H.data(), 1);

    // Compute the Energy matrix, 0.5*(F+H)*D
    LAWrap::gemm('N', 'N', 7, 7, 7, 0.5, H.data(), 7, D.data(), 7, 0.0, E.data(), 7);

    // Compute the Energy, Tr(E_mat)
    double E_scf = 0;
    for(int i=0; i< nel*nel; i+=8) E_scf += E[i];

    //for(int i=0; i < nel*nel; i++) cout << D[i] << endl; 

    cout << "The SCF energy is: " <<  E_scf << endl;

    F_stream.close();
    C_stream.close();
    S_stream.close();
    F_stream.close();

    return 0;
}
