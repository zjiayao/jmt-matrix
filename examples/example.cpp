#include <cstdio>
#include "../cmatrix.hpp"
// #include "../dmatrix.hpp"

typedef jmt::mat mat;
typedef jmt::cmat cmat;
using jmt::complex;

int main()
{

	// generate a random matrix
	cmat cm = cmat::getNRand(6, 3);
	cmat U, S, V, Q, R;


	cm.print();
	clock_t t = clock();
	cm.SVD(U, S, V, 1e-10);
	t = clock() - t;
	printf("SVD: %lf sec\n", float(t) / CLOCKS_PER_SEC);
	printf("U:\n");
	U.print();
	printf("S:\n");
	S.print();
	printf("V:\n");
	V.print();

	// check result
	printf("USV^T:\n");
	(U * S * V.hermitian()).print();
	printf("L2 Error: %.3f\n\n", (cm - U * S * V.hermitian()).norm());

	U.QR(Q, R);
	U = Q;
	S = R * S;
	printf("QR decomposition of matrix:\n");
	U.print();

	// multiple printing style
	printf("Q: \n");
	Q.print(jmt::MATLAB);
	printf("R: \n");
	R.print(stdout, jmt::NUMPY);

	return 0;
 }
