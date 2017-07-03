#include <cmath>
#include <cstdio>
#include "matrix.hpp"

#ifndef DMATRIX_HPP
#define DMATRIX_HPP

namespace jmt
{

#ifndef FLOAT_WIDTH
#define FLOAT_WIDTH 8
#endif

#ifndef FLOAT_PRE
#define FLOAT_PRE 3
#endif

#ifndef SEP_SPACE
#define SEP_SPACE 2
#endif


typedef matrix<double> mat;

template <>
void mat::print(FILE *stream, size_type r, size_type c, int precision) const
{
	fprintf(stream, "%*.*f%*c", FLOAT_WIDTH, precision,
							 data[idx(r,c)], SEP_SPACE, ' ');
}
template <>
void mat::read(FILE *stream, size_type i)
{
	fscanf(stream, "%lf", &data[i]);
}

template <>
double mat::sq_norm() const
{
	double res = 0;
	for(size_type i = 0; i < nrow * ncol; ++i)
	{
		res += data[i] * data[i];
	}
	return res;
}

template <>
double mat::sq_norm(size_type r1, size_type c1, size_type r2, size_type c2) const
{
	double res = 0;
	if(bound_check(r1, c1) && bound_check(r2, c2))
	{

		for(size_type i = r1; i <= r2; ++i)
		{
			for(size_type j = c1; j <= c2; ++j)
			{
				res += data[idx(i,j)] * data[idx(i,j)];
			}
		}
	} else { fprintf(STDERR, "sq_norm failed: invalid indices\n");}
	return res;
}

template <>
mat mat::hermitian() const
{
	return transpose();
}
template <>
void mat::QR(mat &Q, mat &R, bool _clear_zero_) const
{
		/* if use Gram-Schmit */
		//this->norm_basis(Q);
		//R = Q.transpose() * (*this);

		/* use householder's reflection */
		Q.eye(nrow, nrow);
		R = this->clone();
		mat Qt = getEye(nrow, nrow);
		size_type n = nrow - 1 < ncol ? nrow - 1 : ncol;

		mat x, vi, Hi;

		for(size_type i = 0; i < n; ++i)
		{
			vi = R.subMat(i, i, nrow - 1, i);
			x = R.subMat(i, i, nrow - 1, i);

			vi(0, 0) -= std::copysign(vi.norm(), vi(0,0));

			mat Hi = getEye(nrow - i, nrow - i) - (2.0 / vi.sq_norm()) * vi * vi.transpose();
			Qt.eye();
			Qt.fill(i, i, nrow - 1, nrow - 1, Hi.data);

			Q = Q * Qt.transpose();
			R = Qt * R;

		}
		if(_clear_zero_)
		{	
			Q.clear_zero();
			R.clear_zero();
		}
}

template <>
void mat::jacobi_svd(const mat &A, mat &U,
								mat &S, mat &V, double criterion)
{
	U = A.clone();
	size_type m = U.nrow, n = U.ncol; // m >= n
	V.eye(n, n);
	mat  Ui(m, 1), Uj(m, 1),
			Vi(n, 1), Vj(n, 1);
			// Ri(n, 1), Rj(n, 1); 
	double conv = criterion + 1.0, conv_tmp;	
	while(conv > criterion)
	{
		conv = 0.0;
		for(size_type j = 1; j < n; ++j)
		{
			for(size_type i = 0; i < j; ++i)
			{
				Ui.fill(0, i, m - 1, i, U);
				Uj.fill(0, j, m - 1, j, U);
				double 	alpha(Ui.sq_norm()),
							beta(Uj.sq_norm()),
							gamma(inner_prod(Ui, Uj));

				if(!feq(gamma, 0.0))
				{					
					conv_tmp = fabs(gamma) / sqrt(alpha * beta);
					if(conv < conv_tmp) { conv = conv_tmp; }

					double	zeta = (beta - alpha) / (2 * gamma) ,
								t = 1.0 / (zeta + std::copysign(1.0, zeta) * sqrt(1.0 + zeta * zeta)),
								c = 1 / sqrt(1 + t * t),
								s = t * c;

					U.fillCol(i, c * Ui - s * Uj);
					U.fillCol(j, s * Ui + c * Uj);
					
					// Ri.fill(i, 0, i, n - 1, U); Rj.fill(j, 0, j, n - 1, U);
					// U.fillRow(i, c * Ri - s * Rj);
					// U.fillRow(j, s * Ri + c * Rj);
					
					Vi.fill(0, i, n - 1, i, V); Vj.fill(0, j, n - 1, j, V);
					V.fillCol(i, c * Vi - s * Vj);
					V.fillCol(j, s * Vi + c * Vj);
				}
			}
		}
	}

	S.zeros(n, n); double tmp_norm = 0.0;
	for(size_type i = 0; i < n; ++i)
	{
		tmp_norm = U.norm(0, i, m - 1, i);
		S.data[S.idx(i,i)] = tmp_norm;
		// printf("Singular value: %lf\n", tmp_norm);
		if(!feq(tmp_norm, 0.0))
		{
			for(size_type j = 0; j < m; ++j)
			{
				U(j, i) /= tmp_norm;
			}				
		}
		// U.getCol(i).print();
		// printf("%u-th norm: %lf - %lf\n", i, U.norm(0, i, n - 1, i), U.getCol(i).norm());

	}

	// selection sort
	for(size_type i = 0; i < n - 1; ++i)
	{
		size_type max = i;
		for(size_type j = i + 1; j < n; ++j)
		{
			if( S(max, max) < S(j, j) ) { max = j; }
		}
		double tmp = S(max, max);
		S(max, max) = S(i ,i);
		S(i, i) = tmp;
		U.exCol(i, max);
		V.exCol(i, max);
	}


	U.clear_zero();
	S.clear_zero();
	V.clear_zero();
} 



}

#endif