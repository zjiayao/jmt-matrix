#include <cmath>
#include <cstdio>
#include <ctime>
#include "complex.hpp"
#include "matrix.hpp"
#include "dmatrix.hpp"

#ifndef CMATRIX_HPP
#define CMATRIX_HPP

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

typedef matrix<complex> cmat;


template <>
void cmat::print(FILE *stream, size_type r, size_type c, int precision) const
{
	const complex &p = data[idx(r,c)];
	fprintf(stream, "%-*.*f", FLOAT_WIDTH, precision, p.re);
	if(feq(p.im, 0.0))
	{
		fprintf(stream, "  %*c", FLOAT_WIDTH + SEP_SPACE + 2, ' ');
	} else { fprintf(stream, " %c %*.*fi%*c", p.im > 0 ? '+' : '-', FLOAT_WIDTH
		, precision, fabs(p.im), SEP_SPACE, ' ');}
}
template <>
void cmat::read(FILE *stream, size_type i)
{
	char sign;
	fscanf(stream, "%lf %c %lf%*c", &(data[i].re), &sign, &(data[i].im));
	if(sign == '-') { data[i].im = -data[i].im; }
}

template <>
double cmat::sq_norm() const
{
	double res = 0;
	for(size_type i = 0; i < nrow * ncol; ++i)
	{
		res += data[i].sq_norm();
	}
	return res;
}

template <>
cmat& cmat::nrand(double mu, double sigma)
{
	unsigned seed = std::clock();
	std::default_random_engine gen ( seed );
	std::normal_distribution<double> normal(mu, sigma);
	for(size_type i = 0; i < ncol * nrow; ++i)
	{
		data[i] = complex(normal(gen), normal(gen));
	}
	return *this;
}



template <>
double cmat::sq_norm(size_type r1, size_type c1, size_type r2, size_type c2) const
{
	double res = 0;
	if(bound_check(r1, c1) && bound_check(r2, c2))
	{

		for(size_type i = r1; i <= r2; ++i)
		{
			for(size_type j = c1; j <= c2; ++j)
			{
				const complex &c = data[idx(i,j)];
				res += jmt::sq_norm(c);
			}
		}
	} else { fprintf(STDERR, "sq_norm failed: invalid indices\n");}
	return res;
}

template <>
complex inner_prod(const cmat &m1, const cmat &m2)
{
	complex res(0.0, 0.0);
	if( (m1.nrow == m2.nrow) && (m1.ncol == 1) && (m2.ncol == 1))
	{
		for(size_type i = 0; i < m1.nrow; ++i)
		{
			res += conjugate(m1(i, 0)) * m2(i, 0);
		}

	} else {
		fscanf(STDERR, "inner_prod only valid for column vectors\n");
	}
	return res;
}

template <>
cmat cmat::hermitian() const
{
	matrix res(this->ncol, this->nrow);
	for(size_type i = 0; i < this->nrow; ++i)
	{
		for(size_type j = 0; j < this->ncol; ++j)
		{
			res(j,i) = conjugate(data[idx(i,j)]);
		}
	}
	return res;
}
template <>
cmat& cmat::clear_zero()
{
	for(size_type i = 0; i < nrow * ncol; ++i)
	{
		if(feq(real(data[i]), 0.0)) { data[i].re = 0.0; }
		if(feq(im(data[i]), 0.0)) { data[i].im = 0.0; }

	}
	return *this;
}

template <>
void cmat::QR(cmat &Q, cmat &R, bool _clear_zero_) const
	{
		/* if use Gram-Schmit */
		//this->norm_basis(Q);
		//R = Q.transpose() * (*this);

		/* use householder's reflection */
		Q.eye(nrow, nrow);
		R = this->clone();
		cmat Qt = getEye(nrow, nrow);
		size_type n = nrow - 1 < ncol ? nrow - 1: ncol;

		cmat x, vi, Hi;

		for(size_type i = 0; i < n; ++i)
		{
			vi = R.subMat(i, i, nrow - 1, i);
			x = R.subMat(i, i, nrow - 1, i);
			complex sign = jmt::normalize(vi(0, 0));
			vi(0, 0) += vi.norm() * sign;// vi(0,0) > 0 ? -vi.norm() : vi.norm();

			cmat Hi = (getEye(nrow - i, nrow - i) - complex(2.0 / vi.sq_norm()) * vi * vi.hermitian());
			Qt.eye();
			Qt.fill(i, i, nrow - 1, nrow - 1, Hi.data);
			Q = Q * Qt.hermitian();
			R = Qt * R;

		}
		if(_clear_zero_)
		{	
			Q.clear_zero();
			R.clear_zero();
		}
	}


template <>
void cmat::jacobi_svd(const cmat &A, cmat &U,
					cmat &S, cmat &V, double criterion)
{
	// m >= n
	size_type m = A.nrow, n = A.ncol;
	U.zeros(m, n); S.eye(n, n); V.eye(n, n);

	mat R(m, n), I(m, n), Ii(m, n);
	for(size_type i = 0; i < m * n; ++i)
	{
		R.data[i] = A.data[i].re;
		I.data[i] = A.data[i].im;
		Ii.data[i] = -A.data[i].im;
	}

	mat K(2 * m, 2 * n);
	K.fill(0, 0, m - 1, n - 1, R.data);
	K.fill(0, n, m - 1, 2 * n - 1, Ii.data);
	K.fill(m, 0, 2 * m - 1, n - 1, I.data);
	K.fill(m, n, 2 * m - 1, 2 * n - 1, R.data);
	mat u, s, v;
	// K.print();
	K.SVD(u, s, v, criterion);
	// u.print();s.print();
	for(size_type i = 0; i < m; ++i)
	{
		for(size_type j = 0; j < n; ++j)
		{
			U(i, j).re = u(i, 2 * j);
			U(i, j).im = u(i + m, 2 * j);
			if(i < n)
			{
				V(i, j).re = v(i, 2 * j);
				V(i, j).im = v(i + n, 2 * j);
			}
			// printf("%u:%u\n", i, j);
		}
		if(i < n)
		{ S(i, i) = s(2 * i, 2 * i); }
	}

	/*

	U = A.clone();
	size_type m = U.nrow, n = U.ncol; // m >= n
	V.eye(n, n);
	cmat  Ui(m, 1), Uj(m, 1),
			Vi(n, 1), Vj(n, 1);

	double conv = criterion + 1.0, conv_tmp;
	unsigned counter = 0;	

	while(conv > criterion)
	{
		conv = 0.0;
		for(size_type j = 1; j < n; ++j)
		{
			for(size_type i = 0; i < j; ++i)
			{
				Ui.fill(0, i, m - 1, i, U);
				Uj.fill(0, j, m - 1, j, U);
				complex 	alpha(Ui.sq_norm()),
							beta(Uj.sq_norm()),
							gamma(inner_prod(Ui, Uj));

				// if(!feq(gamma.re, 0.0) && !feq(alpha.re - beta.re, 0.0))
				// {					
					conv_tmp = jmt::norm(gamma) / jmt::norm(alpha * beta);
					if(conv < conv_tmp) { conv = conv_tmp; }


					// cubersome to work out math
					// but faster indeed
					
					 // phi_2, theta_2 
					double 		sp2, 
								cp2,
								denom = real(alpha - beta);

					if(feq(denom, 0.0)) // pi / 2
					{
						sp2 = 1.0; cp2 = 0.0;

						if(denom < 0.0)
						{ sp2 = -1; }
					} else {
						double 	tg_phi_2 = 2.0 * jmt::norm(gamma) / denom;
								cp2 = 1.0 / sqrt(1.0 + tg_phi_2 * tg_phi_2);
								sp2 = cp2 * tg_phi_2;

					}

					 // phi_1, theta_1 
					double		sp1,
								cp1;
								denom = gamma.re;

					if(feq(denom, 0.0))
					{
						sp1 = -1.0; cp1 = 0.0;
						if((denom >= 0.0 && gamma.im >= 0.0) || (denom < 0.0 && gamma.im < 0.0))
						{
							sp1 = 1.0;

						}
					} else {
						double  tg_phi_1 = gamma.im / denom;
								cp1 = 1.0 / sqrt(1.0 + tg_phi_1 * tg_phi_1);
								sp1 = cp1 * tg_phi_1;
					}
									

					double	csum = cp1 * cp2 - sp1 * sp2,
							ssum = cp1 * sp2 + cp2 * sp1,
							thsum = (1 - csum) / ssum,
							cdif = cp1 * cp2 + sp1 * sp2,
							sdif = sp1 * cp2 - sp2 * cp1,
							thdif = (1 - cdif) / sdif,
							tsum = 1.0 - 2.0 / (1.0 + thsum),
							tdif = 1.0 - 2.0 / (1.0 + thdif),
							ctdif = 1.0 / sqrt(1.0 + tdif * tdif),
							ctsum = 1.0 / sqrt(1.0 + tsum * tsum),
							B = 0.5 * (ctdif + ctsum),
							A = 0.5 * (ctsum - ctdif),
							stdif = ctdif * tdif,
							stsum = ctsum * tsum,
							C = 0.5 * (stsum - stdif),
							D = 0.5 * (stsum + stdif);

							complex s1(-A, -C),
									c1(-B, -D),
									c2(B, -D),
									s2(-A, C);

					// Givens rotaion matrix:
						// s1 c1
						// c2 s2
								
					U.fillCol(i, s1 * (Ui) + c2 * (Uj));
					U.fillCol(j, c1 * (Ui) + s2 * (Uj));

					Vi.fill(0, i, n - 1, i, V); Vj.fill(0, j, n - 1, j, V);
					V.fillCol(i, s1 * Vi + c2 * Vj);
					V.fillCol(j, c1 * Vi + s2 * Vj);

				// }

			}
		}
							printf("%u-th iteration: %lf\n", ++counter, conv);
							// U.print();
							(U.hermitian() * U ).print();

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

	U.clear_zero();
	S.clear_zero();
	V.clear_zero();*/
} 

// a one week's time on calibrating a
// complex SVD is not even as a fraction of 
// accuracy of simply doing svd on equivalent
// real matrix.
// let me cry for the humanity.
template <>
void cmat::sym_svd(cmat &U, cmat &S, double criterion)
{
	mat R(nrow, ncol), I(nrow, ncol), Ii(nrow, ncol);
	for(size_type i = 0; i < nrow * ncol; ++i)
	{
		R.data[i] = this->data[i].re;
		I.data[i] = this->data[i].im;
		Ii.data[i] = -this->data[i].im;
	}
	mat K(2 * nrow, 2 * ncol);
	K.fill(0, 0, nrow - 1, ncol - 1, R.data);
	K.fill(0, ncol, nrow - 1, 2 * ncol - 1, Ii.data);
	K.fill(nrow, 0, 2 * nrow - 1, ncol - 1, I.data);
	K.fill(nrow, ncol, 2 * nrow - 1, 2 * ncol - 1, R.data);
	mat u, s, v;
	K.SVD(u, s, v, criterion);
	U.zeros(nrow, ncol); S.eye(nrow, ncol);


	for(size_type i = 0; i < nrow; ++i)
	{
		for(size_type j = 0; j < ncol; ++j)
		{
			U(i, j).re = u(i, 2 * j);
			U(i, j).im = u(i + nrow, 2 * j);
			// printf("%u:%u\n", i, j);
		}
		S(i, i) = s(2 * i, 2 * i);
	}

	/*
	S = this->clone();
	size_type m = S.nrow, n = S.ncol; // m >= n
	U.eye(n, n);
	cmat  	Ui(m, 1), Uj(m, 1),
			Vi(1, n), Vj(1, n);
			unsigned counter = 0;

	double conv = criterion + 1.0, conv_tmp;
	while(conv > criterion)
	{
		conv = 0.0;
		for(size_type j = 1; j < n; ++j)
		{
			for(size_type i = 0; i < j; ++i)
			{
				complex 	alpha(S(i, i)),
							beta(S(j, j)),
							gamma(S(i, j));
				printf("%u-th itr: %lf\n", counter++, conv);
				// if(!feq(gamma, complex(0.0)))
				// {					
					conv_tmp = fabs(gamma) / jmt::norm(alpha * beta);
					if(conv < conv_tmp) { conv = conv_tmp; }

							// cubersome to work out math
							// but faster indeed
					double  tg_phi_1 = gamma.im / gamma.re,
							tg_phi_2 = 2 * jmt::norm(gamma) / real(alpha - beta),
							cp1 = 1.0 / sqrt(1 + tg_phi_1 * tg_phi_1),
							cp2 = 1.0 / sqrt(1 + tg_phi_2 * tg_phi_2),
							sp1 = cp1 * tg_phi_1,
							sp2 = cp2 * tg_phi_2,
							csum = cp1 * cp2 - sp1 * sp2,
							ssum = cp1 * sp2 + cp2 * sp1,
							thsum = (1 - csum) / ssum,
							cdif = cp1 * cp2 + sp1 * sp2,
							sdif = sp1 * cp2 - sp2 * cp1,
							thdif = (1 - cdif) / sdif,
							tsum = 1.0 - 2.0 / (1.0 + thsum),
							tdif = 1.0 - 2.0 / (1.0 + thdif),
							ctdif = 1.0 / sqrt(1.0 + tdif * tdif),
							ctsum = 1.0 / sqrt(1.0 + tsum * tsum),
							B = 0.5 * (ctdif + ctsum),
							A = 0.5 * (ctsum - ctdif),
							stdif = ctdif * tdif,
							stsum = ctsum * tsum,
							C = 0.5 * (stsum - stdif),
							D = 0.5 * (stsum + stdif);

							complex s1(-A, -C),
									c1(-B, -D),
									c2(B, -D),
									s2(-A, C);


					// Givens rotaion matrix:
						// s1 c1
						// c2 s2

						// s1  ~c2
						// ~c1 s2
					Ui.fill(0, i, m - 1, i, S);
					Uj.fill(0, j, m - 1, j, S);								
					S.fillCol(i, s1 * (Ui) + c2 * (Uj));
					S.fillCol(j, c1 * (Ui) + s2 * (Uj));

					Ui.fill(0, i, m - 1, i, U);
					Uj.fill(0, j, m - 1, j, U);								
					U.fillCol(i, s1 * (Ui) + c2 * (Uj));
					U.fillCol(j, c1 * (Ui) + s2 * (Uj));

					Vi.fill(i, 0, i, n - 1, S);
					Vj.fill(j, 0, j, n - 1, S);
					S.fillRow(i, s1.conj() * Vi + c2.conj() * Vj);
					S.fillRow(j, c1.conj() * Vi + s2.conj() * Vj);
				// }
			}
		}
	}
	
	U.clear_zero();
	S.clear_zero();*/
}

} // namespace

#endif