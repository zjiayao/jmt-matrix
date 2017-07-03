//////////////////////////////////
//	   Jiayao's Math Toolbox	//
//	  Linear Algebra Library 	//	
// 	 		Matrix Class		//
// 								//
//								//
// 	(C) 2017 Jiayao Zhang		//
//								//
// v0.1 27-Apr-2017				//
//		Basic functionality		//
//		for real matrices.		//
//								//
// v0.2 08-May-2017				//
//		Complex support on		//
//		- Cholesky 				//
//		- Inverse				//
//		- QR					//
//		- SVD					// 
//////////////////////////////////
// matrix.hpp

#include <cmath>
#include <cstdio>
#include <random>
#include "complex.hpp"

#ifndef MATRIX_HPP
#define MATRIX_HPP
#ifndef STDERR
#define STDERR stderr
#endif

#ifndef FLOAT_WIDTH
#define FLOAT_WIDTH 8
#endif

#ifndef FLOAT_PRE
#define FLOAT_PRE 3
#endif

#ifndef SEP_SPACE
#define SEP_SPACE 2
#endif


namespace jmt
{

// typedef double mat_type;
typedef unsigned size_type;

enum mat_style
{
	C, MATLAB, NUMPY, UNFORMATTED
};

// forward declarations
template <class mat_type>
class matrix;


template <typename mat_type>
matrix<mat_type> operator+(const matrix<mat_type> &m1, const matrix<mat_type> &m2);
template <typename mat_type>
matrix<mat_type> operator-(const matrix<mat_type> &m1, const matrix<mat_type> &m2);
template <typename mat_type>
matrix<mat_type> operator*(const matrix<mat_type> &m1, const matrix<mat_type> &m2);
template <typename mat_type>
matrix<mat_type> operator+(const matrix<mat_type> &mat, mat_type p);
template <typename mat_type>
matrix<mat_type> operator+(mat_type p, const matrix<mat_type> &mat);
template <typename mat_type>
matrix<mat_type> operator-(const matrix<mat_type> &mat, mat_type p);
template <typename mat_type>
matrix<mat_type> operator-(mat_type p, const matrix<mat_type> &mat);
template <typename mat_type>
matrix<mat_type> operator*(const matrix<mat_type> &mat, mat_type p);
template <typename mat_type>
matrix<mat_type> operator*(mat_type p, const matrix<mat_type> &mat);
template <typename mat_type>
matrix<mat_type> operator/(const matrix<mat_type> &mat, mat_type p);
template <typename mat_type>
matrix<mat_type> operator/(mat_type p, const matrix<mat_type> &mat);
template <typename mat_type>
mat_type inner_prod(const matrix<mat_type> &m1, const matrix<mat_type> &m2);


template <class mat_type>
class matrix
{
protected:
	mat_type *data;
public:
	size_type nrow, ncol;
private:
	size_type *ref;

protected:
	inline virtual void clean_data()
	{ if(data) { delete [] data; } data = nullptr; }

	inline size_type idx(size_type r, size_type c) const
	{ return r * ncol + c; }


private:

	inline bool bound_check(size_type r, size_type c) const
	{ return (r < nrow) && (c < ncol); }

	/**
	 * recycle
	 *
	 * deallocate dynamic memory,
	 * reset reference counter
	 */
	void recycle()
	{
		if( !ref || !(*ref) )
		{
			this->clean_data();
			if(ref) { delete ref; }
			ref = nullptr;
			ncol = nrow = 0;
		} else if(*ref) {
			--(*ref);
			if(!(*ref))
			{
				this->recycle();
			}
		}
	}

	/**
	 * fill(r1, c1, r2, c2, src)
	 *
	 * fill (r1, c1) to (r2, c2)
	 * in row-major order from array src
	 */
	void fill(size_type r1, size_type c1, size_type r2, size_type c2, mat_type *src)
	{
		size_type k = 0;
		for(size_type i = r1; i <= r2; ++i)
		{
			for(size_type j = c1; j <= c2; ++j)
			{
				data[idx(i, j)] = mat_type(src[k++]);
			}
		}
	}

	/**
	 * fill(idx, n, src)
	 *
	 * fill n elements from the idx-th
	 * entry of data from array src
	 */
	void fill(size_type idx, size_type n, mat_type *src)
	{ 
		for(size_type i = 0; i < n; ++i)
		{
			data[idx + i] = mat_type(src[i]);
		}
	}

	/**
	 * fill(r1, c1, r2, c2, mat)
	 *
	 * fill from begining of this matrix
	 * from (r1, c1) to (r2, c2) of mat
	 * in row-major order from matrix mat
	 */
	void fill(size_type r1, size_type c1, size_type r2, size_type c2, const matrix &mat)
	{
		size_type k = 0;
		for(size_type i = r1; i <= r2; ++i)
		{
			for(size_type j = c1; j <= c2; ++j)
			{
				data[k++] = mat_type(mat.data[mat.idx(i, j)]);
			}
		}
	}

	/**
	 * fill(idx, n, mat)
	 *
	 * fill n elements to idx-th entry of this matrix
	 * from the idx-th entry of mat
	 */
	void fill(size_type idx, size_type n, const matrix &mat)
	{
		for(size_type i = 0; i < n; ++i)
		{
			data[idx + i] = mat_type(mat.data[idx + i]);
		}
	}

	/**
	 * elem_set(e)
	 *
	 * set all element to e
	 */
	void elem_set(const mat_type &e)
	{ 
		for(size_type i = 0; i < nrow * ncol; ++i)
		{ data[i] = mat_type(e); }
	}


public:
	/* constructor, copy, move and destructor */
	// constructor
	matrix() : data(nullptr), nrow(0),  ncol(0), ref(nullptr) {}
	matrix(size_type row, size_type col, mat_type e) : data(new mat_type[row * col]), nrow(row), ncol(col), ref(new size_type(1))
	{ elem_set(e); }
	matrix(size_type row, size_type col) : matrix(row, col, mat_type(0.0)) {}
	matrix(size_type row, size_type col, mat_type *ee) : data(new mat_type[row * col]), nrow(row), ncol(col), ref(new size_type(1))
	{ fill(0, nrow * ncol, ee); }



	// destructor
	~matrix()
	{ this->recycle(); }

	// copy constructor, shallow copy
	matrix(const matrix &mat) : data(mat.data), nrow(mat.nrow), ncol(mat.ncol), ref(mat.ref)
	{ 
		if(!ref)
		{
			fprintf(STDERR, "copy constructor: copy from uninitilized matrix\n");
			return;
		}
		++(*ref); 
	}

	// copy assignment
	virtual matrix& operator=(const matrix &mat) 
	{
		this->recycle();
		this->data = mat.data;
		this->nrow = mat.nrow; this->ncol = mat.ncol;
		this->ref = mat.ref;

		if(ref)
		{ ++(*(this->ref)); }
		else {
			fprintf(STDERR, "copy assignment: assign from uninitilized matrix\n");
		}

		return *this;
	}

	// move construcor, shallow copy
	matrix(matrix &&mat) : data(mat.data), nrow(mat.nrow), ncol(mat.ncol), ref(mat.ref)
	{ 
		mat.data = nullptr;
		mat.recycle();
	}

	// move assignment
	virtual matrix& operator=(matrix &&mat)
	{
		this->recycle();
		this->data = mat.data;
		this->nrow = mat.nrow; this->ncol = mat.ncol;

		if(mat.ref)
		{ this->ref = new size_type(*(mat.ref)); }
		else { this->ref = nullptr; }

		mat.data = nullptr;
		mat.recycle();

		return *this;
	}

	/* matrix access & common matrices */

	/**
	 * clone()
	 *
	 * return a clone of the matrix
	 */
	inline virtual matrix clone() const
	{
		matrix ret(nrow, ncol);
		ret.fill(0, nrow * ncol, *this);
		return ret;
	}

	/**
	 * clone(mat)
	 *
	 * set mat to a clone of the matrix
	 */
	inline matrix clone(matrix &mat) const
	{
		mat = this->clone();
		return mat;
	}

	/**
	 * set(nrow, ncol)
	 *
	 * clean the matrix and reset it with
	 * dimension nrow * ncol
	 */
	virtual void set(size_type nrow, size_type ncol)
	{ 
		this->recycle();
		this->nrow = nrow; this->ncol = ncol; 
		this->ref = new size_type(1);
		this->data = new mat_type[nrow * ncol];			
	}

	/**
	 * fillRow(r, row)
	 *
	 * fill the r-th row from a row matrix
	 * bound checking is performed
	 */
	virtual matrix& fillRow(size_type r, const matrix& row)
	{
		if((	(row.ncol == 1 && row.nrow == this->ncol)
			|| 	(row.nrow == 1 && row.ncol == this->ncol))
			&& bound_check(r, 0))
		{
			this->fill(r * ncol, ncol, row.data);
		}
		return *this;
	}

	/**
	 * fillRow(r, row)
	 *
	 * fill the r-th row from a row array
	 * bound checking is performed
	 */
	virtual matrix& fillRow(size_type r, mat_type *row)
	{
		if(bound_check(r, 0))
		{
			this->fill(r * ncol, ncol, row);
		}
		return *this;
	}

	/**
	 * fillCol(c, col)
	 *
	 * fill the c-th column from a column matrix
	 * bound checking is performed
	 */
	virtual matrix& fillCol(size_type c, const matrix &col)
	{
		if((	(col.ncol == 1 && col.nrow == this->nrow)
			|| 	(col.nrow == 1 && col.ncol == this->nrow))
			&& bound_check(0, c))
		{
			for(size_type i = 0; i < nrow; ++i)
			{
				data[idx(i,c)] = col.data[i];
			}
		}	
		return *this;
	}

	/**
	 * fillCol(c, col)
	 *
	 * fill the c-th column from a column array
	 * bound checking is performed
	 */
	virtual matrix& fillCol(size_type c, mat_type *col)
	{
		if(bound_check(0, c))
		{
			for(size_type i = 0; i < nrow; ++i)
			{
				data[idx(i,c)] = col[i];
			}
		}	
		return *this;
	}

	/**
	 * zeros()
	 *
	 * set a matrix with zeros
	 * matrix should be non-null
	 */
	inline virtual matrix& zeros()
	{ elem_set(mat_type(0)); return *this; }

	/**
	 * zeros(nrow, ncol)
	 *
	 * replace this matrix 
	 * with a nrow * ncol zero matrix
	 */
	inline virtual matrix& zeros(size_type nrow, size_type ncol)
	{
		this->set(nrow, ncol);
		return this->zeros();
	}

	/**
	 * diag(nrow, ncol, dat)
	 *
	 * replace this matrix with a diagonal matrix
	 * filling from array dat
	 */
	virtual matrix& diag(size_type nrow, size_type ncol, mat_type *dat)
	{
		this->zeros(nrow, ncol);
		size_type n = nrow < ncol ? nrow : ncol;
		for(size_type i = 0; i < n; ++i)
		{
			data[idx(i, i)] = dat[i];
		}
		return *this;
	}

	/**
	 * diag(nrow, ncol, mat)
	 *
	 * replace this matrix with a diagonal matrix
	 * filling from matrix mat in row major order
	 */
	inline virtual matrix& diag(size_type nrow, size_type ncol, const matrix &mat)
	{ return this->diag(nrow, ncol, mat.data); }

	/**
	 * ones()
	 *
	 * set a matrix with ones
	 * matrix should be non-null
	 */
	inline virtual matrix& ones()
	{ elem_set(mat_type(1)); return *this; }

	/**
	 * ones(nrow, ncol)
	 *
	 * replace this matrix 
	 * with a nrow * ncol one-matrix
	 */
	inline virtual matrix& ones(size_type nrow, size_type ncol)
	{ 
		this->set(nrow, ncol);
		return this->ones();
	}

	/**
	 * eye()
	 *
	 * set a matrix with identity
	 * matrix should be non-null
	 * but not necessarily be square
	 */
	virtual matrix& eye()
	{
		for(size_type i = 0; i < nrow; ++i)
		{
			for(size_type j = 0; j < ncol; ++j)
			{
				if(i == j)
				{ data[idx(i,j)] = mat_type(1.0); }
				else { data[idx(i,j)] = mat_type(0.0); }
			}
		}	
		return *this;
	}

	/**
	 * eye(nrow, ncol)
	 *
	 * replace this matrix 
	 * with a nrow * ncol identity matrix
	 */
	inline virtual matrix& eye(size_type nrow, size_type ncol)
	{
		this->set(nrow, ncol);
		return this->eye();
	}

	/**
	 * seq(init = 0, step = 1)
	 *
	 * filling the matrix in row major
	 * order, from init with step
	 */
	virtual matrix& seq(mat_type init = mat_type(), mat_type step = double(1))
	{
		for(size_type i = 0; i < nrow * ncol; ++i)
		{
			data[i] = init + step * i;
		}

		return *this;
	}

	/**
	 * seq(nrow, ncol, init = 0, step = 1)
	 *
	 * replace this matrix with one nrow * ncol
	 * matrix then invoke seq(init, step)
	 */
	inline virtual matrix& seq(size_type nrow, size_type ncol, mat_type init = mat_type(), mat_type step = double(1))
	{
		this->set(nrow, ncol);
		return this->seq(init, step);
	}

	virtual matrix& nrand(double mu = 0.0, double sigma = 1.0)
	{
		std::default_random_engine gen;
		std::normal_distribution<double> normal(mu, sigma);
		for(size_type i = 0; i < ncol * nrow; ++i)
		{
			data[i] = mat_type(normal(gen));
		}
		return *this;
	}

	inline virtual matrix& nrand(size_type nrow, size_type ncol, double mu = 0.0, double sigma = 1.0)
	{
		this->set(nrow, ncol);
		return this->nrand();
	}

	inline static matrix getNRand(size_type nrow, size_type ncol, double mu = 0.0, double sigma = 1.0)
	{
		matrix ret(nrow, ncol);
		return ret.nrand(mu, sigma);
	}

	inline static matrix getSeq(size_type nrow, size_type ncol, mat_type init = mat_type())
	{
		matrix ret(nrow, ncol);
		return ret.seq(init);
	}

	inline static matrix getOnes(size_type nrow, size_type ncol)
	{
		matrix ret(nrow, ncol);
		return ret.ones();
	}

	inline static matrix getZeros(size_type nrow, size_type ncol)
	{
		matrix ret(nrow, ncol);
		return ret.zeros();
	}

	inline static matrix getEye(size_type nrow, size_type ncol)
	{
		matrix ret(nrow, ncol);
		ret.eye();
		return ret;
	}



	// IOs
// protected:
	// to be overriden / specialized
	virtual void read(FILE *stream, size_type i) {}
	
public:
	void read(FILE *stream = stdin)
	{
		if(ref && *ref)
		{
			for(size_type i = 0; i < nrow * ncol; ++i)
			{
				this->read(stream, i);
			}

		} else {
			if(!ref) { ref = new size_type; }
			*ref = 1;
			if(!(ncol && nrow))
			{ fscanf(stream, "%u%u", &nrow, &ncol); }
			this->read(stream);
		}
	}

// protected:
	// to be overriden / specialized
	virtual void print(FILE *stream, size_type r, size_type c, int precision) const {}
public:
	void print(FILE *stream, mat_style style, int precision = FLOAT_PRE) const
	{
		char lb = '[', rb = ']', sep = ' ', emp = ' ';
		if(style == C)
		{
			lb = '{'; rb = '}'; sep = ',';
		} else if(style == NUMPY) {
			lb = '['; rb = ']'; sep = ',';
		} else if(style == UNFORMATTED) {
			lb = rb = sep = ' ';
		}

		if(ref && (*ref))
		{
			for(size_type i = 0; i < nrow; ++i)
			{
				fprintf(stream, "%c%c ", (i == 0 ? lb : emp), lb);
				for(size_type j = 0; j < ncol; ++j)
				{
					print(stream, i, j, precision);
					fprintf(stream, "%c", (j == ncol - 1 ? emp : sep) );
				}
				fprintf(stream, " %c%c\n", rb, i == nrow - 1 ? rb : sep);
			}
		} else {
			fprintf(STDERR, "matrix not printable\n");
		}
	}

	inline virtual void print(mat_style style, int precision = FLOAT_PRE) const
	{ this->print(stdout, style, precision); }

	inline virtual void print(FILE *stream) const
	{ this->print(stream, MATLAB); }

	inline virtual void print() const
	{ this->print(stdout, MATLAB); }

	/* element access */
	virtual const mat_type& operator()(size_type r, size_type c) const
	{ return data[idx(r, c)]; }
	virtual mat_type& operator()(size_type r, size_type c)
	{ return data[idx(r, c)]; }
	virtual const mat_type* operator[](size_type r) const
	{ return &data[r * ncol]; }
	virtual mat_type* operator[](size_type r)
	{ return &data[r * ncol]; }

	/* submatrix / row / col */
	virtual matrix& exRow(size_type i, size_type j)
	{
		mat_type tmp;
		for(size_type n = 0; n < ncol; ++n)
		{
			tmp = data[idx(i ,n)];
			data[idx(i, n)] = data[idx(j, n)];
			data[idx(j ,n)] = tmp;
		}
		return *this;
	}

	virtual matrix& exCol(size_type i, size_type j)
	{
		mat_type tmp;
		for(size_type n = 0; n < nrow; ++n)
		{
			tmp = data[idx(n ,i)];
			data[idx(n, i)] = data[idx(n, j)];
			data[idx(n ,j)] = tmp;
		}
		return *this;
	}

	virtual matrix& operator*=(double d)
	{ 
		for(size_type i = 0; i < nrow * ncol; ++i)
		{
			data[i] *= d;
		}

		return *this;
	}
	// return new matrix 
	virtual matrix getCol(size_type col) const
	{
		matrix res(nrow, 1);
		if(bound_check(0, col))
		{
			for(size_type i = 0; i < nrow; ++i)
			{
				res(i, 0) = data[idx(i, col)];
			}

		} else {
			fscanf(STDERR, "failed obtain column: index out of bound\n");
		}
		return res;
	}

	// return new matrix
	virtual matrix getRow(size_type r) const
	{
		matrix res(1, ncol);
		if(bound_check(r, 0))
		{
			res.fill(0, ncol, (data + r * ncol));
		} else {
			fscanf(STDERR, "failed obtain row: index out of bound\n");
		}
		return res;
	}

	// return new matrix
	virtual matrix subMat(size_type r1, size_type c1, size_type r2, size_type c2) const
	{
		matrix res;
		if(bound_check(r1, c1) && bound_check(r2, c2) && (r1 <= r2) && (c1 <= c2))
		{
			res.set(r2 - r1 + 1, c2 - c1 + 1);
			size_type k = 0;
			for(size_type i = r1; i <= r2; ++i)
			{
				for(size_type j = c1; j <= c2; ++j)
				{
					res.data[k++] = data[idx(i, j)];
				}
			}
		} else {
			fscanf(STDERR, "failed obtain submatrix: invalid indices\n");
		}
		return res;
	}

	/* arithmatics */
	// to be specialized
	virtual double sq_norm(size_type r1, size_type c1, size_type r2, size_type c2) const
	{}
	// to be specialized
	virtual double sq_norm() const
	{}

	virtual double norm(size_type r1, size_type c1, size_type r2, size_type c2) const
	{ return std::sqrt(this->sq_norm(r1, c1, r2, c2)); }

	virtual double norm() const
	{ return std::sqrt(this->sq_norm()); }

	virtual matrix& add(const matrix &mat)
	{
		if(this->nrow == mat.nrow && this->ncol == mat.ncol)
		{
			for(size_type i = 0; i < nrow; ++i)
			{
				for(size_type j = 0; j < ncol; ++j)
				{
					data[idx(i,j)] += mat(i,j);
				}
			}
		} else { fprintf(STDERR, "operation add failed, dimension not match\n"); }
		return *this;
	}

	virtual matrix& operator+=(const matrix &mat)
	{ return this->add(mat); }

	virtual matrix& sub(const matrix &mat)
	{
		if(this->nrow == mat.nrow && this->ncol == mat.ncol)
		{
			for(size_type i = 0; i < nrow; ++i)
			{
				for(size_type j = 0; j < ncol; ++j)
				{
					data[idx(i,j)] -= mat(i,j);
				}
			}
		} else { fprintf(STDERR, "operation sub failed, dimension not match\n"); }
		return *this;
	}

	virtual matrix& operator-=(const matrix &mat)
	{ return this->sub(mat); }


	/* transformations */
	virtual matrix transpose() const
	{
		matrix res(this->ncol, this->nrow);
		for(size_type i = 0; i < this->nrow; ++i)
		{
			for(size_type j = 0; j < this->ncol; ++j)
			{
				res(j,i) = data[idx(i,j)];
			}
		}
		return res;
	}

	// to be specialized for Complex
	virtual matrix hermitian() const
	{ }
// protected:
	virtual matrix& clear_zero()
	{
		for(size_type i = 0; i < nrow * ncol; ++i)
		{
			if(feq(data[i], mat_type(0.0))) { data[i] = mat_type(0.0); }
		}
		return *this;
	}
public:
	virtual matrix& normalizeCol()
	{
		for(size_type i = 0; i < ncol; ++i)
		{
			matrix tmp = this->getCol(i);
			double coef = tmp.norm();
			for(size_type j = 0; j < nrow; ++j)
			{
				if(!feq(coef, 0.0) && !feq(data[idx(i,j)], mat_type(0.0)))
				{ data[idx(j,i)] /= coef; }
				else { data[idx(j,i)] = 0.0; }
			}
		}
		return *this;
	}

	inline virtual matrix& normalize()
	{
		double coef = this->norm();
		for(size_type i = 0; i < nrow * ncol; ++i)
		{
			data[i] /= coef;
		}
		return *this;
	}

	inline virtual mat_type trace() const
	{
		size_type n = fmin(nrow, ncol);
		mat_type res(0);
		for(size_type i = 0; i < n; ++i)
		{
			res += data[idx(i,i)];
		}
		return res;
	}

	// use Cholesky decomposition
	inline virtual matrix inverse() const
	{
		matrix 	res = this->cholesky(), 
				inv = res.linv();
		res = inv.hermitian() * inv;
		res.clear_zero();
		return res;
	}


	// use Gram-Schmit Process
	virtual void basis(matrix &B) const
	{
		const size_type &r = this->nrow, &c = this->nrow;
		B.set(r, c);
		for(size_type i = 0; i < c; ++i)
		{
			matrix vi = this->getCol(i),
				sum(r, 1);
			for(size_type j = 0; j < i; ++j)
			{
				matrix tmp = B.getCol(j);
				if(!feq(tmp.sq_norm(), 0))
				{ sum += (inner_prod(vi, tmp) / tmp.sq_norm()) * tmp; }
				else { sum = vi; break; }
			}
			vi -= sum;

			B.fillCol(i, vi);
		}
	}

	virtual void norm_basis(matrix &B) const
	{
		this->basis(B);
		B.normalizeCol();
	}

	// use Householder Relection
	virtual void QR(matrix &Q, matrix &R, bool _clear_zero_ = true) const
	{ /* specialized */}

	// m * n matrix =>
	// U: m * n orthogonal
	// S: n * n diagonal
	// V: n * n orthogonal
	virtual void SVD(matrix &U, matrix &S, matrix &V, double criterion=1e-8) const
	{
		bool trans = false;
		matrix res, q, r;
		if(nrow < ncol) { trans = true; res = this->hermitian(); }
		else { res = this->clone(); }

		jacobi_svd(res, U, S, V, criterion);
		
		if(trans)
		{
			res = U;
			U = V;
			V = res;
		}
	}

	virtual void sym_svd(matrix &U, matrix &S, double criterion=1e-8)
	{ /* specialized */ }


	virtual matrix cholesky() const
	{
		matrix L;
		if(nrow != ncol)
		{
			fprintf(STDERR, "cholesky failed: matrix not square\n");
		} else { matrix::cholesky(*this, L); }
		return L;
	}

	virtual matrix linv() const
	{
		matrix res;
		matrix::linv(*this, res);
		return res;
	}

private:
	// inverse of a lower triangular matrix L
	// S = L^-1
	static void linv(const matrix &L, matrix &S)
	{
		S.zeros(L.nrow, L.ncol);
		mat_type tmp(0);
		for(size_type j = 0; j < L.ncol; ++j)
		{
			S(j, j) = 1.0 / L(j, j);

			for(size_type i = j + 1; i < L.nrow; ++i)
			{
				// compute S(i, j)
				tmp = mat_type(0);
				for(size_type k = j; k < i; ++k)
				{
					tmp += L(i, k) * S(k,j);
				}
				S(i, j) = (- tmp) / L(i, i);
			}
		}

	}

	// single-sided Jacobi SVD algorithm
	// for m * n matrix with m >= n
	static void jacobi_svd(const matrix &A, matrix &U,
									matrix &S, matrix &V, double criterion)
	{ /* specialized */ } 

	// A is hermitian
	// L is lower triangular
	// A = L * (L*)
	static void cholesky(const matrix &A, matrix &L)
	{
		size_type n = A.nrow;
		L.zeros(n ,n);
		mat_type tmp(0);
		double res = 0.0;
		for(size_type i = 0; i < n; ++i)
		{
			for(size_type j = 0; j < i; ++j)
			{
				tmp = mat_type(0.0);
				// compute L(i, j)
				for(size_type k = 0; k < j; ++k)
				{
					tmp += L(i, k) * conjugate(L(j, k));
				}
				L(i, j) = (A(i, j) - tmp) / L(j, j);

				// printf("L(%u,%u) = %lf\n", i, j, L(i, j));
			}

			// compute L(i, i)
			res = 0.0;
			for(size_type k = 0; k < i; ++k)
			{
				res += jmt::sq_norm(L(i, k));
			}

			res = real(A(i, i)) - res;

			if(res < 0 || !feq(im(A(i,i)), 0.0) )
			{
				fprintf(STDERR, "cholesky failed: matrix must be hermitian and positive definite\n");
				L.zeros();
				return;
			}

			L(i , i) = mat_type(::sqrt(res));

			// printf("L(%u,%u) = %lf\n", i, i, L(i, i));

		}
	}


	/* friend operators */
	friend matrix operator+<mat_type>(const matrix &m1, const matrix &m2);
	friend matrix operator-<mat_type>(const matrix &m1, const matrix &m2);
	friend matrix operator*<mat_type>(const matrix &m1, const matrix &m2);
	friend matrix operator+<mat_type>(const matrix &mat, mat_type p);
	friend matrix operator+<mat_type>(mat_type p, const matrix &mat);
	friend matrix operator-<mat_type>(const matrix &mat, mat_type p);
	friend matrix operator*<mat_type>(const matrix &mat, mat_type p);
	friend matrix operator*<mat_type>(mat_type p, const matrix &mat);
	friend matrix operator/<mat_type>(const matrix &mat, mat_type p);
	friend matrix operator/<mat_type>(mat_type p, const matrix &mat);
	friend matrix operator/<mat_type>(mat_type p, const matrix &mat);
	friend mat_type inner_prod<mat_type>(const matrix &m1, const matrix &m2);
};

// operators
template <typename mat_type>
matrix<mat_type> operator+(const matrix<mat_type> &m1, const matrix<mat_type> &m2)
{
	matrix<mat_type> res;
	if(m1.nrow == m2.nrow && m1.ncol == m2.ncol)
	{
		res.set(m1.nrow, m1.ncol);
		for(size_type i = 0; i < m1.nrow; ++i)
		{
			for(size_type j = 0; j < m1.ncol; ++j)
			{
				res(i,j) = m1(i, j) + m2(i,j);
			}
		}
	} else { fprintf(STDERR, "operation+ failed, dimension not match\n"); }
	return res;
}

template <typename mat_type>
matrix<mat_type> operator-(const matrix<mat_type> &m1, const matrix<mat_type> &m2)
{
	matrix<mat_type> res;
	if(m1.nrow == m2.nrow && m1.ncol == m2.ncol)
	{
		res.set(m1.nrow, m1.ncol);
		for(size_type i = 0; i < m1.nrow; ++i)
		{
			for(size_type j = 0; j < m1.ncol; ++j)
			{
				res(i,j) = m1(i, j) - m2(i,j);
			}
		}
	} else { fprintf(STDERR, "operation- failed, dimension not match\n"); }
	return res;
}

template <typename mat_type>
matrix<mat_type> operator*(const matrix<mat_type> &m1, const matrix<mat_type> &m2)
{
	matrix<mat_type> res;
	if(m1.ncol == m2.nrow)
	{
		res.zeros(m1.nrow, m2.ncol);

		for(size_type k = 0; k < m1.ncol; ++k)
		{
			for(size_type i = 0; i < m1.nrow; ++i)
			{
				const mat_type &r = m1(i, k);
				for(size_type j = 0; j < m2.ncol; ++j)
				{
					res(i, j) += r * m2(k, j);
				}
			}
		}
		/*
		for(size_type i = 0; i < m1.nrow; ++i)
		{
			for(size_type j = 0; j < m2.ncol; ++j)
			{
				for(size_type k = 0; k < m1.ncol; ++k)
				{
					res(i,j) += m1(i,k) * m2(k,j);
				}
			}
		}
		*/
	} else { fprintf(STDERR, "operation* failed, dimension not match\n"); }
	return res;
}

template <typename mat_type>
matrix<mat_type> operator+(const matrix<mat_type> &mat, mat_type p)
{
	matrix<mat_type> res;
	res.set(mat.nrow, mat.ncol);
	for(size_type i = 0; i < mat.nrow; ++i)
	{
		for(size_type j = 0; j < mat.ncol; ++j)
		{
			res(i,j) = mat(i,j) + p;
		}
	}
	return res;
}

template <typename mat_type>
matrix<mat_type> operator+(mat_type p, const matrix<mat_type> &mat)
{ return mat + p; }

template <typename mat_type>
matrix<mat_type> operator-(const matrix<mat_type> &mat, mat_type p)
{
	matrix<mat_type> res;
	res.set(mat.nrow, mat.ncol);
	for(size_type i = 0; i < mat.nrow; ++i)
	{
		for(size_type j = 0; j < mat.ncol; ++j)
		{
			res(i,j) = mat(i,j) - p;
		}
	}
	return res;
}

template <typename mat_type>
matrix<mat_type> operator-(mat_type p, const matrix<mat_type> &mat)
{ return mat - p; }


template <typename mat_type>
matrix<mat_type> operator*(const matrix<mat_type> &mat, mat_type p)
{
	matrix<mat_type> res;
	res.set(mat.nrow, mat.ncol);
	for(size_type i = 0; i < mat.nrow; ++i)
	{
		for(size_type j = 0; j < mat.ncol; ++j)
		{
			res(i,j) = mat(i,j) * p;
		}
	}
	return res;
}

template <typename mat_type>
matrix<mat_type> operator*(mat_type p, const matrix<mat_type> &mat)
{ return mat * p; }

template <typename mat_type>
matrix<mat_type> operator/(const matrix<mat_type> &mat, mat_type p)
{
	matrix<mat_type> res;
	res.set(mat.nrow, mat.ncol);
	for(size_type i = 0; i < mat.nrow; ++i)
	{
		for(size_type j = 0; j < mat.ncol; ++j)
		{
			res(i,j) = mat(i,j) / p;
		}
	}
	return res;
}

template <typename mat_type>
matrix<mat_type> operator/(mat_type p, const matrix<mat_type> &mat)
{ return mat / p; }

template <typename mat_type>
mat_type inner_prod(const matrix<mat_type> &m1, const matrix<mat_type> &m2)
{
	mat_type res(0.0);
	if( (m1.nrow == m2.nrow) && (m1.ncol == 1) && (m2.ncol == 1))
	{
		for(size_type i = 0; i < m1.nrow; ++i)
		{
			res += m1(i, 0) * m2(i, 0);
		}

	} else {
		fscanf(STDERR, "inner_prod only valid for column vectors\n");
	}
	return res;
}

} // namespace

#endif
