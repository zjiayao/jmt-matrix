#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#ifndef FERR
#define FERR 1e-10
#endif

namespace jmt
{

#ifndef  M_PI
#define M_PI 3.141592653589793238L
#endif
struct complex
{
	double re, im;
	complex(double re, double im) : re(re), im(im) {}
	// not allow implicit conversion
	explicit complex(double re) : complex(re, 0.0) {}
	complex() : complex(0.0, 0.0) {}
	complex(const complex &c) : re(c.re), im(c.im) {}
	complex& operator=(const complex &c)
	{ re = c.re; im = c.im; return *this; }

	complex& operator=(double d)
	{ re = d; im = 0.0; return *this; }

	inline complex operator+=(const complex &c)
	{ re += c.re; im += c.im; return *this; }
	inline complex operator+=(double d)
	{ re += d; return *this; }
	inline complex operator-=(const complex &c)
	{ re -= c.re; im -= c.im; return *this; }
	inline complex operator-=(double d)
	{ re -= d; return *this; }
	inline complex operator*=(const complex &c)
	{ 
		double rr = re * c.re - im * c.im, ri = re * c.im + im * c.re;
		re = rr; im = ri;
		return *this;
	}
	inline complex operator*=(double d)
	{ re *= d; return *this; }

	inline complex operator/=(double d)
	{ re /= d; im /= d; return *this; }

	inline complex conj()
	{ return complex(re, -im);  }

	inline double sq_norm() const
	{ return re * re + im * im; }
	inline double norm() const
	{ return std::sqrt(this->sq_norm()); }
};


const complex M_I(0.0, 1.0);

inline complex euler(double norm, double arg)
{ return complex(norm * std::cos(arg), norm * std::sin(arg)); }

inline complex exp(const complex &c)
{
	return euler(std::exp(c.re), c.im);
}


inline double sq_norm(double d)
{
	return d * d;
}

inline double sq_norm(const complex &c)
{
	return c.re * c.re + c.im * c.im;
}

inline double norm(double d)
{ return d; }

inline double norm(const complex &c)
{ return sqrt(sq_norm(c)); }

inline complex normalize(const complex &c)
{ return complex(c.re / norm(c), c.im / norm(c)); }

inline complex conjugate(const complex &c)
{ return complex(c.re, -c.im); }

inline double conjugate(double c)
{ return c; }

inline complex operator+(const complex &c1, const complex &c2)
{ return complex(c1.re + c2.re, c1.im + c2.im); }

inline complex operator+(const complex &c1, double d)
{ return complex(c1.re + d); }

inline complex operator+(double d, const complex &c)
{ return c + d; }

inline complex operator-(const complex &c1, const complex &c2)
{ return complex(c1.re - c2.re, c1.im - c2.im); }

inline complex operator-(const complex &c1, double d)
{ return complex(c1.re - d); }

inline complex operator-(double d, const complex &c)
{ return complex(d - c.re, -c.im); }

inline complex operator*(const complex &c1, const complex &c2)
{ return complex(c1.re * c2.re - c1.im * c2.im,  c1.re * c2.im + c1.im * c2.re);}

inline complex operator*(const complex &c, double d)
{ return complex(c.re * d, c.im * d); }

inline complex operator*(double d, const complex &c)
{ return c * d; }

inline complex operator/(const complex &c, const complex &d)
{ return c * complex(d.re / sq_norm(d), -d.im / sq_norm(d)); }

inline complex operator/(double c, const complex &d)
{ return c * complex(d.re / sq_norm(d), -d.im / sq_norm(d)); }

inline complex operator/(const complex &c, double d)
{ return complex(c.re / d, c.im / d); }

inline bool operator==(const complex &c1, const complex &c2)
{ return (c1.re == c2.re) && (c1.im == c2.im); }

inline bool operator<(const complex &c1, const complex &c2)
{ return sq_norm(c1) < sq_norm(c2); }

inline bool operator<(const complex &c, double d)
{ return sq_norm(c) < d; }

inline bool operator<(double d, const complex &c)
{ return c < d; }

inline bool operator<=(const complex &c1, const complex &c2)
{ return sq_norm(c1) <= sq_norm(c2); }

inline bool operator<=(const complex &c, double d)
{ return sq_norm(c) <= d; }

inline bool operator<=(double d, const complex &c)
{ return c <= d; }

inline bool operator>(const complex &c1, const complex &c2)
{ return sq_norm(c1) > sq_norm(c2); }

inline bool operator>(const complex &c, double d)
{ return sq_norm(c) > d; }

inline bool operator>(double d, complex &c)
{ return c > d; }

inline bool operator>=(const complex &c1, const complex &c2)
{ return sq_norm(c1) >= sq_norm(c2); }

inline bool operator>=(const complex &c, double d)
{ return sq_norm(c) >= d; }

inline bool operator>=(double d, const complex &c)
{ return c >= d; }

inline complex operator+(const complex &c)
{ return c; }

inline complex operator-(const complex &c)
{ return complex(-c.re, -c.im); }


inline double fabs(const complex &c)
{
	return norm(c);
}

inline double fabs(double c)
{
	return std::fabs(c);
}

// to be done
inline double sqrt(const complex &c)
{
	return norm(c);
}

inline double sqrt(double c)
{
	return std::sqrt(c);
}

inline double copysign(double c, double d)
{ return std::copysign(c, d);}


inline bool feq(const complex &l, const complex &r)
{ return (fabs(l - r) < FERR); }

inline bool feq(double l, double r)
{ return (fabs(l - r) < FERR); }

inline double real(const complex &c)
{ return c.re; }

inline double real(complex &c)
{ return c.re; }

inline double real(const double &d)
{ return d; }

inline double real(double &d)
{ return d; }

inline double im(const double d)
{ return 0.0; }

inline double im(const complex &c)
{ return c.im; }

inline double im(complex &c)
{ return c.im; }

inline complex c(double re, double im)
{ return complex(re, im); }


} // namespace
#endif

