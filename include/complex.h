/*
complex.h

複素数型と演算 (OpenFDTD include/complex.h の d_complex_t を移植)
MSVC の <complex.h> は C99 非準拠のため自前実装を使う。
*/

#ifndef _PEEC_COMPLEX_H_
#define _PEEC_COMPLEX_H_

#include <math.h>

typedef struct {double r, i;} d_complex_t;

static inline d_complex_t d_complex(double r, double i)
{
	d_complex_t z;

	z.r = r;
	z.i = i;

	return z;
}

static inline d_complex_t d_add(d_complex_t a, d_complex_t b)
{
	return d_complex(a.r + b.r, a.i + b.i);
}

static inline d_complex_t d_sub(d_complex_t a, d_complex_t b)
{
	return d_complex(a.r - b.r, a.i - b.i);
}

static inline d_complex_t d_mul(d_complex_t a, d_complex_t b)
{
	return d_complex((a.r * b.r) - (a.i * b.i), (a.r * b.i) + (a.i * b.r));
}

static inline d_complex_t d_div(d_complex_t a, d_complex_t b)
{
	if ((fabs(b.r) <= 0) && (fabs(b.i) <= 0)) return d_complex(0, 0);
	return d_complex(((a.r * b.r) + (a.i * b.i)) / ((b.r * b.r) + (b.i * b.i)),
	                 ((a.i * b.r) - (a.r * b.i)) / ((b.r * b.r) + (b.i * b.i)));
}

static inline d_complex_t d_rmul(double r, d_complex_t z)
{
	return d_complex(r * z.r, r * z.i);
}

static inline d_complex_t d_inv(d_complex_t z)
{
	const double n2 = (z.r * z.r) + (z.i * z.i);
	if (n2 <= 0) return d_complex(0, 0);
	return d_complex(z.r / n2, -z.i / n2);
}

static inline d_complex_t d_conj(d_complex_t z)
{
	return d_complex(z.r, -z.i);
}

static inline double d_abs(d_complex_t z)
{
	return sqrt((z.r * z.r) + (z.i * z.i));
}

static inline double d_norm(d_complex_t z)
{
	return ((z.r * z.r) + (z.i * z.i));
}

// r * e^(j*deg)
static inline d_complex_t d_polar_deg(double r, double deg)
{
	const double rad = deg * atan(1.0) / 45.0;
	return d_complex(r * cos(rad), r * sin(rad));
}

#endif		// _PEEC_COMPLEX_H_
