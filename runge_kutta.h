/*

Copyright 2024 Simon Peter

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H


#include <stdint.h>
#include <stddef.h>
#include <assert.h>


 /* It is assumed that an initial value problem, of the (explicit) form
		 dx/dt = f ( t, x, u )
		 x(t0) = x0
	 is being solved.
	 x and f are vector entities
	 u can be any entitie, e.g. a struct entitie to handle different kind of input signals
 */
template <typename T, typename Tu> class RkVec {
public:
	RkVec(void(*_f)(T, const T*, const Tu*, T*), const size_t _m, const size_t _s, const T* _a, const  T* _b, const T* _c)
		: s(_s), a(_a), b(_b), c(_c), m(_m), f(_f)
	{
		xn = new T[m];
		clearState();
		buf0 = new T[m];
		buf1 = new T[m];
		buf2 = new T[m];

		ki = new T*[s];
		for (size_t i = 0; i < s; ++i)
			ki[i] = new T[m];

		// Check method requirements
		T sum;
		constexpr T eps = T(1e-8);
		sum = T(0);		// sum_{i=1}^{s}(bi) = 1
		for (size_t i = 0; i < s; ++i)
			sum = sum + b[i];
		assert((sum - T(1) > T(0) ? sum : -sum) < eps);
	}

	~RkVec()
	{
		delete[] xn;
		delete[] buf0;
		delete[] buf1;
		delete[] buf2;
		for (size_t i = 0; i < s; ++i)
			delete[] ki[i];
		delete[] ki;
	}

	// Set system state: State variables with a time step
	void setState(const T* xiNew, T tNew) {
		tn = tNew;
		for (size_t i = 0; i < m; ++i)
			xn[i] = xiNew[i];
	}

	// Set system state: One state variables at the index position
	void setState(const T xNew, size_t i) {
		xn[i] = xNew;
	}

	// Set system state: Time step
	void setState(T tNew) {
		tn = tNew;
	}

	// Sets the system to the original state (state variables and time)
	void clearState() {
		tn = T(0);
		vecClr(xn, m);
	}

	const T* getState() const { return xn; }
	T getTime() const { return tn; }

	// Execute an iteration (calculate the next system state with the passed step size) 
	//  The passed control variable u has to be constat during the step size
	void solve(const Tu* u, T h)
	{
		/*  The family of explicit Runge-Kutta methods is given by
			   xn+1 = xn + h * sum_{i=1}^{s}(bi * ki)
			with
			   k1 = f(tn, xn)
			   k2 = f(tn + c2*h, xn+(a21*k1)*h)
			   k3 = f(tn + c3*h, xn+(a31*k1 + a32*k2)*h)
				 ...
			   ks = f(tn + cs*h, xn+(as1*k1 + as2*k2 + ... + ass-1*ks-1)*h)

		     f is expanded by u, to handle input signals (e.g. control variables)
		*/
		for (size_t i = 0; i < s; ++i) {
			T tEval = tn + c[i] * h;
			vecAssign(buf0, xn, m);
			vecClr(buf2, m);

			for (size_t j = 0; j < i; ++j) {
				vecScalMul(buf1, ki[j], a[i*s + j], m);
				vecAdd(buf2, buf1, m);
			}

			if (i > 0) {
				vecScalMul(buf2, buf2, h, m);
				vecAdd(buf0, buf2, m);
			}

			// Calculate i'th stage dx/dt = f ( t, x, u )
			f(tEval, buf0, u, ki[i]);
		}

		vecClr(buf0, m);
		for (size_t i = 0; i < s; ++i) {
			vecScalMul(buf1, ki[i], b[i], m);
			vecAdd(buf0, buf1, m);
		}
		vecScalMul(buf0, buf0, h, m);
		vecAdd(buf0, xn, m);

		vecAssign(xn, buf0, m);
		tn += h;
	}

private:
	// v = 0
	inline void vecClr(/*out*/ T* v, size_t m) {
		for (size_t i = 0; i < m; ++i)
			v[i] = T(0);
	}

	// dst = src
	inline void vecAssign(/*out*/ T* dst, const T* src, size_t m) {
		for (size_t i = 0; i < m; ++i)
			dst[i] = src[i];
	}

	// v = src * s
	inline void vecScalMul(/*out*/ T* dst, const T* src, T s, size_t m) {
		for (size_t i = 0; i < m; ++i)
			dst[i] = src[i] * s;
	}

	// v = v + src;
	inline void vecAdd(/*inOut*/ T* v, const T* src, size_t m) {
		for (size_t i = 0; i < m; ++i)
			v[i] += src[i];
	}

	/* Runge-Kutta tableau (also Butcher scheme)
		 c |  A         Dimension:   aij  for 1 <= j <= i <= s
		---|----                     bi   for i = 1,2,..,s
		   | b^T                     ci   for i = 1,2,..,s    with c[1] = 0
		e.g.
		  Forward Euler
			c = [0]  b^T = [1]  A = [0]
		  Midpoint method
			c = [0 1/2]  b^T = [0 1]  A = [0, 0; 1/2, 0]
	*/
	const size_t s;	// Number of stages
	const T* a;		// [s,s]
	const T* b;		// [s]
	const T* c;		// [s]

	const size_t m;	// Dimension
	T tn;			// Timestamp
	T* xn;			// [m] State variables
	T** ki;			// [s,m] Solution of the stages
	T* buf0;		// Buffers for solving (Prevents dynamic memory acquisition)
	T* buf1;
	T* buf2;

	// Explicit vectorial nonlinear system
	void(*f)(T t, const T* x, const Tu* u, /*out*/ T* dxdt);
};



#endif
