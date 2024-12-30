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

#define WITH_IOS 1

#include <iostream>
#if WITH_IOS
	#include <iomanip>
	#include <sstream>
	#include <iostream>
	#include <fstream>
#endif

#include "runge_kutta.h"

typedef float calcT;

int main()
{
	auto f = [](calcT t, const calcT* x, const void* u, /*out*/ calcT* dxdt) -> void {
		dxdt[0] = sin(t)*sin(t)*x[0];
	};

	// Start value
	const calcT t0{ calcT(0.0) };
	const calcT x0[1]{ calcT(2.0) };

	// Forward euler
	const calcT* aEuler = nullptr;
	const calcT bEuler[1]{ calcT(1.0) };
	const calcT cEuler[1]{ calcT(0.0) };
	const size_t sEuler = sizeof(bEuler) / sizeof(bEuler[0]);

	// Heun solver matrix
	const calcT aHeun[2][2]{ { calcT(0.0), calcT(0.0) },
								{ calcT(1.0), calcT(0.0) } };
	const calcT bHeun[2]{ calcT(0.5), calcT(0.5) };
	const calcT cHeun[2]{ calcT(0.0), calcT(1.0) };
	const size_t sHeun = sizeof(bHeun) / sizeof(bHeun[0]);

	// RK4 solver matrix
	const calcT aRk4[4][4]{ {calcT(0.0), calcT(0.0), calcT(0.0), calcT(0.0)},
								{calcT(0.5), calcT(0.0), calcT(0.0), calcT(0.0)},
								{calcT(0.0), calcT(0.5), calcT(0.0), calcT(0.0)},
								{calcT(0.0), calcT(0.0), calcT(1.0), calcT(0.0)} };
	const calcT bRk4[4]{ calcT(1.0 / 6.0), calcT(1.0 / 3.0), calcT(1.0 / 3.0), calcT(1.0 / 6.0) };
	const calcT cRk4[4]{ calcT(0.0), calcT(0.5), calcT(0.5), calcT(1.0) };
	const size_t sRk4 = sizeof(bRk4) / sizeof(bRk4[0]);

	const size_t solvers = 3;
	const char* names[solvers + 1]{ "Explicit Euler", "2nd order Runge-Kutta",  "4th order Runge-Kutta", nullptr };
	RkVec<calcT, void> rkVec[solvers]{ RkVec<calcT, void>(f, 1, sEuler, aEuler, &bEuler[0], &cEuler[0]),
											RkVec<calcT, void>(f, 1, sHeun, &aHeun[0][0], &bHeun[0], &cHeun[0]),
											RkVec<calcT, void>(f, 1, sRk4, &aRk4[0][0], &bRk4[0], &cRk4[0]) };

	calcT end{ calcT(5.0) };
	calcT dt{ calcT(0.5) };
	size_t runs = size_t(end / dt) + 1;
	calcT* ti = new calcT[runs];
	calcT* xi = new calcT[solvers*runs];

	for (size_t k = 0; k < solvers; ++k) {
		rkVec[k].setState(&x0[0], t0);
		for (size_t j = 0; j < runs; ++j) {
			const calcT* xiLast = rkVec[k].getState();
			const calcT tiLast = rkVec[k].getTime();
			if (k == 0)
				ti[j] = tiLast;
			xi[k*runs + j] = *xiLast;
			rkVec[k].solve(nullptr, dt);
		}
	}

#if WITH_IOS
	std::stringstream foutPlt;
	foutPlt << "\"t\",";
	for (size_t k = 0; k < solvers; ++k) {
		foutPlt << "\"" << names[k] << "\"";
		if (k < solvers - 1)
			foutPlt << ",";
		else
			foutPlt << "\n";
	}
	for (size_t j = 0; j < runs; ++j) {
		foutPlt << ti[j] << ",";
		for (size_t k = 0; k < solvers; ++k) {
			foutPlt << xi[k*runs + j];
			if (k < solvers - 1)
				foutPlt << ",";
			else
				foutPlt << "\n";
		}
	}
	std::ofstream fout1;
	fout1.open("wikiPlt.csv");
	fout1 << foutPlt.rdbuf();
	fout1.close();
#endif
}