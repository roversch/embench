#include "mpcexport.h"

#define POW(x) ((x)*(x))
#define SQRT(x) sqrt(x)

int main ()
{	
	const double dt = 0.1;

	/* dynamic states */
	DifferentialState x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20;

	/* controls */
	Control u0, u1, u2;

  /* continuous time model */
	DifferentialEquation f;
	f << dot(x0) == x12;
	f << dot(x1) == x13;
	f << dot(x2) == x14;
	f << dot(x3) == x15;
	f << dot(x4) == x16;
	f << dot(x5) == x17;
	f << dot(x6) == x18;
	f << dot(x7) == x19;
	f << dot(x8) == x20;
	f << dot(x9) == u0;
	f << dot(x10) == u1;
	f << dot(x11) == u2;
	f << dot(x12) == 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x0) + POW(x1) + POW(x2)))*x0 + 0.1*(4. - 0.55 / SQRT(POW(x0 - x3) + POW(x1 - x4) + POW(x2 - x5)))*(-x0 + x3));
	f << dot(x13) == 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x0) + POW(x1) + POW(x2)))*x1 + 0.1*(4. - 0.55 / SQRT(POW(x0 - x3) + POW(x1 - x4) + POW(x2 - x5)))*(-x1 + x4));
	f << dot(x14) == -9.81 + 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x0) + POW(x1) + POW(x2)))*x2 + 0.1*(4. - 0.55 / SQRT(POW(x0 - x3) + POW(x1 - x4) + POW(x2 - x5)))*(-x2 + x5));
	f << dot(x15) == 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x0 - x3) + POW(x1 - x4) + POW(x2 - x5)))*(-x0 + x3) + 0.1*(4. - 0.55 / SQRT(POW(x3 - x6) + POW(x4 - x7) + POW(x5 - x8)))*(-x3 + x6));
	f << dot(x16) == 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x0 - x3) + POW(x1 - x4) + POW(x2 - x5)))*(-x1 + x4) + 0.1*(4. - 0.55 / SQRT(POW(x3 - x6) + POW(x4 - x7) + POW(x5 - x8)))*(-x4 + x7));
	f << dot(x17) == -9.81 + 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x0 - x3) + POW(x1 - x4) + POW(x2 - x5)))*(-x2 + x5) + 0.1*(4. - 0.55 / SQRT(POW(x3 - x6) + POW(x4 - x7) + POW(x5 - x8)))*(-x5 + x8));
	f << dot(x18) == 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x3 - x6) + POW(x4 - x7) + POW(x5 - x8)))*(-x3 + x6) + 0.1*(4. - 0.55 / SQRT(POW(x6 - x9) + POW(x7 - x10) + POW(x8 - x11)))*(-x6 + x9));
	f << dot(x19) == 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x3 - x6) + POW(x4 - x7) + POW(x5 - x8)))*(-x4 + x7) + 0.1*(4. - 0.55 / SQRT(POW(x6 - x9) + POW(x7 - x10) + POW(x8 - x11)))*(-x7 + x10));
	f << dot(x20) == -9.81 + 8.88888888888889*(-0.1*(4. - 0.55 / SQRT(POW(x3 - x6) + POW(x4 - x7) + POW(x5 - x8)))*(-x5 + x8) + 0.1*(4. - 0.55 / SQRT(POW(x6 - x9) + POW(x7 - x10) + POW(x8 - x11)))*(-x8 + x11));

	/* costfunction weights */
	Matrix Q = eye(15);
	Matrix P = eye(3);
	Function intCost;
	Function terCost;
	intCost << x9 -7.5;
	intCost << x10;
	intCost << x11;
	intCost << x12;
	intCost << x13;
	intCost << x14;
	intCost << x15;
	intCost << x16;
	intCost << x17;
	intCost << x18;
	intCost << x19;
	intCost << x20;
	intCost << u0;
	intCost << u1;
	intCost << u2;
	terCost << x9 -7.5;
	terCost << x10;
	terCost << x11;
	Q(0, 0) = 2.5;
	Q(1, 1) = 2.5;
	Q(2, 2) = 2.5;
	Q(3, 3) = 25;
	Q(4, 4) = 25;
	Q(5, 5) = 25;
	Q(6, 6) = 25;
	Q(7, 7) = 25;
	Q(8, 8) = 25;
	Q(9, 9) = 25;
	Q(10, 10) = 25;
	Q(11, 11) = 25;
	Q(12, 12) = 0.01;
	Q(13, 13) = 0.01;
	Q(14, 14) = 0.01;
	P(0, 0) = 2.5;
	P(1, 1) = 2.5;
	P(2, 2) = 2.5;

	MPCexport mpc(0.0, 8.0, 40);
	mpc.minimizeLSQ(Q, intCost);
	mpc.minimizeLSQEndTerm(P, terCost);
	mpc.subjectTo(f);
	mpc.subjectTo(-1.0 <= u0 <= 1.0);
	mpc.subjectTo(-1.0 <= u1 <= 1.0);
	mpc.subjectTo(-1.0 <= u2 <= 1.0);

	mpc.set(INTEGRATOR_TYPE, INT_RK4);
	mpc.set(EXPORT_CPP, YES);
	mpc.set(EXPORT_SFUNCTION, YES);

	mpc.exportCode("_viatoc");

	return 0;
}

