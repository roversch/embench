/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * Developed at the Institute of Measurement, Control, and Microtechnology,
 * Ulm University. All rights reserved.
 *
 * GRAMPC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * GRAMPC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>
 *
 *                                           _T
 *                                          /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *             .
 *      s.t.   x(t) = f(t0+t,x(t),u(t),p), x(0) = x0
 *             h(x)  <= 0
 *             u_min <= u(t) <= u_max
 *
 */

#include "probfct.h"

#define POW2(a) ((a)*(a))

/* Casadi generated functions */

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int int
#endif

static casadi_real casadi_sq(casadi_real x) {
    return x*x;
}

/* pendulum:(i0[4],i1)->(o0[4]) */
static int pendulum(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1, a2, a3, a4, a5, a6, a7;
  a0=arg[0] ? arg[0][2] : 0;
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[0] ? arg[0][3] : 0;
  if (res[0]!=0) res[0][1]=a0;
  a1=-8.0000000000000016e-02;
  a2=arg[0] ? arg[0][1] : 0;
  a3=sin(a2);
  a3=(a1*a3);
  a4=casadi_sq(a0);
  a3=(a3*a4);
  a4=9.8100000000000009e-01;
  a5=cos(a2);
  a4=(a4*a5);
  a5=sin(a2);
  a4=(a4*a5);
  a3=(a3+a4);
  a4=arg[1] ? arg[1][0] : 0;
  a3=(a3+a4);
  a5=1.1000000000000001e+00;
  a6=1.0000000000000001e-01;
  a7=cos(a2);
  a7=casadi_sq(a7);
  a6=(a6*a7);
  a5=(a5-a6);
  a3=(a3/a5);
  if (res[0]!=0) res[0][2]=a3;
  a3=cos(a2);
  a1=(a1*a3);
  a3=sin(a2);
  a1=(a1*a3);
  a0=casadi_sq(a0);
  a1=(a1*a0);
  a0=cos(a2);
  a4=(a4*a0);
  a1=(a1+a4);
  a4=1.0791000000000002e+01;
  a2=sin(a2);
  a4=(a4*a2);
  a1=(a1+a4);
  a4=8.0000000000000004e-01;
  a4=(a4*a5);
  a1=(a1/a4);
  if (res[0]!=0) res[0][3]=a1;
  return 0;
}

/* adj_df_dx:(i0[4],i1,i2[4])->(o0[4]) */
static int adj_df_dx(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a2, a3, a4, a5, a6, a7, a8, a9;
  a0=0.;
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[0] ? arg[0][1] : 0;
  a1=cos(a0);
  a2=1.0791000000000002e+01;
  a3=arg[2] ? arg[2][3] : 0;
  a4=8.0000000000000004e-01;
  a5=1.1000000000000001e+00;
  a6=1.0000000000000001e-01;
  a7=cos(a0);
  a8=casadi_sq(a7);
  a8=(a6*a8);
  a5=(a5-a8);
  a8=(a4*a5);
  a9=(a3/a8);
  a10=(a2*a9);
  a1=(a1*a10);
  a10=sin(a0);
  a11=arg[1] ? arg[1][0] : 0;
  a12=(a11*a9);
  a10=(a10*a12);
  a1=(a1-a10);
  a10=cos(a0);
  a12=-8.0000000000000016e-02;
  a13=cos(a0);
  a13=(a12*a13);
  a14=arg[0] ? arg[0][3] : 0;
  a15=casadi_sq(a14);
  a16=(a15*a9);
  a17=(a13*a16);
  a10=(a10*a17);
  a1=(a1+a10);
  a10=sin(a0);
  a17=sin(a0);
  a16=(a17*a16);
  a16=(a12*a16);
  a10=(a10*a16);
  a1=(a1-a10);
  a10=sin(a0);
  a7=(a7+a7);
  a13=(a13*a17);
  a15=(a13*a15);
  a17=cos(a0);
  a17=(a11*a17);
  a15=(a15+a17);
  a17=sin(a0);
  a2=(a2*a17);
  a15=(a15+a2);
  a15=(a15/a8);
  a15=(a15/a8);
  a15=(a15*a3);
  a4=(a4*a15);
  a15=sin(a0);
  a15=(a12*a15);
  a3=casadi_sq(a14);
  a8=(a15*a3);
  a2=9.8100000000000009e-01;
  a17=cos(a0);
  a17=(a2*a17);
  a16=sin(a0);
  a18=(a17*a16);
  a8=(a8+a18);
  a8=(a8+a11);
  a8=(a8/a5);
  a8=(a8/a5);
  a11=arg[2] ? arg[2][2] : 0;
  a8=(a8*a11);
  a4=(a4+a8);
  a6=(a6*a4);
  a7=(a7*a6);
  a10=(a10*a7);
  a1=(a1-a10);
  a10=cos(a0);
  a11=(a11/a5);
  a17=(a17*a11);
  a10=(a10*a17);
  a1=(a1+a10);
  a10=sin(a0);
  a16=(a16*a11);
  a2=(a2*a16);
  a10=(a10*a2);
  a1=(a1-a10);
  a0=cos(a0);
  a3=(a3*a11);
  a12=(a12*a3);
  a0=(a0*a12);
  a1=(a1+a0);
  if (res[0]!=0) res[0][1]=a1;
  a1=arg[2] ? arg[2][0] : 0;
  if (res[0]!=0) res[0][2]=a1;
  a1=(a14+a14);
  a13=(a13*a9);
  a1=(a1*a13);
  a14=(a14+a14);
  a15=(a15*a11);
  a14=(a14*a15);
  a1=(a1+a14);
  a14=arg[2] ? arg[2][1] : 0;
  a1=(a1+a14);
  if (res[0]!=0) res[0][3]=a1;
  return 0;
}

/* adj_df_du:(i0[4],i1,i2[4])->(o0) */
static int adj_df_du(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1, a2, a3, a4, a5;
  a0=arg[0] ? arg[0][1] : 0;
  a1=cos(a0);
  a2=arg[2] ? arg[2][3] : 0;
  a3=8.0000000000000004e-01;
  a4=1.1000000000000001e+00;
  a5=1.0000000000000001e-01;
  a0=cos(a0);
  a0=casadi_sq(a0);
  a5=(a5*a0);
  a4=(a4-a5);
  a3=(a3*a4);
  a2=(a2/a3);
  a1=(a1*a2);
  a2=arg[2] ? arg[2][2] : 0;
  a2=(a2/a4);
  a1=(a1+a2);
  if (res[0]!=0) res[0][0]=a1;
  return 0;
}

/** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = 4;
	*Nu = 1;
	*Np = 0;
	*Nh = 0;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}


/** System function f(t,x,u,p,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
    
    const casadi_real *arg[2];
    arg[0] = x;
    arg[1] = u;
    
    casadi_real *res[1];
    res[0] = out;
    
	pendulum(arg, res, 0, 0, 0);
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	    
    const casadi_real *arg[3];
    arg[0] = x;
    arg[1] = u;
    arg[2] = vec;
    
    casadi_real *res[1];
    res[0] = out;
    
    adj_df_dx(arg, res, 0, 0, 0);
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{

    const casadi_real *arg[3];
    arg[0] = x;
    arg[1] = u;
    arg[2] = vec;
    
    casadi_real *res[1];
    res[0] = out;
    
    adj_df_du(arg, res, 0, 0, 0);
    
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] =
        pCost[0] * POW2(x[0] - xdes[0]) +
		pCost[1] * POW2(x[1] - xdes[1]) +
		pCost[2] * POW2(x[2] - xdes[2]) +
		pCost[3] * POW2(x[3] - xdes[3]) +
		pCost[4] * POW2(u[0] - udes[0]);
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] = 2 * pCost[0] * (x[0] - xdes[0]);
	out[1] = 2 * pCost[1] * (x[1] - xdes[1]);
	out[2] = 2 * pCost[2] * (x[2] - xdes[2]);
	out[3] = 2 * pCost[3] * (x[3] - xdes[3]);
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] = 2 * pCost[4] * (u[0] - udes[0]);
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

    out[0] =
        pCost[5] * POW2(x[0] - xdes[0]) +
		pCost[6] * POW2(x[1] - xdes[1]) +
		pCost[7] * POW2(x[2] - xdes[2]) +
		pCost[8] * POW2(x[3] - xdes[3]);
}

/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum *pCost = (ctypeRNum*)userparam;

	out[0] = 2 * pCost[5] * (x[0] - xdes[0]);
	out[1] = 2 * pCost[6] * (x[1] - xdes[1]);
	out[2] = 2 * pCost[7] * (x[2] - xdes[2]);
	out[3] = 2 * pCost[8] * (x[3] - xdes[3]);
}
/** Gradient dV/dp **/
void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dT **/
void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}


/** Equality constraints g(t,x(t),u(t),p,uperparam) = 0 
    --------------------------------------------------- **/
void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal equality constraints gT(T,x(T),p,uperparam) = 0 
    -------------------------------------------------------- **/
void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal inequality constraints hT(T,x(T),p,uperparam) <= 0 
    ----------------------------------------------------------- **/
void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Additional functions required for semi-implicit systems 
    M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS 
    ------------------------------------------------------- **/
/** Jacobian df/dx in vector form (column-wise) **/
void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dx in vector form (column-wise) **/
void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dt **/
void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian d(dH/dx)/dt  **/
void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mfct(typeRNum *out, typeUSERPARAM *userparam)
{
}
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
{
}
