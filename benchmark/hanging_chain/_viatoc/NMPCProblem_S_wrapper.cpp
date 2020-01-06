#include "NMPCProblem.h"

#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif

void NMPCProblem_S_Initialize_wrapper(void **ptr)
{
    *ptr = (void *) 0;
}

void NMPCProblem_S_Outputs_wrapper(const real_T *x, const real_T *x_ref, const real_T *u_ref, const real_T *h_ref, const real_T *h_ref_end, const real_T *params, const real_T *Q, const real_T *R, const real_T *P, const real_T *S, const real_T *S_END, real_T *u, void **ptr)
{
   int i,j;

   NMPCProblem *nmpc = (NMPCProblem *) *ptr;

   //Initialize if pointer is zero!
   if(nmpc == NULL)
   {
        //allocate memory (FIXME not freed enywhere!)
        nmpc = new NMPCProblem();

        //save pointer to pointer
        *ptr = (void *) nmpc;

       //set initial control trajectory
       for(i = 0; i < nmpc->numControls*nmpc->numSteps; i++)
           nmpc->u[i] = 0.0;

       //set initial state trajectory
       for(j = 0; j < nmpc->numSteps; j++)
           for(i = 0; i < nmpc->numStates; i++)
               nmpc->x[j*nmpc->numStates+i] = x[i];

   }

   //set current state
   for(i = 0; i < nmpc->numStates; i++)
       nmpc->x[i] = x[i];

   //set x_ref
   if(x_ref != 0 && nmpc->x_ref != 0)
       for(i = 0; i < nmpc->numStates*nmpc->numSteps; i++)
           nmpc->x_ref[i] = x_ref[i];

   //set u_ref
   if(u_ref != 0 && nmpc->u_ref != 0)
       for(i = 0; i < nmpc->numControls*nmpc->numSteps; i++)
           nmpc->u_ref[i] = u_ref[i];

   //set h_ref
   if(h_ref != 0 && nmpc->h_ref != 0)
       for(i = 0; i < nmpc->numObjectiveStates*nmpc->numSteps; i++)
           nmpc->h_ref[i] = h_ref[i];

   //set h_ref_end
   if(h_ref_end != 0 && nmpc->h_ref_end != 0)
       for(i = 0; i < nmpc->numObjectiveEndStates; i++)
           nmpc->h_ref_end[i] = h_ref_end[i];

   //set params
   if(params != 0 && nmpc->p != 0)
       for(i = 0; i < nmpc->numParameters; i++)
           nmpc->p[i] = params[i];

   //set Q
   if(Q != 0 && nmpc->Q != 0)
       for(i = 0; i < nmpc->numStates*nmpc->numStates; i++)
           nmpc->Q[i] = Q[i];

   //set R
   if(R != 0 && nmpc->R != 0)
       for(i = 0; i < nmpc->numControls*nmpc->numControls; i++)
           nmpc->R[i] = R[i];

   //set P
   if(P != 0 && nmpc->P != 0)
       for(i = 0; i < nmpc->numStates*nmpc->numStates; i++)
           nmpc->P[i] = P[i];

   //set S
   if(S != 0 && nmpc->S != 0)
       for(i = 0; i < nmpc->numObjectiveStates*nmpc->numObjectiveStates; i++)
           nmpc->S[i] = S[i];

   //set S_END
   if(S_END != 0 && nmpc->S_END != 0)
       for(i = 0; i < nmpc->numObjectiveEndStates*nmpc->numObjectiveEndStates; i++)
           nmpc->S_END[i] = S_END[i];

   //run NMPC
   nmpc->optimize(20);

   //copy controls as outputs (only first control vector!)
   for(i = 0; i < nmpc->numControls; i++)
       u[i] = nmpc->u[i];
}
