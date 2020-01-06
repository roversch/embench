#define S_FUNCTION_NAME NMPCProblem_S
#define INITIALIZE_WRAPPER NMPCProblem_S_Initialize_wrapper
#define OUTPUTS_WRAPPER NMPCProblem_S_Outputs_wrapper
#define NUM_CONTROLS 3
#define NUM_STATES 21
#define NUM_PARAMETERS 0
#define SAMPLE_TIME 0.2
#define NUM_STEPS 40
#define NUM_OBJECTIVE_STATES 15
#define NUM_OBJECTIVE_END_STATES 3

#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#include <time.h>

extern void INITIALIZE_WRAPPER(void **ptr);

extern void OUTPUTS_WRAPPER(const real_T *x, const real_T *x_ref, const real_T *u_ref, const real_T *h_ref, const real_T *h_ref_end,
    const real_T *parmas, const real_T *Q, const real_T *R, const real_T *P, const real_T *S, const real_T *S_END,
    real_T *u, void **ptr);

static void mdlInitializeSizes(SimStruct *S)
{
    DECL_AND_INIT_DIMSINFO(inputDimsInfo);
    DECL_AND_INIT_DIMSINFO(outputDimsInfo);
    ssSetNumSFcnParams(S, 0);
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S))
        return;

    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    if (!ssSetNumInputPorts(S, 3)) return;

    //Input: x
    ssSetInputPortWidth(S,  0, NUM_STATES);
    ssSetInputPortDataType(S, 0, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 0, COMPLEX_NO);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    ssSetInputPortRequiredContiguous(S, 0, 1);

    //Input: h_ref
    //inputDimsInfo.width = NUM_STEPS;
    //ssSetInputPortDimensionInfo(S, 1, &inputDimsInfo);
    ssSetInputPortMatrixDimensions(  S ,1, NUM_OBJECTIVE_STATES, NUM_STEPS);
    ssSetInputPortFrameData(S, 1, FRAME_NO);
    ssSetInputPortDataType(S, 1, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 1, COMPLEX_NO);
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    ssSetInputPortRequiredContiguous(S, 1, 1);

    //Input: h_ref_end
    ssSetInputPortWidth(S,  2, NUM_OBJECTIVE_END_STATES);
    ssSetInputPortDataType(S, 2, SS_DOUBLE);
    ssSetInputPortComplexSignal(S, 2, COMPLEX_NO);
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    ssSetInputPortRequiredContiguous(S, 2, 1);

    if (!ssSetNumOutputPorts(S, 2)) return;

    //Output: u
    ssSetOutputPortWidth(S, 0, NUM_CONTROLS);
    ssSetOutputPortDataType(S, 0, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 0, COMPLEX_NO);
    
    ssSetOutputPortWidth(S, 1, 1);
    ssSetOutputPortDataType(S, 1, SS_DOUBLE);
    ssSetOutputPortComplexSignal(S, 1, COMPLEX_NO);

    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);

    //one D-work vector (pointer)
    ssSetNumDWork(S, 1);
    ssSetDWorkWidth(S, 0, 1);
    ssSetDWorkDataType(S, 0, SS_POINTER);
    ssSetDWorkName(S, 0, "MY_DWORK_VECTOR");

    //execption free code
    ssSetOptions(S, (SS_OPTION_EXCEPTION_FREE_CODE | SS_OPTION_USE_TLC_WITH_ACCELERATOR | SS_OPTION_WORKS_WITH_CODE_REUSE));
}

#define MDL_SET_INPUT_PORT_FRAME_DATA
static void mdlSetInputPortFrameData(SimStruct  *S, int_T port, Frame_T frameData)
{
    ssSetInputPortFrameData(S, port, frameData);
}

static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_START
static void mdlStart(SimStruct *S)
{
   void **ptrs = (void **) ssGetDWork(S,0); /*get a pointer to a pointer*/
   INITIALIZE_WRAPPER(ptrs);
}

#define MDL_SET_INPUT_PORT_DATA_TYPE
static void mdlSetInputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetInputPortDataType( S, 0, dType);
}

#define MDL_SET_OUTPUT_PORT_DATA_TYPE
static void mdlSetOutputPortDataType(SimStruct *S, int port, DTypeId dType)
{
    ssSetOutputPortDataType(S, 0, dType);
}

#define MDL_SET_DEFAULT_PORT_DATA_TYPES
static void mdlSetDefaultPortDataTypes(SimStruct *S)
{
  ssSetInputPortDataType( S, 0, SS_DOUBLE);
  ssSetOutputPortDataType(S, 0, SS_DOUBLE);
}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    const real_T   *x  = (const real_T*) ssGetInputPortSignal(S,0);
    const real_T   *x_ref  = 0;
    const real_T   *u_ref  = 0;
    const real_T   *h_ref  = (const real_T*) ssGetInputPortSignal(S,1);
    const real_T   *h_ref_end  = (const real_T*) ssGetInputPortSignal(S,2);
    const real_T   *param  = 0;
    const real_T   *Q  = 0;
    const real_T   *R  = 0;
    const real_T   *P  = 0;
    const real_T   *S_  = 0;
    const real_T   *S_END  = 0;
    real_T         *u  = (real_T *)  ssGetOutputPortRealSignal(S,0);
    real_T         *cpu_time  = (real_T *)  ssGetOutputPortRealSignal(S,1);

    
    void **ptrs = (void **) ssGetDWork(S,0); /*get a pointer to a pointer*/

    struct timespec tic, toc;
    
    clock_gettime(CLOCK_MONOTONIC, &tic);
    OUTPUTS_WRAPPER(x, x_ref, u_ref, h_ref, h_ref_end, param, Q, R, P, S_, S_END, u, ptrs);
    clock_gettime(CLOCK_MONOTONIC, &toc);

    
    *cpu_time = (toc.tv_sec - tic.tv_sec) + (toc.tv_nsec - tic.tv_nsec) / 1.0e9;
}

static void mdlTerminate(SimStruct *S)
{
   /*FIXME free memory!*/
}

#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
