/* Include files */

#include "Control_linear_sfun.h"
#include "c1_Control_linear.h"
#include "c2_Control_linear.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
uint32_T _Control_linearMachineNumber_;
real_T _sfTime_;

/* Function Declarations */

/* Function Definitions */
void Control_linear_initializer(void)
{
}

void Control_linear_terminator(void)
{
}

/* SFunction Glue Code */
unsigned int sf_Control_linear_method_dispatcher(SimStruct *simstructPtr,
  unsigned int chartFileNumber, const char* specsCksum, int_T method, void *data)
{
  if (chartFileNumber==1) {
    c1_Control_linear_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  if (chartFileNumber==2) {
    c2_Control_linear_method_dispatcher(simstructPtr, method, data);
    return 1;
  }

  return 0;
}

unsigned int sf_Control_linear_process_testpoint_info_call( int nlhs, mxArray *
  plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char machineName[128];
  if (nrhs < 3 || !mxIsChar(prhs[0]) || !mxIsChar(prhs[1]))
    return 0;

  /* Possible call to get testpoint info. */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_testpoint_info"))
    return 0;
  mxGetString(prhs[1], machineName, sizeof(machineName)/sizeof(char));
  machineName[(sizeof(machineName)/sizeof(char)-1)] = '\0';
  if (!strcmp(machineName, "Control_linear")) {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
    switch (chartFileNumber) {
     case 1:
      {
        extern mxArray *sf_c1_Control_linear_get_testpoint_info(void);
        plhs[0] = sf_c1_Control_linear_get_testpoint_info();
        break;
      }

     case 2:
      {
        extern mxArray *sf_c2_Control_linear_get_testpoint_info(void);
        plhs[0] = sf_c2_Control_linear_get_testpoint_info();
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }

    return 1;
  }

  return 0;

#else

  return 0;

#endif

}

unsigned int sf_Control_linear_process_check_sum_call( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[20];
  if (nrhs<1 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the checksum */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"sf_get_check_sum"))
    return 0;
  plhs[0] = mxCreateDoubleMatrix( 1,4,mxREAL);
  if (nrhs>1 && mxIsChar(prhs[1])) {
    mxGetString(prhs[1], commandName,sizeof(commandName)/sizeof(char));
    commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
    if (!strcmp(commandName,"machine")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1996610623U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2988766787U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3145317866U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(4226067495U);
    } else if (!strcmp(commandName,"exportedFcn")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0U);
    } else if (!strcmp(commandName,"makefile")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3211485813U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(353647284U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1361384603U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(149900454U);
    } else if (nrhs==3 && !strcmp(commandName,"chart")) {
      unsigned int chartFileNumber;
      chartFileNumber = (unsigned int)mxGetScalar(prhs[2]);
      switch (chartFileNumber) {
       case 1:
        {
          extern void sf_c1_Control_linear_get_check_sum(mxArray *plhs[]);
          sf_c1_Control_linear_get_check_sum(plhs);
          break;
        }

       case 2:
        {
          extern void sf_c2_Control_linear_get_check_sum(mxArray *plhs[]);
          sf_c2_Control_linear_get_check_sum(plhs);
          break;
        }

       default:
        ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(0.0);
        ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(0.0);
      }
    } else if (!strcmp(commandName,"target")) {
      ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(784597580U);
      ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(508478511U);
      ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1608615269U);
      ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1082679159U);
    } else {
      return 0;
    }
  } else {
    ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(2174650185U);
    ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2579002285U);
    ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3728464309U);
    ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2040521739U);
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Control_linear_autoinheritance_info( int nlhs, mxArray * plhs[],
  int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[32];
  char aiChksum[64];
  if (nrhs<3 || !mxIsChar(prhs[0]) )
    return 0;

  /* Possible call to get the autoinheritance_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_autoinheritance_info"))
    return 0;
  mxGetString(prhs[2], aiChksum,sizeof(aiChksum)/sizeof(char));
  aiChksum[(sizeof(aiChksum)/sizeof(char)-1)] = '\0';

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        if (strcmp(aiChksum, "N7OrUlUXiL1BFg3aiY2VxF") == 0) {
          extern mxArray *sf_c1_Control_linear_get_autoinheritance_info(void);
          plhs[0] = sf_c1_Control_linear_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     case 2:
      {
        if (strcmp(aiChksum, "OdG4ZtwXJT9qhAnC30WsMF") == 0) {
          extern mxArray *sf_c2_Control_linear_get_autoinheritance_info(void);
          plhs[0] = sf_c2_Control_linear_get_autoinheritance_info();
          break;
        }

        plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

unsigned int sf_Control_linear_get_eml_resolved_functions_info( int nlhs,
  mxArray * plhs[], int nrhs, const mxArray * prhs[] )
{

#ifdef MATLAB_MEX_FILE

  char commandName[64];
  if (nrhs<2 || !mxIsChar(prhs[0]))
    return 0;

  /* Possible call to get the get_eml_resolved_functions_info */
  mxGetString(prhs[0], commandName,sizeof(commandName)/sizeof(char));
  commandName[(sizeof(commandName)/sizeof(char)-1)] = '\0';
  if (strcmp(commandName,"get_eml_resolved_functions_info"))
    return 0;

  {
    unsigned int chartFileNumber;
    chartFileNumber = (unsigned int)mxGetScalar(prhs[1]);
    switch (chartFileNumber) {
     case 1:
      {
        extern const mxArray
          *sf_c1_Control_linear_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c1_Control_linear_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     case 2:
      {
        extern const mxArray
          *sf_c2_Control_linear_get_eml_resolved_functions_info(void);
        mxArray *persistentMxArray = (mxArray *)
          sf_c2_Control_linear_get_eml_resolved_functions_info();
        plhs[0] = mxDuplicateArray(persistentMxArray);
        mxDestroyArray(persistentMxArray);
        break;
      }

     default:
      plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
    }
  }

  return 1;

#else

  return 0;

#endif

}

void Control_linear_debug_initialize(void)
{
  _Control_linearMachineNumber_ = sf_debug_initialize_machine("Control_linear",
    "sfun",0,2,0,0,0);
  sf_debug_set_machine_event_thresholds(_Control_linearMachineNumber_,0,0);
  sf_debug_set_machine_data_thresholds(_Control_linearMachineNumber_,0);
}

void Control_linear_register_exported_symbols(SimStruct* S)
{
}

static mxArray* sRtwOptimizationInfoStruct= NULL;
mxArray* load_Control_linear_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct==NULL) {
    sRtwOptimizationInfoStruct = sf_load_rtw_optimization_info("Control_linear",
      "Control_linear");
    mexMakeArrayPersistent(sRtwOptimizationInfoStruct);
  }

  return(sRtwOptimizationInfoStruct);
}

void unload_Control_linear_optimization_info(void)
{
  if (sRtwOptimizationInfoStruct!=NULL) {
    mxDestroyArray(sRtwOptimizationInfoStruct);
    sRtwOptimizationInfoStruct = NULL;
  }
}
