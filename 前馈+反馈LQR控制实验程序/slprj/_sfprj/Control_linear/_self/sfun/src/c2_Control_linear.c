/* Include files */

#include "blascompat32.h"
#include "Control_linear_sfun.h"
#include "c2_Control_linear.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "Control_linear_sfun_debug_macros.h"

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c2_debug_family_names[6] = { "nargin", "nargout", "f1", "f2",
  "v1", "v2" };

/* Function Declarations */
static void initialize_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance);
static void initialize_params_c2_Control_linear
  (SFc2_Control_linearInstanceStruct *chartInstance);
static void enable_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance);
static void disable_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance);
static void c2_update_debugger_state_c2_Control_linear
  (SFc2_Control_linearInstanceStruct *chartInstance);
static void ext_mode_exec_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c2_Control_linear
  (SFc2_Control_linearInstanceStruct *chartInstance);
static void set_sim_state_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_st);
static void finalize_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance);
static void sf_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance);
static void initSimStructsc2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber);
static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData);
static real_T c2_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_v2, const char_T *c2_identifier);
static real_T c2_b_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static real_T c2_mpower(SFc2_Control_linearInstanceStruct *chartInstance, real_T
  c2_a);
static void c2_eml_error(SFc2_Control_linearInstanceStruct *chartInstance);
static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData);
static int32_T c2_c_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData);
static uint8_T c2_d_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_is_active_c2_Control_linear, const char_T
  *c2_identifier);
static uint8_T c2_e_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId);
static void init_dsm_address_info(SFc2_Control_linearInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance)
{
  int32_T *c2_sfEvent;
  uint8_T *c2_is_active_c2_Control_linear;
  c2_is_active_c2_Control_linear = (uint8_T *)ssGetDWork(chartInstance->S, 3);
  c2_sfEvent = (int32_T *)ssGetDWork(chartInstance->S, 0);
  *c2_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  *c2_is_active_c2_Control_linear = 0U;
}

static void initialize_params_c2_Control_linear
  (SFc2_Control_linearInstanceStruct *chartInstance)
{
}

static void enable_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c2_update_debugger_state_c2_Control_linear
  (SFc2_Control_linearInstanceStruct *chartInstance)
{
}

static void ext_mode_exec_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance)
{
  c2_update_debugger_state_c2_Control_linear(chartInstance);
}

static const mxArray *get_sim_state_c2_Control_linear
  (SFc2_Control_linearInstanceStruct *chartInstance)
{
  const mxArray *c2_st;
  const mxArray *c2_y = NULL;
  real_T c2_hoistedGlobal;
  real_T c2_u;
  const mxArray *c2_b_y = NULL;
  real_T c2_b_hoistedGlobal;
  real_T c2_b_u;
  const mxArray *c2_c_y = NULL;
  uint8_T c2_c_hoistedGlobal;
  uint8_T c2_c_u;
  const mxArray *c2_d_y = NULL;
  real_T *c2_v1;
  real_T *c2_v2;
  uint8_T *c2_is_active_c2_Control_linear;
  c2_v2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_v1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_is_active_c2_Control_linear = (uint8_T *)ssGetDWork(chartInstance->S, 3);
  c2_st = NULL;
  c2_st = NULL;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_createcellarray(3), FALSE);
  c2_hoistedGlobal = *c2_v1;
  c2_u = c2_hoistedGlobal;
  c2_b_y = NULL;
  sf_mex_assign(&c2_b_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 0, c2_b_y);
  c2_b_hoistedGlobal = *c2_v2;
  c2_b_u = c2_b_hoistedGlobal;
  c2_c_y = NULL;
  sf_mex_assign(&c2_c_y, sf_mex_create("y", &c2_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 1, c2_c_y);
  c2_c_hoistedGlobal = *c2_is_active_c2_Control_linear;
  c2_c_u = c2_c_hoistedGlobal;
  c2_d_y = NULL;
  sf_mex_assign(&c2_d_y, sf_mex_create("y", &c2_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c2_y, 2, c2_d_y);
  sf_mex_assign(&c2_st, c2_y, FALSE);
  return c2_st;
}

static void set_sim_state_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_st)
{
  const mxArray *c2_u;
  boolean_T *c2_doneDoubleBufferReInit;
  real_T *c2_v1;
  real_T *c2_v2;
  uint8_T *c2_is_active_c2_Control_linear;
  c2_v2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_v1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_is_active_c2_Control_linear = (uint8_T *)ssGetDWork(chartInstance->S, 3);
  c2_doneDoubleBufferReInit = (boolean_T *)ssGetDWork(chartInstance->S, 2);
  *c2_doneDoubleBufferReInit = TRUE;
  c2_u = sf_mex_dup(c2_st);
  *c2_v1 = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 0)),
    "v1");
  *c2_v2 = c2_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c2_u, 1)),
    "v2");
  *c2_is_active_c2_Control_linear = c2_d_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c2_u, 2)), "is_active_c2_Control_linear");
  sf_mex_destroy(&c2_u);
  c2_update_debugger_state_c2_Control_linear(chartInstance);
  sf_mex_destroy(&c2_st);
}

static void finalize_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance)
{
}

static void sf_c2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance)
{
  real_T c2_hoistedGlobal;
  real_T c2_b_hoistedGlobal;
  real_T c2_f1;
  real_T c2_f2;
  uint32_T c2_debug_family_var_map[6];
  real_T c2_nargin = 2.0;
  real_T c2_nargout = 2.0;
  real_T c2_v1;
  real_T c2_v2;
  real_T c2_b;
  real_T c2_y;
  real_T c2_b_b;
  real_T c2_b_y;
  int32_T *c2_sfEvent;
  real_T *c2_b_f1;
  real_T *c2_b_v1;
  real_T *c2_b_f2;
  real_T *c2_b_v2;
  c2_b_v2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c2_b_f2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c2_b_v1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c2_b_f1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c2_sfEvent = (int32_T *)ssGetDWork(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 1U, *c2_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c2_b_f1, 0U);
  _SFD_DATA_RANGE_CHECK(*c2_b_v1, 1U);
  _SFD_DATA_RANGE_CHECK(*c2_b_f2, 2U);
  _SFD_DATA_RANGE_CHECK(*c2_b_v2, 3U);
  *c2_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 1U, *c2_sfEvent);
  c2_hoistedGlobal = *c2_b_f1;
  c2_b_hoistedGlobal = *c2_b_f2;
  c2_f1 = c2_hoistedGlobal;
  c2_f2 = c2_b_hoistedGlobal;
  sf_debug_symbol_scope_push_eml(0U, 6U, 6U, c2_debug_family_names,
    c2_debug_family_var_map);
  sf_debug_symbol_scope_add_eml_importable(&c2_nargin, 0U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  sf_debug_symbol_scope_add_eml_importable(&c2_nargout, 1U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  sf_debug_symbol_scope_add_eml(&c2_f1, 2U, c2_sf_marshallOut);
  sf_debug_symbol_scope_add_eml(&c2_f2, 3U, c2_sf_marshallOut);
  sf_debug_symbol_scope_add_eml_importable(&c2_v1, 4U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  sf_debug_symbol_scope_add_eml_importable(&c2_v2, 5U, c2_sf_marshallOut,
    c2_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, *c2_sfEvent, 3);
  c2_b = c2_mpower(chartInstance, c2_f1);
  c2_y = 1.677 * c2_b;
  c2_v1 = c2_y - 0.3538;
  _SFD_EML_CALL(0U, *c2_sfEvent, 4);
  c2_b_b = c2_mpower(chartInstance, c2_f2);
  c2_b_y = 1.677 * c2_b_b;
  c2_v2 = c2_b_y - 0.3538;
  _SFD_EML_CALL(0U, *c2_sfEvent, -4);
  sf_debug_symbol_scope_pop();
  *c2_b_v1 = c2_v1;
  *c2_b_v2 = c2_v2;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 1U, *c2_sfEvent);
  sf_debug_check_for_state_inconsistency(_Control_linearMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc2_Control_linear(SFc2_Control_linearInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c2_machineNumber, uint32_T
  c2_chartNumber)
{
}

static const mxArray *c2_sf_marshallOut(void *chartInstanceVoid, void *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  real_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc2_Control_linearInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(real_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static real_T c2_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_v2, const char_T *c2_identifier)
{
  real_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_v2), &c2_thisId);
  sf_mex_destroy(&c2_v2);
  return c2_y;
}

static real_T c2_b_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  real_T c2_y;
  real_T c2_d0;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_d0, 1, 0, 0U, 0, 0U, 0);
  c2_y = c2_d0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_v2;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  real_T c2_y;
  SFc2_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc2_Control_linearInstanceStruct *)chartInstanceVoid;
  c2_v2 = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_v2), &c2_thisId);
  sf_mex_destroy(&c2_v2);
  *(real_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

const mxArray *sf_c2_Control_linear_get_eml_resolved_functions_info(void)
{
  const mxArray *c2_nameCaptureInfo;
  c2_ResolvedFunctionInfo c2_info[8];
  c2_ResolvedFunctionInfo (*c2_b_info)[8];
  const mxArray *c2_m0 = NULL;
  int32_T c2_i0;
  c2_ResolvedFunctionInfo *c2_r0;
  c2_nameCaptureInfo = NULL;
  c2_nameCaptureInfo = NULL;
  c2_b_info = (c2_ResolvedFunctionInfo (*)[8])c2_info;
  (*c2_b_info)[0].context = "";
  (*c2_b_info)[0].name = "mpower";
  (*c2_b_info)[0].dominantType = "double";
  (*c2_b_info)[0].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  (*c2_b_info)[0].fileTimeLo = 1286797242U;
  (*c2_b_info)[0].fileTimeHi = 0U;
  (*c2_b_info)[0].mFileTimeLo = 0U;
  (*c2_b_info)[0].mFileTimeHi = 0U;
  (*c2_b_info)[1].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m";
  (*c2_b_info)[1].name = "power";
  (*c2_b_info)[1].dominantType = "double";
  (*c2_b_info)[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m";
  (*c2_b_info)[1].fileTimeLo = 1336500496U;
  (*c2_b_info)[1].fileTimeHi = 0U;
  (*c2_b_info)[1].mFileTimeLo = 0U;
  (*c2_b_info)[1].mFileTimeHi = 0U;
  (*c2_b_info)[2].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  (*c2_b_info)[2].name = "eml_scalar_eg";
  (*c2_b_info)[2].dominantType = "double";
  (*c2_b_info)[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  (*c2_b_info)[2].fileTimeLo = 1286797196U;
  (*c2_b_info)[2].fileTimeHi = 0U;
  (*c2_b_info)[2].mFileTimeLo = 0U;
  (*c2_b_info)[2].mFileTimeHi = 0U;
  (*c2_b_info)[3].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  (*c2_b_info)[3].name = "eml_scalexp_alloc";
  (*c2_b_info)[3].dominantType = "double";
  (*c2_b_info)[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  (*c2_b_info)[3].fileTimeLo = 1330583234U;
  (*c2_b_info)[3].fileTimeHi = 0U;
  (*c2_b_info)[3].mFileTimeLo = 0U;
  (*c2_b_info)[3].mFileTimeHi = 0U;
  (*c2_b_info)[4].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  (*c2_b_info)[4].name = "floor";
  (*c2_b_info)[4].dominantType = "double";
  (*c2_b_info)[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  (*c2_b_info)[4].fileTimeLo = 1286797142U;
  (*c2_b_info)[4].fileTimeHi = 0U;
  (*c2_b_info)[4].mFileTimeLo = 0U;
  (*c2_b_info)[4].mFileTimeHi = 0U;
  (*c2_b_info)[5].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  (*c2_b_info)[5].name = "eml_scalar_floor";
  (*c2_b_info)[5].dominantType = "double";
  (*c2_b_info)[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  (*c2_b_info)[5].fileTimeLo = 1286797126U;
  (*c2_b_info)[5].fileTimeHi = 0U;
  (*c2_b_info)[5].mFileTimeLo = 0U;
  (*c2_b_info)[5].mFileTimeHi = 0U;
  (*c2_b_info)[6].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower";
  (*c2_b_info)[6].name = "eml_error";
  (*c2_b_info)[6].dominantType = "char";
  (*c2_b_info)[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m";
  (*c2_b_info)[6].fileTimeLo = 1305296400U;
  (*c2_b_info)[6].fileTimeHi = 0U;
  (*c2_b_info)[6].mFileTimeLo = 0U;
  (*c2_b_info)[6].mFileTimeHi = 0U;
  (*c2_b_info)[7].context = "";
  (*c2_b_info)[7].name = "mtimes";
  (*c2_b_info)[7].dominantType = "double";
  (*c2_b_info)[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  (*c2_b_info)[7].fileTimeLo = 1289494492U;
  (*c2_b_info)[7].fileTimeHi = 0U;
  (*c2_b_info)[7].mFileTimeLo = 0U;
  (*c2_b_info)[7].mFileTimeHi = 0U;
  sf_mex_assign(&c2_m0, sf_mex_createstruct("nameCaptureInfo", 1, 8), FALSE);
  for (c2_i0 = 0; c2_i0 < 8; c2_i0++) {
    c2_r0 = &c2_info[c2_i0];
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c2_r0->context)), "context", "nameCaptureInfo",
                    c2_i0);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c2_r0->name)), "name", "nameCaptureInfo", c2_i0);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c2_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c2_i0);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", c2_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c2_r0->resolved)), "resolved", "nameCaptureInfo",
                    c2_i0);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c2_i0);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c2_i0);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c2_i0);
    sf_mex_addfield(c2_m0, sf_mex_create("nameCaptureInfo", &c2_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c2_i0);
  }

  sf_mex_assign(&c2_nameCaptureInfo, c2_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c2_nameCaptureInfo);
  return c2_nameCaptureInfo;
}

static real_T c2_mpower(SFc2_Control_linearInstanceStruct *chartInstance, real_T
  c2_a)
{
  real_T c2_b_a;
  real_T c2_c_a;
  real_T c2_ak;
  c2_b_a = c2_a;
  c2_c_a = c2_b_a;
  c2_ak = c2_c_a;
  if (c2_ak < 0.0) {
    c2_eml_error(chartInstance);
  }

  return muDoubleScalarPower(c2_ak, 0.4716);
}

static void c2_eml_error(SFc2_Control_linearInstanceStruct *chartInstance)
{
  int32_T c2_i1;
  static char_T c2_varargin_1[31] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o',
    'o', 'l', 'b', 'o', 'x', ':', 'p', 'o', 'w', 'e', 'r', '_', 'd', 'o', 'm',
    'a', 'i', 'n', 'E', 'r', 'r', 'o', 'r' };

  char_T c2_u[31];
  const mxArray *c2_y = NULL;
  for (c2_i1 = 0; c2_i1 < 31; c2_i1++) {
    c2_u[c2_i1] = c2_varargin_1[c2_i1];
  }

  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", c2_u, 10, 0U, 1U, 0U, 2, 1, 31), FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 1U, 14,
    c2_y));
}

static const mxArray *c2_b_sf_marshallOut(void *chartInstanceVoid, void
  *c2_inData)
{
  const mxArray *c2_mxArrayOutData = NULL;
  int32_T c2_u;
  const mxArray *c2_y = NULL;
  SFc2_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc2_Control_linearInstanceStruct *)chartInstanceVoid;
  c2_mxArrayOutData = NULL;
  c2_u = *(int32_T *)c2_inData;
  c2_y = NULL;
  sf_mex_assign(&c2_y, sf_mex_create("y", &c2_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c2_mxArrayOutData, c2_y, FALSE);
  return c2_mxArrayOutData;
}

static int32_T c2_c_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  int32_T c2_y;
  int32_T c2_i2;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_i2, 1, 6, 0U, 0, 0U, 0);
  c2_y = c2_i2;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void c2_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c2_mxArrayInData, const char_T *c2_varName, void *c2_outData)
{
  const mxArray *c2_sfEvent;
  const char_T *c2_identifier;
  emlrtMsgIdentifier c2_thisId;
  int32_T c2_y;
  SFc2_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc2_Control_linearInstanceStruct *)chartInstanceVoid;
  c2_sfEvent = sf_mex_dup(c2_mxArrayInData);
  c2_identifier = c2_varName;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c2_sfEvent), &c2_thisId);
  sf_mex_destroy(&c2_sfEvent);
  *(int32_T *)c2_outData = c2_y;
  sf_mex_destroy(&c2_mxArrayInData);
}

static uint8_T c2_d_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_is_active_c2_Control_linear, const char_T
  *c2_identifier)
{
  uint8_T c2_y;
  emlrtMsgIdentifier c2_thisId;
  c2_thisId.fIdentifier = c2_identifier;
  c2_thisId.fParent = NULL;
  c2_y = c2_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c2_is_active_c2_Control_linear), &c2_thisId);
  sf_mex_destroy(&c2_is_active_c2_Control_linear);
  return c2_y;
}

static uint8_T c2_e_emlrt_marshallIn(SFc2_Control_linearInstanceStruct
  *chartInstance, const mxArray *c2_u, const emlrtMsgIdentifier *c2_parentId)
{
  uint8_T c2_y;
  uint8_T c2_u0;
  sf_mex_import(c2_parentId, sf_mex_dup(c2_u), &c2_u0, 1, 3, 0U, 0, 0U, 0);
  c2_y = c2_u0;
  sf_mex_destroy(&c2_u);
  return c2_y;
}

static void init_dsm_address_info(SFc2_Control_linearInstanceStruct
  *chartInstance)
{
}

/* SFunction Glue Code */
static uint32_T* sf_get_sfun_dwork_checksum();
void sf_c2_Control_linear_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3946366289U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(4071037826U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3472561624U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2872965269U);
}

mxArray *sf_c2_Control_linear_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("OdG4ZtwXJT9qhAnC30WsMF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

static const mxArray *sf_get_sim_state_info_c2_Control_linear(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[5],T\"v1\",},{M[1],M[8],T\"v2\",},{M[8],M[0],T\"is_active_c2_Control_linear\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c2_Control_linear_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc2_Control_linearInstanceStruct *chartInstance;
    chartInstance = (SFc2_Control_linearInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (_Control_linearMachineNumber_,
           2,
           1,
           1,
           4,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_Control_linearMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (_Control_linearMachineNumber_,chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(_Control_linearMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,1,1,0,"f1");
          _SFD_SET_DATA_PROPS(1,2,0,1,"v1");
          _SFD_SET_DATA_PROPS(2,1,1,0,"f2");
          _SFD_SET_DATA_PROPS(3,2,0,1,"v2");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,91);
        _SFD_TRANS_COV_WTS(0,0,0,1,0);
        if (chartAlreadyPresent==0) {
          _SFD_TRANS_COV_MAPS(0,
                              0,NULL,NULL,
                              0,NULL,NULL,
                              1,NULL,NULL,
                              0,NULL,NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c2_sf_marshallOut,(MexInFcnForType)c2_sf_marshallIn);

        {
          real_T *c2_f1;
          real_T *c2_v1;
          real_T *c2_f2;
          real_T *c2_v2;
          c2_v2 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c2_f2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c2_v1 = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          c2_f1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c2_f1);
          _SFD_SET_DATA_VALUE_PTR(1U, c2_v1);
          _SFD_SET_DATA_VALUE_PTR(2U, c2_f2);
          _SFD_SET_DATA_VALUE_PTR(3U, c2_v2);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(_Control_linearMachineNumber_,
        chartInstance->chartNumber,chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization()
{
  return "4A0XmNvbvTlRXQdasYK3IH";
}

static void sf_check_dwork_consistency(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    const uint32_T *sfunDWorkChecksum = sf_get_sfun_dwork_checksum();
    mxArray *infoStruct = load_Control_linear_optimization_info();
    mxArray* mxRTWDWorkChecksum = sf_get_dwork_info_from_mat_file(S,
      sf_get_instance_specialization(), infoStruct, 2, "dworkChecksum");
    if (mxRTWDWorkChecksum != NULL) {
      double *pr = mxGetPr(mxRTWDWorkChecksum);
      if ((uint32_T)pr[0] != sfunDWorkChecksum[0] ||
          (uint32_T)pr[1] != sfunDWorkChecksum[1] ||
          (uint32_T)pr[2] != sfunDWorkChecksum[2] ||
          (uint32_T)pr[3] != sfunDWorkChecksum[3]) {
        sf_mex_error_message("Code generation and simulation targets registered different sets of persistent variables for the block. "
                             "External or Rapid Accelerator mode simulation requires code generation and simulation targets to "
                             "register the same set of persistent variables for this block. "
                             "This discrepancy is typically caused by MATLAB functions that have different code paths for "
                             "simulation and code generation targets where these code paths define different sets of persistent variables. "
                             "Please identify these code paths in the offending block and rewrite the MATLAB code so that "
                             "the set of persistent variables is the same between simulation and code generation.");
      }
    }
  }
}

static void sf_opaque_initialize_c2_Control_linear(void *chartInstanceVar)
{
  sf_check_dwork_consistency(((SFc2_Control_linearInstanceStruct*)
    chartInstanceVar)->S);
  chart_debug_initialization(((SFc2_Control_linearInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c2_Control_linear((SFc2_Control_linearInstanceStruct*)
    chartInstanceVar);
  initialize_c2_Control_linear((SFc2_Control_linearInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c2_Control_linear(void *chartInstanceVar)
{
  enable_c2_Control_linear((SFc2_Control_linearInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c2_Control_linear(void *chartInstanceVar)
{
  disable_c2_Control_linear((SFc2_Control_linearInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c2_Control_linear(void *chartInstanceVar)
{
  sf_c2_Control_linear((SFc2_Control_linearInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_ext_mode_exec_c2_Control_linear(void *chartInstanceVar)
{
  ext_mode_exec_c2_Control_linear((SFc2_Control_linearInstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c2_Control_linear(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c2_Control_linear
    ((SFc2_Control_linearInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_Control_linear();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c2_Control_linear(SimStruct* S, const
  mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c2_Control_linear();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c2_Control_linear((SFc2_Control_linearInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c2_Control_linear(SimStruct* S)
{
  return sf_internal_get_sim_state_c2_Control_linear(S);
}

static void sf_opaque_set_sim_state_c2_Control_linear(SimStruct* S, const
  mxArray *st)
{
  sf_internal_set_sim_state_c2_Control_linear(S, st);
}

static void sf_opaque_terminate_c2_Control_linear(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc2_Control_linearInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
    }

    finalize_c2_Control_linear((SFc2_Control_linearInstanceStruct*)
      chartInstanceVar);
    free((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }

  unload_Control_linear_optimization_info();
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc2_Control_linear((SFc2_Control_linearInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c2_Control_linear(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c2_Control_linear((SFc2_Control_linearInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

mxArray *sf_c2_Control_linear_get_testpoint_info(void)
{
  const char *infoEncStr[] = {
    "100 S'varName','path'{{T\"is_active_c2_Control_linear\",T\"is_active_c2_Control_linear\"}}"
  };

  mxArray *mxTpInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 1, 10);
  return mxTpInfo;
}

static void sf_set_sfun_dwork_info(SimStruct *S)
{
  const char *dworkEncStr[] = {
    "100 S1x4'type','isSigned','wordLength','bias','slope','exponent','isComplex','size'{{T\"int32\",,,,,,M[0],M[]},{T\"boolean\",,,,,,M[0],M[]},{T\"boolean\",,,,,,M[0],M[]},{T\"uint8\",,,,,,M[0],M[]}}"
  };

  sf_set_encoded_dwork_info(S, dworkEncStr, 4, 10);
}

static uint32_T* sf_get_sfun_dwork_checksum()
{
  static uint32_T checksum[4] = { 3851270630U, 3363230343U, 1651207761U,
    946165807U };

  return checksum;
}

static void mdlSetWorkWidths_c2_Control_linear(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Control_linear_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      2);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,2,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,2,
      "gatewayCannotBeInlinedMultipleTimes"));
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,2,2);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,2,2);
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,2);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
    sf_set_sfun_dwork_info(S);
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(843013872U));
  ssSetChecksum1(S,(1682906029U));
  ssSetChecksum2(S,(1367129730U));
  ssSetChecksum3(S,(1873113020U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c2_Control_linear(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c2_Control_linear(SimStruct *S)
{
  SFc2_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc2_Control_linearInstanceStruct *)malloc(sizeof
    (SFc2_Control_linearInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc2_Control_linearInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c2_Control_linear;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c2_Control_linear;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c2_Control_linear;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c2_Control_linear;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c2_Control_linear;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c2_Control_linear;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c2_Control_linear;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c2_Control_linear;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c2_Control_linear;
  chartInstance->chartInfo.mdlStart = mdlStart_c2_Control_linear;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c2_Control_linear;
  chartInstance->chartInfo.extModeExec =
    sf_opaque_ext_mode_exec_c2_Control_linear;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c2_Control_linear_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c2_Control_linear(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c2_Control_linear(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c2_Control_linear(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c2_Control_linear_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
