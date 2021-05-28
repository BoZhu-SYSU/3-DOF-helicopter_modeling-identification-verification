/* Include files */

#include "blascompat32.h"
#include "Control_linear_sfun.h"
#include "c1_Control_linear.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "Control_linear_sfun_debug_macros.h"

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c1_debug_family_names[13] = { "theta", "epsilon", "nargin",
  "nargout", "ele_d", "theta_d", "a1", "a2", "b1", "a3", "b2", "fs_d", "fd_d" };

/* Function Declarations */
static void initialize_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance);
static void initialize_params_c1_Control_linear
  (SFc1_Control_linearInstanceStruct *chartInstance);
static void enable_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance);
static void disable_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance);
static void c1_update_debugger_state_c1_Control_linear
  (SFc1_Control_linearInstanceStruct *chartInstance);
static void ext_mode_exec_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance);
static const mxArray *get_sim_state_c1_Control_linear
  (SFc1_Control_linearInstanceStruct *chartInstance);
static void set_sim_state_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_st);
static void finalize_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance);
static void sf_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance);
static void initSimStructsc1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber);
static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData);
static real_T c1_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_fd_d, const char_T *c1_identifier);
static real_T c1_b_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData);
static int32_T c1_c_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData);
static uint8_T c1_d_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_is_active_c1_Control_linear, const char_T
  *c1_identifier);
static uint8_T c1_e_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId);
static void init_dsm_address_info(SFc1_Control_linearInstanceStruct
  *chartInstance);

/* Function Definitions */
static void initialize_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance)
{
  int32_T *c1_sfEvent;
  uint8_T *c1_is_active_c1_Control_linear;
  c1_is_active_c1_Control_linear = (uint8_T *)ssGetDWork(chartInstance->S, 3);
  c1_sfEvent = (int32_T *)ssGetDWork(chartInstance->S, 0);
  *c1_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  *c1_is_active_c1_Control_linear = 0U;
}

static void initialize_params_c1_Control_linear
  (SFc1_Control_linearInstanceStruct *chartInstance)
{
}

static void enable_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c1_update_debugger_state_c1_Control_linear
  (SFc1_Control_linearInstanceStruct *chartInstance)
{
}

static void ext_mode_exec_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance)
{
  c1_update_debugger_state_c1_Control_linear(chartInstance);
}

static const mxArray *get_sim_state_c1_Control_linear
  (SFc1_Control_linearInstanceStruct *chartInstance)
{
  const mxArray *c1_st;
  const mxArray *c1_y = NULL;
  real_T c1_hoistedGlobal;
  real_T c1_u;
  const mxArray *c1_b_y = NULL;
  real_T c1_b_hoistedGlobal;
  real_T c1_b_u;
  const mxArray *c1_c_y = NULL;
  uint8_T c1_c_hoistedGlobal;
  uint8_T c1_c_u;
  const mxArray *c1_d_y = NULL;
  real_T *c1_fd_d;
  real_T *c1_fs_d;
  uint8_T *c1_is_active_c1_Control_linear;
  c1_fd_d = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_fs_d = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c1_is_active_c1_Control_linear = (uint8_T *)ssGetDWork(chartInstance->S, 3);
  c1_st = NULL;
  c1_st = NULL;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_createcellarray(3), FALSE);
  c1_hoistedGlobal = *c1_fd_d;
  c1_u = c1_hoistedGlobal;
  c1_b_y = NULL;
  sf_mex_assign(&c1_b_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 0, c1_b_y);
  c1_b_hoistedGlobal = *c1_fs_d;
  c1_b_u = c1_b_hoistedGlobal;
  c1_c_y = NULL;
  sf_mex_assign(&c1_c_y, sf_mex_create("y", &c1_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 1, c1_c_y);
  c1_c_hoistedGlobal = *c1_is_active_c1_Control_linear;
  c1_c_u = c1_c_hoistedGlobal;
  c1_d_y = NULL;
  sf_mex_assign(&c1_d_y, sf_mex_create("y", &c1_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c1_y, 2, c1_d_y);
  sf_mex_assign(&c1_st, c1_y, FALSE);
  return c1_st;
}

static void set_sim_state_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_st)
{
  const mxArray *c1_u;
  boolean_T *c1_doneDoubleBufferReInit;
  real_T *c1_fd_d;
  real_T *c1_fs_d;
  uint8_T *c1_is_active_c1_Control_linear;
  c1_fd_d = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_fs_d = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c1_is_active_c1_Control_linear = (uint8_T *)ssGetDWork(chartInstance->S, 3);
  c1_doneDoubleBufferReInit = (boolean_T *)ssGetDWork(chartInstance->S, 2);
  *c1_doneDoubleBufferReInit = TRUE;
  c1_u = sf_mex_dup(c1_st);
  *c1_fd_d = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u,
    0)), "fd_d");
  *c1_fs_d = c1_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c1_u,
    1)), "fs_d");
  *c1_is_active_c1_Control_linear = c1_d_emlrt_marshallIn(chartInstance,
    sf_mex_dup(sf_mex_getcell(c1_u, 2)), "is_active_c1_Control_linear");
  sf_mex_destroy(&c1_u);
  c1_update_debugger_state_c1_Control_linear(chartInstance);
  sf_mex_destroy(&c1_st);
}

static void finalize_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance)
{
}

static void sf_c1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance)
{
  real_T c1_hoistedGlobal;
  real_T c1_b_hoistedGlobal;
  real_T c1_c_hoistedGlobal;
  real_T c1_d_hoistedGlobal;
  real_T c1_e_hoistedGlobal;
  real_T c1_f_hoistedGlobal;
  real_T c1_g_hoistedGlobal;
  real_T c1_ele_d;
  real_T c1_theta_d;
  real_T c1_a1;
  real_T c1_a2;
  real_T c1_b1;
  real_T c1_a3;
  real_T c1_b2;
  uint32_T c1_debug_family_var_map[13];
  real_T c1_theta;
  real_T c1_epsilon;
  real_T c1_nargin = 7.0;
  real_T c1_nargout = 2.0;
  real_T c1_fs_d;
  real_T c1_fd_d;
  real_T c1_x;
  real_T c1_b_x;
  real_T c1_a;
  real_T c1_b;
  real_T c1_y;
  real_T c1_c_x;
  real_T c1_d_x;
  real_T c1_b_a;
  real_T c1_b_b;
  real_T c1_b_y;
  real_T c1_e_x;
  real_T c1_f_x;
  real_T c1_c_a;
  real_T c1_c_b;
  real_T c1_c_y;
  real_T c1_g_x;
  real_T c1_h_x;
  real_T c1_d_a;
  real_T c1_d_b;
  real_T c1_d_y;
  real_T c1_A;
  real_T c1_B;
  real_T c1_i_x;
  real_T c1_e_y;
  real_T c1_j_x;
  real_T c1_f_y;
  real_T c1_k_x;
  real_T c1_l_x;
  real_T c1_e_a;
  real_T c1_e_b;
  real_T c1_g_y;
  real_T c1_m_x;
  real_T c1_n_x;
  real_T c1_f_a;
  real_T c1_f_b;
  real_T c1_h_y;
  real_T c1_b_A;
  real_T c1_b_B;
  real_T c1_o_x;
  real_T c1_i_y;
  real_T c1_p_x;
  real_T c1_j_y;
  real_T *c1_b_fd_d;
  real_T *c1_b_fs_d;
  real_T *c1_b_b2;
  real_T *c1_b_a3;
  real_T *c1_b_b1;
  real_T *c1_b_a2;
  real_T *c1_b_a1;
  real_T *c1_b_theta_d;
  real_T *c1_b_ele_d;
  int32_T *c1_sfEvent;
  c1_b_b2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
  c1_b_a3 = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c1_b_b1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c1_b_a2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c1_b_a1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c1_b_theta_d = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c1_b_fd_d = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
  c1_b_ele_d = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c1_b_fs_d = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c1_sfEvent = (int32_T *)ssGetDWork(chartInstance->S, 0);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 0U, *c1_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c1_b_fs_d, 0U);
  _SFD_DATA_RANGE_CHECK(*c1_b_ele_d, 1U);
  _SFD_DATA_RANGE_CHECK(*c1_b_fd_d, 2U);
  _SFD_DATA_RANGE_CHECK(*c1_b_theta_d, 3U);
  _SFD_DATA_RANGE_CHECK(*c1_b_a1, 4U);
  _SFD_DATA_RANGE_CHECK(*c1_b_a2, 5U);
  _SFD_DATA_RANGE_CHECK(*c1_b_b1, 6U);
  _SFD_DATA_RANGE_CHECK(*c1_b_a3, 7U);
  _SFD_DATA_RANGE_CHECK(*c1_b_b2, 8U);
  *c1_sfEvent = CALL_EVENT;
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 0U, *c1_sfEvent);
  c1_hoistedGlobal = *c1_b_ele_d;
  c1_b_hoistedGlobal = *c1_b_theta_d;
  c1_c_hoistedGlobal = *c1_b_a1;
  c1_d_hoistedGlobal = *c1_b_a2;
  c1_e_hoistedGlobal = *c1_b_b1;
  c1_f_hoistedGlobal = *c1_b_a3;
  c1_g_hoistedGlobal = *c1_b_b2;
  c1_ele_d = c1_hoistedGlobal;
  c1_theta_d = c1_b_hoistedGlobal;
  c1_a1 = c1_c_hoistedGlobal;
  c1_a2 = c1_d_hoistedGlobal;
  c1_b1 = c1_e_hoistedGlobal;
  c1_a3 = c1_f_hoistedGlobal;
  c1_b2 = c1_g_hoistedGlobal;
  sf_debug_symbol_scope_push_eml(0U, 13U, 13U, c1_debug_family_names,
    c1_debug_family_var_map);
  sf_debug_symbol_scope_add_eml_importable(&c1_theta, 0U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  sf_debug_symbol_scope_add_eml_importable(&c1_epsilon, 1U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  sf_debug_symbol_scope_add_eml_importable(&c1_nargin, 2U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  sf_debug_symbol_scope_add_eml_importable(&c1_nargout, 3U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  sf_debug_symbol_scope_add_eml(&c1_ele_d, 4U, c1_sf_marshallOut);
  sf_debug_symbol_scope_add_eml(&c1_theta_d, 5U, c1_sf_marshallOut);
  sf_debug_symbol_scope_add_eml(&c1_a1, 6U, c1_sf_marshallOut);
  sf_debug_symbol_scope_add_eml(&c1_a2, 7U, c1_sf_marshallOut);
  sf_debug_symbol_scope_add_eml(&c1_b1, 8U, c1_sf_marshallOut);
  sf_debug_symbol_scope_add_eml(&c1_a3, 9U, c1_sf_marshallOut);
  sf_debug_symbol_scope_add_eml(&c1_b2, 10U, c1_sf_marshallOut);
  sf_debug_symbol_scope_add_eml_importable(&c1_fs_d, 11U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  sf_debug_symbol_scope_add_eml_importable(&c1_fd_d, 12U, c1_sf_marshallOut,
    c1_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, *c1_sfEvent, 2);
  c1_theta = c1_theta_d;
  _SFD_EML_CALL(0U, *c1_sfEvent, 3);
  c1_epsilon = c1_ele_d;
  _SFD_EML_CALL(0U, *c1_sfEvent, 5);
  c1_x = c1_theta;
  c1_b_x = c1_x;
  c1_b_x = muDoubleScalarCos(c1_b_x);
  c1_a = c1_a1;
  c1_b = c1_b_x;
  c1_y = c1_a * c1_b;
  c1_c_x = c1_epsilon;
  c1_d_x = c1_c_x;
  c1_d_x = muDoubleScalarSin(c1_d_x);
  c1_b_a = c1_y;
  c1_b_b = c1_d_x;
  c1_b_y = c1_b_a * c1_b_b;
  c1_e_x = c1_epsilon;
  c1_f_x = c1_e_x;
  c1_f_x = muDoubleScalarCos(c1_f_x);
  c1_c_a = c1_a2;
  c1_c_b = c1_f_x;
  c1_c_y = c1_c_a * c1_c_b;
  c1_g_x = c1_theta;
  c1_h_x = c1_g_x;
  c1_h_x = muDoubleScalarCos(c1_h_x);
  c1_d_a = c1_b1;
  c1_d_b = c1_h_x;
  c1_d_y = c1_d_a * c1_d_b;
  c1_A = c1_b_y + c1_c_y;
  c1_B = c1_d_y;
  c1_i_x = c1_A;
  c1_e_y = c1_B;
  c1_j_x = c1_i_x;
  c1_f_y = c1_e_y;
  c1_fs_d = c1_j_x / c1_f_y;
  _SFD_EML_CALL(0U, *c1_sfEvent, 6);
  c1_k_x = c1_epsilon;
  c1_l_x = c1_k_x;
  c1_l_x = muDoubleScalarCos(c1_l_x);
  c1_e_a = c1_a3;
  c1_e_b = c1_l_x;
  c1_g_y = c1_e_a * c1_e_b;
  c1_m_x = c1_theta;
  c1_n_x = c1_m_x;
  c1_n_x = muDoubleScalarSin(c1_n_x);
  c1_f_a = c1_g_y;
  c1_f_b = c1_n_x;
  c1_h_y = c1_f_a * c1_f_b;
  c1_b_A = c1_h_y;
  c1_b_B = c1_b2;
  c1_o_x = c1_b_A;
  c1_i_y = c1_b_B;
  c1_p_x = c1_o_x;
  c1_j_y = c1_i_y;
  c1_fd_d = c1_p_x / c1_j_y;
  _SFD_EML_CALL(0U, *c1_sfEvent, -6);
  sf_debug_symbol_scope_pop();
  *c1_b_fs_d = c1_fs_d;
  *c1_b_fd_d = c1_fd_d;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 0U, *c1_sfEvent);
  sf_debug_check_for_state_inconsistency(_Control_linearMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void initSimStructsc1_Control_linear(SFc1_Control_linearInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c1_machineNumber, uint32_T
  c1_chartNumber)
{
}

static const mxArray *c1_sf_marshallOut(void *chartInstanceVoid, void *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  real_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc1_Control_linearInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(real_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, FALSE);
  return c1_mxArrayOutData;
}

static real_T c1_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_fd_d, const char_T *c1_identifier)
{
  real_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_fd_d), &c1_thisId);
  sf_mex_destroy(&c1_fd_d);
  return c1_y;
}

static real_T c1_b_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  real_T c1_y;
  real_T c1_d0;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_d0, 1, 0, 0U, 0, 0U, 0);
  c1_y = c1_d0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_fd_d;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  real_T c1_y;
  SFc1_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc1_Control_linearInstanceStruct *)chartInstanceVoid;
  c1_fd_d = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_fd_d), &c1_thisId);
  sf_mex_destroy(&c1_fd_d);
  *(real_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

const mxArray *sf_c1_Control_linear_get_eml_resolved_functions_info(void)
{
  const mxArray *c1_nameCaptureInfo;
  c1_ResolvedFunctionInfo c1_info[8];
  c1_ResolvedFunctionInfo (*c1_b_info)[8];
  const mxArray *c1_m0 = NULL;
  int32_T c1_i0;
  c1_ResolvedFunctionInfo *c1_r0;
  c1_nameCaptureInfo = NULL;
  c1_nameCaptureInfo = NULL;
  c1_b_info = (c1_ResolvedFunctionInfo (*)[8])c1_info;
  (*c1_b_info)[0].context = "";
  (*c1_b_info)[0].name = "cos";
  (*c1_b_info)[0].dominantType = "double";
  (*c1_b_info)[0].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  (*c1_b_info)[0].fileTimeLo = 1286797106U;
  (*c1_b_info)[0].fileTimeHi = 0U;
  (*c1_b_info)[0].mFileTimeLo = 0U;
  (*c1_b_info)[0].mFileTimeHi = 0U;
  (*c1_b_info)[1].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  (*c1_b_info)[1].name = "eml_scalar_cos";
  (*c1_b_info)[1].dominantType = "double";
  (*c1_b_info)[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  (*c1_b_info)[1].fileTimeLo = 1286797122U;
  (*c1_b_info)[1].fileTimeHi = 0U;
  (*c1_b_info)[1].mFileTimeLo = 0U;
  (*c1_b_info)[1].mFileTimeHi = 0U;
  (*c1_b_info)[2].context = "";
  (*c1_b_info)[2].name = "mtimes";
  (*c1_b_info)[2].dominantType = "double";
  (*c1_b_info)[2].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  (*c1_b_info)[2].fileTimeLo = 1289494492U;
  (*c1_b_info)[2].fileTimeHi = 0U;
  (*c1_b_info)[2].mFileTimeLo = 0U;
  (*c1_b_info)[2].mFileTimeHi = 0U;
  (*c1_b_info)[3].context = "";
  (*c1_b_info)[3].name = "sin";
  (*c1_b_info)[3].dominantType = "double";
  (*c1_b_info)[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  (*c1_b_info)[3].fileTimeLo = 1286797150U;
  (*c1_b_info)[3].fileTimeHi = 0U;
  (*c1_b_info)[3].mFileTimeLo = 0U;
  (*c1_b_info)[3].mFileTimeHi = 0U;
  (*c1_b_info)[4].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  (*c1_b_info)[4].name = "eml_scalar_sin";
  (*c1_b_info)[4].dominantType = "double";
  (*c1_b_info)[4].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  (*c1_b_info)[4].fileTimeLo = 1286797136U;
  (*c1_b_info)[4].fileTimeHi = 0U;
  (*c1_b_info)[4].mFileTimeLo = 0U;
  (*c1_b_info)[4].mFileTimeHi = 0U;
  (*c1_b_info)[5].context = "";
  (*c1_b_info)[5].name = "mrdivide";
  (*c1_b_info)[5].dominantType = "double";
  (*c1_b_info)[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  (*c1_b_info)[5].fileTimeLo = 1342789344U;
  (*c1_b_info)[5].fileTimeHi = 0U;
  (*c1_b_info)[5].mFileTimeLo = 1319708366U;
  (*c1_b_info)[5].mFileTimeHi = 0U;
  (*c1_b_info)[6].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  (*c1_b_info)[6].name = "rdivide";
  (*c1_b_info)[6].dominantType = "double";
  (*c1_b_info)[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  (*c1_b_info)[6].fileTimeLo = 1286797244U;
  (*c1_b_info)[6].fileTimeHi = 0U;
  (*c1_b_info)[6].mFileTimeLo = 0U;
  (*c1_b_info)[6].mFileTimeHi = 0U;
  (*c1_b_info)[7].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  (*c1_b_info)[7].name = "eml_div";
  (*c1_b_info)[7].dominantType = "double";
  (*c1_b_info)[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  (*c1_b_info)[7].fileTimeLo = 1313326210U;
  (*c1_b_info)[7].fileTimeHi = 0U;
  (*c1_b_info)[7].mFileTimeLo = 0U;
  (*c1_b_info)[7].mFileTimeHi = 0U;
  sf_mex_assign(&c1_m0, sf_mex_createstruct("nameCaptureInfo", 1, 8), FALSE);
  for (c1_i0 = 0; c1_i0 < 8; c1_i0++) {
    c1_r0 = &c1_info[c1_i0];
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->context, 15,
      0U, 0U, 0U, 2, 1, strlen(c1_r0->context)), "context", "nameCaptureInfo",
                    c1_i0);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->name, 15, 0U,
      0U, 0U, 2, 1, strlen(c1_r0->name)), "name", "nameCaptureInfo", c1_i0);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->dominantType,
      15, 0U, 0U, 0U, 2, 1, strlen(c1_r0->dominantType)), "dominantType",
                    "nameCaptureInfo", c1_i0);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", c1_r0->resolved, 15,
      0U, 0U, 0U, 2, 1, strlen(c1_r0->resolved)), "resolved", "nameCaptureInfo",
                    c1_i0);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->fileTimeLo,
      7, 0U, 0U, 0U, 0), "fileTimeLo", "nameCaptureInfo", c1_i0);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->fileTimeHi,
      7, 0U, 0U, 0U, 0), "fileTimeHi", "nameCaptureInfo", c1_i0);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->mFileTimeLo,
      7, 0U, 0U, 0U, 0), "mFileTimeLo", "nameCaptureInfo", c1_i0);
    sf_mex_addfield(c1_m0, sf_mex_create("nameCaptureInfo", &c1_r0->mFileTimeHi,
      7, 0U, 0U, 0U, 0), "mFileTimeHi", "nameCaptureInfo", c1_i0);
  }

  sf_mex_assign(&c1_nameCaptureInfo, c1_m0, FALSE);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c1_nameCaptureInfo);
  return c1_nameCaptureInfo;
}

static const mxArray *c1_b_sf_marshallOut(void *chartInstanceVoid, void
  *c1_inData)
{
  const mxArray *c1_mxArrayOutData = NULL;
  int32_T c1_u;
  const mxArray *c1_y = NULL;
  SFc1_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc1_Control_linearInstanceStruct *)chartInstanceVoid;
  c1_mxArrayOutData = NULL;
  c1_u = *(int32_T *)c1_inData;
  c1_y = NULL;
  sf_mex_assign(&c1_y, sf_mex_create("y", &c1_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c1_mxArrayOutData, c1_y, FALSE);
  return c1_mxArrayOutData;
}

static int32_T c1_c_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  int32_T c1_y;
  int32_T c1_i1;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_i1, 1, 6, 0U, 0, 0U, 0);
  c1_y = c1_i1;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void c1_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c1_mxArrayInData, const char_T *c1_varName, void *c1_outData)
{
  const mxArray *c1_sfEvent;
  const char_T *c1_identifier;
  emlrtMsgIdentifier c1_thisId;
  int32_T c1_y;
  SFc1_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc1_Control_linearInstanceStruct *)chartInstanceVoid;
  c1_sfEvent = sf_mex_dup(c1_mxArrayInData);
  c1_identifier = c1_varName;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c1_sfEvent), &c1_thisId);
  sf_mex_destroy(&c1_sfEvent);
  *(int32_T *)c1_outData = c1_y;
  sf_mex_destroy(&c1_mxArrayInData);
}

static uint8_T c1_d_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_is_active_c1_Control_linear, const char_T
  *c1_identifier)
{
  uint8_T c1_y;
  emlrtMsgIdentifier c1_thisId;
  c1_thisId.fIdentifier = c1_identifier;
  c1_thisId.fParent = NULL;
  c1_y = c1_e_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c1_is_active_c1_Control_linear), &c1_thisId);
  sf_mex_destroy(&c1_is_active_c1_Control_linear);
  return c1_y;
}

static uint8_T c1_e_emlrt_marshallIn(SFc1_Control_linearInstanceStruct
  *chartInstance, const mxArray *c1_u, const emlrtMsgIdentifier *c1_parentId)
{
  uint8_T c1_y;
  uint8_T c1_u0;
  sf_mex_import(c1_parentId, sf_mex_dup(c1_u), &c1_u0, 1, 3, 0U, 0, 0U, 0);
  c1_y = c1_u0;
  sf_mex_destroy(&c1_u);
  return c1_y;
}

static void init_dsm_address_info(SFc1_Control_linearInstanceStruct
  *chartInstance)
{
}

/* SFunction Glue Code */
static uint32_T* sf_get_sfun_dwork_checksum();
void sf_c1_Control_linear_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3547028049U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(1411161142U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3635570390U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3630703071U);
}

mxArray *sf_c1_Control_linear_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("N7OrUlUXiL1BFg3aiY2VxF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,7,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,6,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,6,"type",mxType);
    }

    mxSetField(mxData,6,"complexity",mxCreateDoubleScalar(0));
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

static const mxArray *sf_get_sim_state_info_c1_Control_linear(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[6],T\"fd_d\",},{M[1],M[4],T\"fs_d\",},{M[8],M[0],T\"is_active_c1_Control_linear\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c1_Control_linear_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc1_Control_linearInstanceStruct *chartInstance;
    chartInstance = (SFc1_Control_linearInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (_Control_linearMachineNumber_,
           1,
           1,
           1,
           9,
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
          _SFD_SET_DATA_PROPS(0,2,0,1,"fs_d");
          _SFD_SET_DATA_PROPS(1,1,1,0,"ele_d");
          _SFD_SET_DATA_PROPS(2,2,0,1,"fd_d");
          _SFD_SET_DATA_PROPS(3,1,1,0,"theta_d");
          _SFD_SET_DATA_PROPS(4,1,1,0,"a1");
          _SFD_SET_DATA_PROPS(5,1,1,0,"a2");
          _SFD_SET_DATA_PROPS(6,1,1,0,"b1");
          _SFD_SET_DATA_PROPS(7,1,1,0,"a3");
          _SFD_SET_DATA_PROPS(8,1,1,0,"b2");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,199);
        _SFD_TRANS_COV_WTS(0,0,0,1,0);
        if (chartAlreadyPresent==0) {
          _SFD_TRANS_COV_MAPS(0,
                              0,NULL,NULL,
                              0,NULL,NULL,
                              1,NULL,NULL,
                              0,NULL,NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)c1_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c1_sf_marshallOut,(MexInFcnForType)NULL);

        {
          real_T *c1_fs_d;
          real_T *c1_ele_d;
          real_T *c1_fd_d;
          real_T *c1_theta_d;
          real_T *c1_a1;
          real_T *c1_a2;
          real_T *c1_b1;
          real_T *c1_a3;
          real_T *c1_b2;
          c1_b2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
          c1_a3 = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
          c1_b1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c1_a2 = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c1_a1 = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c1_theta_d = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c1_fd_d = (real_T *)ssGetOutputPortSignal(chartInstance->S, 2);
          c1_ele_d = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          c1_fs_d = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, c1_fs_d);
          _SFD_SET_DATA_VALUE_PTR(1U, c1_ele_d);
          _SFD_SET_DATA_VALUE_PTR(2U, c1_fd_d);
          _SFD_SET_DATA_VALUE_PTR(3U, c1_theta_d);
          _SFD_SET_DATA_VALUE_PTR(4U, c1_a1);
          _SFD_SET_DATA_VALUE_PTR(5U, c1_a2);
          _SFD_SET_DATA_VALUE_PTR(6U, c1_b1);
          _SFD_SET_DATA_VALUE_PTR(7U, c1_a3);
          _SFD_SET_DATA_VALUE_PTR(8U, c1_b2);
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
  return "UjjDYfo4lckhTGHVr7tDAF";
}

static void sf_check_dwork_consistency(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    const uint32_T *sfunDWorkChecksum = sf_get_sfun_dwork_checksum();
    mxArray *infoStruct = load_Control_linear_optimization_info();
    mxArray* mxRTWDWorkChecksum = sf_get_dwork_info_from_mat_file(S,
      sf_get_instance_specialization(), infoStruct, 1, "dworkChecksum");
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

static void sf_opaque_initialize_c1_Control_linear(void *chartInstanceVar)
{
  sf_check_dwork_consistency(((SFc1_Control_linearInstanceStruct*)
    chartInstanceVar)->S);
  chart_debug_initialization(((SFc1_Control_linearInstanceStruct*)
    chartInstanceVar)->S,0);
  initialize_params_c1_Control_linear((SFc1_Control_linearInstanceStruct*)
    chartInstanceVar);
  initialize_c1_Control_linear((SFc1_Control_linearInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_enable_c1_Control_linear(void *chartInstanceVar)
{
  enable_c1_Control_linear((SFc1_Control_linearInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c1_Control_linear(void *chartInstanceVar)
{
  disable_c1_Control_linear((SFc1_Control_linearInstanceStruct*)
    chartInstanceVar);
}

static void sf_opaque_gateway_c1_Control_linear(void *chartInstanceVar)
{
  sf_c1_Control_linear((SFc1_Control_linearInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_ext_mode_exec_c1_Control_linear(void *chartInstanceVar)
{
  ext_mode_exec_c1_Control_linear((SFc1_Control_linearInstanceStruct*)
    chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c1_Control_linear(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c1_Control_linear
    ((SFc1_Control_linearInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_Control_linear();/* state var info */
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

extern void sf_internal_set_sim_state_c1_Control_linear(SimStruct* S, const
  mxArray *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c1_Control_linear();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c1_Control_linear((SFc1_Control_linearInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c1_Control_linear(SimStruct* S)
{
  return sf_internal_get_sim_state_c1_Control_linear(S);
}

static void sf_opaque_set_sim_state_c1_Control_linear(SimStruct* S, const
  mxArray *st)
{
  sf_internal_set_sim_state_c1_Control_linear(S, st);
}

static void sf_opaque_terminate_c1_Control_linear(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc1_Control_linearInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
    }

    finalize_c1_Control_linear((SFc1_Control_linearInstanceStruct*)
      chartInstanceVar);
    free((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }

  unload_Control_linear_optimization_info();
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc1_Control_linear((SFc1_Control_linearInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c1_Control_linear(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c1_Control_linear((SFc1_Control_linearInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

mxArray *sf_c1_Control_linear_get_testpoint_info(void)
{
  const char *infoEncStr[] = {
    "100 S'varName','path'{{T\"is_active_c1_Control_linear\",T\"is_active_c1_Control_linear\"}}"
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

static void mdlSetWorkWidths_c1_Control_linear(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_Control_linear_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      1);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,1,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,1,
      "gatewayCannotBeInlinedMultipleTimes"));
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,1,7);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,1,2);
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,1);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
    sf_set_sfun_dwork_info(S);
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3728179446U));
  ssSetChecksum1(S,(3768152849U));
  ssSetChecksum2(S,(3598050014U));
  ssSetChecksum3(S,(671062562U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c1_Control_linear(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c1_Control_linear(SimStruct *S)
{
  SFc1_Control_linearInstanceStruct *chartInstance;
  chartInstance = (SFc1_Control_linearInstanceStruct *)malloc(sizeof
    (SFc1_Control_linearInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc1_Control_linearInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway =
    sf_opaque_gateway_c1_Control_linear;
  chartInstance->chartInfo.initializeChart =
    sf_opaque_initialize_c1_Control_linear;
  chartInstance->chartInfo.terminateChart =
    sf_opaque_terminate_c1_Control_linear;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c1_Control_linear;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c1_Control_linear;
  chartInstance->chartInfo.getSimState =
    sf_opaque_get_sim_state_c1_Control_linear;
  chartInstance->chartInfo.setSimState =
    sf_opaque_set_sim_state_c1_Control_linear;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c1_Control_linear;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c1_Control_linear;
  chartInstance->chartInfo.mdlStart = mdlStart_c1_Control_linear;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c1_Control_linear;
  chartInstance->chartInfo.extModeExec =
    sf_opaque_ext_mode_exec_c1_Control_linear;
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

void c1_Control_linear_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c1_Control_linear(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c1_Control_linear(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c1_Control_linear(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c1_Control_linear_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
