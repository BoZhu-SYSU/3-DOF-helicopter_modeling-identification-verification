/*
 * Control_linear_dt.h
 *
 * Code generation for model "Control_linear".
 *
 * Model version              : 1.1019
 * Simulink Coder version : 8.3 (R2012b) 20-Jul-2012
 * C source code generated on : Wed Aug 17 19:48:18 2016
 *
 * Target selection: rtwin.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "ext_types.h"

/* data type size table */
static uint_T rtDataTypeSizes[] = {
  sizeof(real_T),
  sizeof(real32_T),
  sizeof(int8_T),
  sizeof(uint8_T),
  sizeof(int16_T),
  sizeof(uint16_T),
  sizeof(int32_T),
  sizeof(uint32_T),
  sizeof(boolean_T),
  sizeof(fcn_call_T),
  sizeof(int_T),
  sizeof(pointer_T),
  sizeof(action_T),
  2*sizeof(uint32_T)
};

/* data type name table */
static const char_T * rtDataTypeNames[] = {
  "real_T",
  "real32_T",
  "int8_T",
  "uint8_T",
  "int16_T",
  "uint16_T",
  "int32_T",
  "uint32_T",
  "boolean_T",
  "fcn_call_T",
  "int_T",
  "pointer_T",
  "action_T",
  "timer_uint32_pair_T"
};

/* data type transitions for block I/O structure */
static DataTypeTransition rtBTransitions[] = {
  { (char_T *)(&Control_linear_B.Clock), 0, 0, 55 }
  ,

  { (char_T *)(&Control_linear_DWork.TransportDelay1_RWORK.modelTStart), 0, 0, 5
  },

  { (char_T *)(&Control_linear_DWork.TransportDelay1_PWORK.TUbufferPtrs[0]), 11,
    0, 56 },

  { (char_T *)(&Control_linear_DWork.sfEvent), 6, 0, 2 },

  { (char_T *)(&Control_linear_DWork.TransportDelay1_IWORK.Tail), 10, 0, 11 },

  { (char_T *)(&Control_linear_DWork.is_active_c2_Control_linear), 3, 0, 2 },

  { (char_T *)(&Control_linear_DWork.isStable), 8, 0, 4 }
};

/* data type transition table for block I/O structure */
static DataTypeTransitionTable rtBTransTable = {
  7U,
  rtBTransitions
};

/* data type transitions for Parameters structure */
static DataTypeTransition rtPTransitions[] = {
  { (char_T *)(&Control_linear_P.GT400SVInitialization1_P1), 0, 0, 91 },

  { (char_T *)(&Control_linear_P.ManualSwitch_CurrentSetting), 3, 0, 14 }
};

/* data type transition table for Parameters structure */
static DataTypeTransitionTable rtPTransTable = {
  2U,
  rtPTransitions
};

/* [EOF] Control_linear_dt.h */
