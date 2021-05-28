/*
 * Control_linear_data.c
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
#include "Control_linear.h"
#include "Control_linear_private.h"

/* Block parameters (auto storage) */
Parameters_Control_linear Control_linear_P = {
  1.0,                                 /* Expression: openloop
                                        * Referenced by: '<Root>/GT400-SV Initialization1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S7>/Transport Delay1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S7>/Transport Delay1'
                                        */
  12.0,                                /* Expression: 12
                                        * Referenced by: '<S7>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S7>/Transport Delay'
                                        */
  0.0046542113386515453,               /* Expression: (4*pi/180)/15
                                        * Referenced by: '<S7>/radtodeg1'
                                        */
  0.069813170079773182,                /* Expression: (4*pi/180)
                                        * Referenced by: '<S7>/Saturation'
                                        */
  -0.069813170079773182,               /* Expression: -(4*pi/180)
                                        * Referenced by: '<S7>/Saturation'
                                        */
  2.0,                                 /* Expression: axis
                                        * Referenced by: '<S10>/Get Current Axis' Position'
                                        */
  0.0026179938779914941,               /* Expression: 2*pi / 2400
                                        * Referenced by: '<S10>/Gain'
                                        */
  -0.545,                              /* Expression: -0.545
                                        * Referenced by: '<Root>/Constant7'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<Root>/radtodeg3'
                                        */
  0.0046542113386515453,               /* Expression: (4*pi/180)/15
                                        * Referenced by: '<S7>/Constant'
                                        */
  12.0,                                /* Expression: 12
                                        * Referenced by: '<S7>/Step1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S7>/Step1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S7>/Step1'
                                        */

  /*  Computed Parameter: TransferFcn_A
   * Referenced by: '<S5>/Transfer Fcn'
   */
  { -180.0, -10000.0 },

  /*  Computed Parameter: TransferFcn_C
   * Referenced by: '<S5>/Transfer Fcn'
   */
  { 10000.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay1'
                                        */
  12.0,                                /* Expression: 12
                                        * Referenced by: '<S9>/Transport Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay'
                                        */
  0.0046542113386515453,               /* Expression: (4*pi/180)/15
                                        * Referenced by: '<S9>/radtodeg1'
                                        */
  0.069813170079773182,                /* Expression: (4*pi/180)
                                        * Referenced by: '<S9>/Saturation'
                                        */
  -0.069813170079773182,               /* Expression: -(4*pi/180)
                                        * Referenced by: '<S9>/Saturation'
                                        */
  3.0,                                 /* Expression: axis
                                        * Referenced by: '<S10>/Get Current Axis' Position2'
                                        */
  -0.001308996938995747,               /* Expression: -0.5 * 2*pi / 2400
                                        * Referenced by: '<S10>/Gain2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<Root>/radtodeg7'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay2'
                                        */
  0.0046542113386515453,               /* Expression: (4*pi/180)/15
                                        * Referenced by: '<S9>/Constant'
                                        */
  12.0,                                /* Expression: 12
                                        * Referenced by: '<S9>/Step1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Step1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S9>/Step1'
                                        */

  /*  Computed Parameter: TransferFcn_A_b
   * Referenced by: '<S4>/Transfer Fcn'
   */
  { -180.0, -10000.0 },

  /*  Computed Parameter: TransferFcn_C_g
   * Referenced by: '<S4>/Transfer Fcn'
   */
  { 10000.0, 0.0 },

  /*  Expression: K_psi2theta
   * Referenced by: '<S4>/k'
   */
  { 1.8, 0.2 },
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S4>/radtodeg3'
                                        */
  1.0,                                 /* Expression: axis
                                        * Referenced by: '<S10>/Get Current Axis' Position1'
                                        */
  -0.0026179938779914941,              /* Expression: -2*pi / 2400
                                        * Referenced by: '<S10>/Gain1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<Root>/radtodeg1'
                                        */

  /*  Computed Parameter: TransferFcn_A_j
   * Referenced by: '<Root>/Transfer Fcn'
   */
  { -180.0, -10000.0 },

  /*  Computed Parameter: TransferFcn_C_o
   * Referenced by: '<Root>/Transfer Fcn'
   */
  { 10000.0, 0.0 },

  /*  Computed Parameter: TransferFcn_A_m
   * Referenced by: '<S6>/Transfer Fcn'
   */
  { -180.0, -10000.0 },

  /*  Computed Parameter: TransferFcn_C_f
   * Referenced by: '<S6>/Transfer Fcn'
   */
  { 10000.0, 0.0 },

  /*  Expression: K
   * Referenced by: '<Root>/radtodeg2'
   */
  { 12.929529681504389, -1.5717798806496006E-15, 7.1156002229318069,
    -1.3445136430558061E-15, -3.04875118463368E-16, 15.758065881561485,
    -1.3396696398383118E-16, 2.5662129581345443 },
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<Root>/radtodeg5'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<Root>/radtodeg4'
                                        */
  0.2634,                              /* Expression: a1
                                        * Referenced by: '<Root>/Constant1'
                                        */
  1.6279,                              /* Expression: a2
                                        * Referenced by: '<Root>/Constant2'
                                        */
  0.531,                               /* Expression: b1
                                        * Referenced by: '<Root>/Constant3'
                                        */
  25.6488,                             /* Expression: a3
                                        * Referenced by: '<Root>/Constant4'
                                        */
  5.3292,                              /* Expression: b2
                                        * Referenced by: '<Root>/Constant5'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<Root>/Gain'
                                        */
  4.0,                                 /* Expression: 4
                                        * Referenced by: '<Root>/Limit1'
                                        */
  0.05,                                /* Expression: 0.05
                                        * Referenced by: '<Root>/Limit1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<Root>/Gain1'
                                        */
  4.0,                                 /* Expression: 4
                                        * Referenced by: '<Root>/Limit2'
                                        */
  0.05,                                /* Expression: 0.05
                                        * Referenced by: '<Root>/Limit2'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<Root>/radtodeg6'
                                        */
  3276.7,                              /* Expression: 32767/10
                                        * Referenced by: '<S1>/Gain2'
                                        */
  15000.0,                             /* Expression: 15000
                                        * Referenced by: '<S1>/Limit1'
                                        */
  -15000.0,                            /* Expression: -15000
                                        * Referenced by: '<S1>/Limit1'
                                        */
  3276.7,                              /* Expression: 32767/10
                                        * Referenced by: '<S1>/Gain1'
                                        */
  15000.0,                             /* Expression: 15000
                                        * Referenced by: '<S1>/Limit2'
                                        */
  -15000.0,                            /* Expression: -15000
                                        * Referenced by: '<S1>/Limit2'
                                        */
  -60.0,                               /* Computed Parameter: TransferFcn2_A
                                        * Referenced by: '<S10>/Transfer Fcn2'
                                        */
  60.0,                                /* Computed Parameter: TransferFcn2_C
                                        * Referenced by: '<S10>/Transfer Fcn2'
                                        */
  1.0,                                 /* Expression: axis
                                        * Referenced by: '<S10>/Set Current Axis' Command1'
                                        */
  -60.0,                               /* Computed Parameter: TransferFcn1_A
                                        * Referenced by: '<S10>/Transfer Fcn1'
                                        */
  60.0,                                /* Computed Parameter: TransferFcn1_C
                                        * Referenced by: '<S10>/Transfer Fcn1'
                                        */
  2.0,                                 /* Expression: axis
                                        * Referenced by: '<S10>/Set Current Axis' Command2'
                                        */
  -0.545,                              /* Expression: -0.545
                                        * Referenced by: '<S5>/Constant7'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S7>/Constant1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S9>/Constant1'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch5_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch5'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch1_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch1'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch4_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch4'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch_CurrentSetting_e
                                        * Referenced by: '<S9>/Manual Switch'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch5_CurrentSetting_d
                                        * Referenced by: '<S9>/Manual Switch5'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch1_CurrentSetting_a
                                        * Referenced by: '<S9>/Manual Switch1'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch4_CurrentSetting_k
                                        * Referenced by: '<S9>/Manual Switch4'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch2_CurrentSetting
                                        * Referenced by: '<S10>/Manual Switch2'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch1_CurrentSetting_f
                                        * Referenced by: '<S10>/Manual Switch1'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch2_CurrentSetting_l
                                        * Referenced by: '<S7>/Manual Switch2'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch3_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch3'
                                        */
  1U,                                  /* Computed Parameter: ManualSwitch2_CurrentSetting_m
                                        * Referenced by: '<S9>/Manual Switch2'
                                        */
  1U                                   /* Computed Parameter: ManualSwitch3_CurrentSetting_d
                                        * Referenced by: '<S9>/Manual Switch3'
                                        */
};
