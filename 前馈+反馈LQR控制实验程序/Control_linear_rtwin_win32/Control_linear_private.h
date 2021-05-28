/*
 * Control_linear_private.h
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
#ifndef RTW_HEADER_Control_linear_private_h_
#define RTW_HEADER_Control_linear_private_h_
#include "rtwtypes.h"
#define SET_PRFL_VEL                   0x4
#define CLR_STATUS                     0x6
#define AXIS_ON                        0x9
#define AXIS_OFF                       0xa
#define OPEN_LOOP                      0xc
#define LMTS_ON                        0xd
#define LMTS_OFF                       0xe
#define ABRUPT_STOP                    0x12
#define RESET                          0x13
#define ZERO_POS                       0x15
#define UPDATE                         0x16
#define SET_LMT_SENSE                  0x2c
#define SET_ENCODER_SENSE              0x2d
#define SET_1                          0x30
#define SET_2                          0x31
#define SET_3                          0x32
#define SET_4                          0x33
#define GET_STATUS                     0x36
#define SET_POS                        0x41
#define SET_VEL                        0x42
#define SET_ACC                        0x43
#define SET_MTR_CMD                    0x49
#define SET_KP                         0x4c
#define SET_KI                         0x4d
#define SET_KD                         0x4e
#define GET_ACTL_POS                   0xc2
#define GET_ENC_POS                    0xc6
#define CHECK_NUM                      0x475410b5
#define DATA_PORT0                     0x0
#define COMMAND_PORT                   0x4
#define STATUS_PORT                    0x4
#define READY_BUSY_BIT_MASK            0x8000
#define PI                             3.1415927

uint_T baseaddress, icount, findcard_tag, positive_limit[4], negative_limit[4];
long positive_limit_value[4], negative_limit_value[4];

#define DELAY_TIME                     icount=0; while((!(inpd(baseaddress+STATUS_PORT)&READY_BUSY_BIT_MASK))&&(icount<6000)) icount++;
#define WritePort(command, data)       DELAY_TIME; outpd(baseaddress+DATA_PORT0,data); outpw(baseaddress+COMMAND_PORT,command); DELAY_TIME;
#define ReadPort(command, data)        DELAY_TIME; outpw(baseaddress+COMMAND_PORT,command); DELAY_TIME; data=inpd(baseaddress+DATA_PORT0);
#define ReadData(data)                 DELAY_TIME; data=inpd(baseaddress + DATA_PORT0);

/* Used by FromWorkspace Block: '<S11>/From Workspace' */
#ifndef rtInterpolate
# define rtInterpolate(v1,v2,f1,f2)    (((v1)==(v2))?((double)(v1)): (((f1)*((double)(v1)))+((f2)*((double)(v2)))))
#endif

#ifndef rtRound
# define rtRound(v)                    ( ((v) >= 0) ? floor((v) + 0.5) : ceil((v) - 0.5) )
#endif

#ifndef __RTWTYPES_H__
#error This file requires rtwtypes.h to be included
#else
#ifdef TMWTYPES_PREVIOUSLY_INCLUDED
#error This file requires rtwtypes.h to be included before tmwtypes.h
#else

/* Check for inclusion of an incorrect version of rtwtypes.h */
#ifndef RTWTYPES_ID_C08S16I32L32N32F1
#error This code was generated with a different "rtwtypes.h" than the file included
#endif                                 /* RTWTYPES_ID_C08S16I32L32N32F1 */
#endif                                 /* TMWTYPES_PREVIOUSLY_INCLUDED */
#endif                                 /* __RTWTYPES_H__ */

extern real_T rt_powd_snf(real_T u0, real_T u1);

/* Exported functions */
real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                  /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *tBuf,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
  ;

/* private model entry point functions */
extern void Control_linear_derivatives(void);

#endif                                 /* RTW_HEADER_Control_linear_private_h_ */
