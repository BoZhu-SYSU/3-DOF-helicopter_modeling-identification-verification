/*
 * Control_linear.h
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
#ifndef RTW_HEADER_Control_linear_h_
#define RTW_HEADER_Control_linear_h_
#ifndef Control_linear_COMMON_INCLUDES_
# define Control_linear_COMMON_INCLUDES_
#include <math.h>
#include <float.h>
#include <string.h>
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "dt_info.h"
#include "ext_work.h"
#include "rtwintgt.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"
#include "rtGetNaN.h"
#include "rt_defines.h"
#endif                                 /* Control_linear_COMMON_INCLUDES_ */

#include "Control_linear_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetBlkStateChangeFlag
# define rtmGetBlkStateChangeFlag(rtm) ((rtm)->ModelData.blkStateChange)
#endif

#ifndef rtmSetBlkStateChangeFlag
# define rtmSetBlkStateChangeFlag(rtm, val) ((rtm)->ModelData.blkStateChange = (val))
#endif

#ifndef rtmGetBlockIO
# define rtmGetBlockIO(rtm)            ((rtm)->ModelData.blockIO)
#endif

#ifndef rtmSetBlockIO
# define rtmSetBlockIO(rtm, val)       ((rtm)->ModelData.blockIO = (val))
#endif

#ifndef rtmGetChecksums
# define rtmGetChecksums(rtm)          ((rtm)->Sizes.checksums)
#endif

#ifndef rtmSetChecksums
# define rtmSetChecksums(rtm, val)     ((rtm)->Sizes.checksums = (val))
#endif

#ifndef rtmGetConstBlockIO
# define rtmGetConstBlockIO(rtm)       ((rtm)->ModelData.constBlockIO)
#endif

#ifndef rtmSetConstBlockIO
# define rtmSetConstBlockIO(rtm, val)  ((rtm)->ModelData.constBlockIO = (val))
#endif

#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->ModelData.contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->ModelData.contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->ModelData.contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->ModelData.contStates = (val))
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        ()
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)   ()
#endif

#ifndef rtmGetDefaultParam
# define rtmGetDefaultParam(rtm)       ((rtm)->ModelData.defaultParam)
#endif

#ifndef rtmSetDefaultParam
# define rtmSetDefaultParam(rtm, val)  ((rtm)->ModelData.defaultParam = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->ModelData.derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->ModelData.derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetDirectFeedThrough
# define rtmGetDirectFeedThrough(rtm)  ((rtm)->Sizes.sysDirFeedThru)
#endif

#ifndef rtmSetDirectFeedThrough
# define rtmSetDirectFeedThrough(rtm, val) ((rtm)->Sizes.sysDirFeedThru = (val))
#endif

#ifndef rtmGetErrorStatusFlag
# define rtmGetErrorStatusFlag(rtm)    ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatusFlag
# define rtmSetErrorStatusFlag(rtm, val) ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetFinalTime
# define rtmSetFinalTime(rtm, val)     ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetFirstInitCondFlag
# define rtmGetFirstInitCondFlag(rtm)  ()
#endif

#ifndef rtmSetFirstInitCondFlag
# define rtmSetFirstInitCondFlag(rtm, val) ()
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->ModelData.intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->ModelData.intgData = (val))
#endif

#ifndef rtmGetMdlRefGlobalTID
# define rtmGetMdlRefGlobalTID(rtm)    ()
#endif

#ifndef rtmSetMdlRefGlobalTID
# define rtmSetMdlRefGlobalTID(rtm, val) ()
#endif

#ifndef rtmGetMdlRefTriggerTID
# define rtmGetMdlRefTriggerTID(rtm)   ()
#endif

#ifndef rtmSetMdlRefTriggerTID
# define rtmSetMdlRefTriggerTID(rtm, val) ()
#endif

#ifndef rtmGetModelMappingInfo
# define rtmGetModelMappingInfo(rtm)   ((rtm)->SpecialInfo.mappingInfo)
#endif

#ifndef rtmSetModelMappingInfo
# define rtmSetModelMappingInfo(rtm, val) ((rtm)->SpecialInfo.mappingInfo = (val))
#endif

#ifndef rtmGetModelName
# define rtmGetModelName(rtm)          ((rtm)->modelName)
#endif

#ifndef rtmSetModelName
# define rtmSetModelName(rtm, val)     ((rtm)->modelName = (val))
#endif

#ifndef rtmGetNonInlinedSFcns
# define rtmGetNonInlinedSFcns(rtm)    ()
#endif

#ifndef rtmSetNonInlinedSFcns
# define rtmSetNonInlinedSFcns(rtm, val) ()
#endif

#ifndef rtmGetNumBlockIO
# define rtmGetNumBlockIO(rtm)         ((rtm)->Sizes.numBlockIO)
#endif

#ifndef rtmSetNumBlockIO
# define rtmSetNumBlockIO(rtm, val)    ((rtm)->Sizes.numBlockIO = (val))
#endif

#ifndef rtmGetNumBlockParams
# define rtmGetNumBlockParams(rtm)     ((rtm)->Sizes.numBlockPrms)
#endif

#ifndef rtmSetNumBlockParams
# define rtmSetNumBlockParams(rtm, val) ((rtm)->Sizes.numBlockPrms = (val))
#endif

#ifndef rtmGetNumBlocks
# define rtmGetNumBlocks(rtm)          ((rtm)->Sizes.numBlocks)
#endif

#ifndef rtmSetNumBlocks
# define rtmSetNumBlocks(rtm, val)     ((rtm)->Sizes.numBlocks = (val))
#endif

#ifndef rtmGetNumContStates
# define rtmGetNumContStates(rtm)      ((rtm)->Sizes.numContStates)
#endif

#ifndef rtmSetNumContStates
# define rtmSetNumContStates(rtm, val) ((rtm)->Sizes.numContStates = (val))
#endif

#ifndef rtmGetNumDWork
# define rtmGetNumDWork(rtm)           ((rtm)->Sizes.numDwork)
#endif

#ifndef rtmSetNumDWork
# define rtmSetNumDWork(rtm, val)      ((rtm)->Sizes.numDwork = (val))
#endif

#ifndef rtmGetNumInputPorts
# define rtmGetNumInputPorts(rtm)      ((rtm)->Sizes.numIports)
#endif

#ifndef rtmSetNumInputPorts
# define rtmSetNumInputPorts(rtm, val) ((rtm)->Sizes.numIports = (val))
#endif

#ifndef rtmGetNumNonSampledZCs
# define rtmGetNumNonSampledZCs(rtm)   ((rtm)->Sizes.numNonSampZCs)
#endif

#ifndef rtmSetNumNonSampledZCs
# define rtmSetNumNonSampledZCs(rtm, val) ((rtm)->Sizes.numNonSampZCs = (val))
#endif

#ifndef rtmGetNumOutputPorts
# define rtmGetNumOutputPorts(rtm)     ((rtm)->Sizes.numOports)
#endif

#ifndef rtmSetNumOutputPorts
# define rtmSetNumOutputPorts(rtm, val) ((rtm)->Sizes.numOports = (val))
#endif

#ifndef rtmGetNumSFcnParams
# define rtmGetNumSFcnParams(rtm)      ((rtm)->Sizes.numSFcnPrms)
#endif

#ifndef rtmSetNumSFcnParams
# define rtmSetNumSFcnParams(rtm, val) ((rtm)->Sizes.numSFcnPrms = (val))
#endif

#ifndef rtmGetNumSFunctions
# define rtmGetNumSFunctions(rtm)      ((rtm)->Sizes.numSFcns)
#endif

#ifndef rtmSetNumSFunctions
# define rtmSetNumSFunctions(rtm, val) ((rtm)->Sizes.numSFcns = (val))
#endif

#ifndef rtmGetNumSampleTimes
# define rtmGetNumSampleTimes(rtm)     ((rtm)->Sizes.numSampTimes)
#endif

#ifndef rtmSetNumSampleTimes
# define rtmSetNumSampleTimes(rtm, val) ((rtm)->Sizes.numSampTimes = (val))
#endif

#ifndef rtmGetNumU
# define rtmGetNumU(rtm)               ((rtm)->Sizes.numU)
#endif

#ifndef rtmSetNumU
# define rtmSetNumU(rtm, val)          ((rtm)->Sizes.numU = (val))
#endif

#ifndef rtmGetNumY
# define rtmGetNumY(rtm)               ((rtm)->Sizes.numY)
#endif

#ifndef rtmSetNumY
# define rtmSetNumY(rtm, val)          ((rtm)->Sizes.numY = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->ModelData.odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->ModelData.odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->ModelData.odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->ModelData.odeY = (val))
#endif

#ifndef rtmGetOffsetTimeArray
# define rtmGetOffsetTimeArray(rtm)    ((rtm)->Timing.offsetTimesArray)
#endif

#ifndef rtmSetOffsetTimeArray
# define rtmSetOffsetTimeArray(rtm, val) ((rtm)->Timing.offsetTimesArray = (val))
#endif

#ifndef rtmGetOffsetTimePtr
# define rtmGetOffsetTimePtr(rtm)      ((rtm)->Timing.offsetTimes)
#endif

#ifndef rtmSetOffsetTimePtr
# define rtmSetOffsetTimePtr(rtm, val) ((rtm)->Timing.offsetTimes = (val))
#endif

#ifndef rtmGetOptions
# define rtmGetOptions(rtm)            ((rtm)->Sizes.options)
#endif

#ifndef rtmSetOptions
# define rtmSetOptions(rtm, val)       ((rtm)->Sizes.options = (val))
#endif

#ifndef rtmGetParamIsMalloced
# define rtmGetParamIsMalloced(rtm)    ()
#endif

#ifndef rtmSetParamIsMalloced
# define rtmSetParamIsMalloced(rtm, val) ()
#endif

#ifndef rtmGetPath
# define rtmGetPath(rtm)               ((rtm)->path)
#endif

#ifndef rtmSetPath
# define rtmSetPath(rtm, val)          ((rtm)->path = (val))
#endif

#ifndef rtmGetPerTaskSampleHits
# define rtmGetPerTaskSampleHits(rtm)  ()
#endif

#ifndef rtmSetPerTaskSampleHits
# define rtmSetPerTaskSampleHits(rtm, val) ()
#endif

#ifndef rtmGetPerTaskSampleHitsArray
# define rtmGetPerTaskSampleHitsArray(rtm) ((rtm)->Timing.perTaskSampleHitsArray)
#endif

#ifndef rtmSetPerTaskSampleHitsArray
# define rtmSetPerTaskSampleHitsArray(rtm, val) ((rtm)->Timing.perTaskSampleHitsArray = (val))
#endif

#ifndef rtmGetPerTaskSampleHitsPtr
# define rtmGetPerTaskSampleHitsPtr(rtm) ((rtm)->Timing.perTaskSampleHits)
#endif

#ifndef rtmSetPerTaskSampleHitsPtr
# define rtmSetPerTaskSampleHitsPtr(rtm, val) ((rtm)->Timing.perTaskSampleHits = (val))
#endif

#ifndef rtmGetPrevZCSigState
# define rtmGetPrevZCSigState(rtm)     ((rtm)->ModelData.prevZCSigState)
#endif

#ifndef rtmSetPrevZCSigState
# define rtmSetPrevZCSigState(rtm, val) ((rtm)->ModelData.prevZCSigState = (val))
#endif

#ifndef rtmGetRTWExtModeInfo
# define rtmGetRTWExtModeInfo(rtm)     ((rtm)->extModeInfo)
#endif

#ifndef rtmSetRTWExtModeInfo
# define rtmSetRTWExtModeInfo(rtm, val) ((rtm)->extModeInfo = (val))
#endif

#ifndef rtmGetRTWGeneratedSFcn
# define rtmGetRTWGeneratedSFcn(rtm)   ((rtm)->Sizes.rtwGenSfcn)
#endif

#ifndef rtmSetRTWGeneratedSFcn
# define rtmSetRTWGeneratedSFcn(rtm, val) ((rtm)->Sizes.rtwGenSfcn = (val))
#endif

#ifndef rtmGetRTWLogInfo
# define rtmGetRTWLogInfo(rtm)         ()
#endif

#ifndef rtmSetRTWLogInfo
# define rtmSetRTWLogInfo(rtm, val)    ()
#endif

#ifndef rtmGetRTWRTModelMethodsInfo
# define rtmGetRTWRTModelMethodsInfo(rtm) ()
#endif

#ifndef rtmSetRTWRTModelMethodsInfo
# define rtmSetRTWRTModelMethodsInfo(rtm, val) ()
#endif

#ifndef rtmGetRTWSfcnInfo
# define rtmGetRTWSfcnInfo(rtm)        ((rtm)->sfcnInfo)
#endif

#ifndef rtmSetRTWSfcnInfo
# define rtmSetRTWSfcnInfo(rtm, val)   ((rtm)->sfcnInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfo
# define rtmGetRTWSolverInfo(rtm)      ((rtm)->solverInfo)
#endif

#ifndef rtmSetRTWSolverInfo
# define rtmSetRTWSolverInfo(rtm, val) ((rtm)->solverInfo = (val))
#endif

#ifndef rtmGetRTWSolverInfoPtr
# define rtmGetRTWSolverInfoPtr(rtm)   ((rtm)->solverInfoPtr)
#endif

#ifndef rtmSetRTWSolverInfoPtr
# define rtmSetRTWSolverInfoPtr(rtm, val) ((rtm)->solverInfoPtr = (val))
#endif

#ifndef rtmGetReservedForXPC
# define rtmGetReservedForXPC(rtm)     ((rtm)->SpecialInfo.xpcData)
#endif

#ifndef rtmSetReservedForXPC
# define rtmSetReservedForXPC(rtm, val) ((rtm)->SpecialInfo.xpcData = (val))
#endif

#ifndef rtmGetRootDWork
# define rtmGetRootDWork(rtm)          ((rtm)->Work.dwork)
#endif

#ifndef rtmSetRootDWork
# define rtmSetRootDWork(rtm, val)     ((rtm)->Work.dwork = (val))
#endif

#ifndef rtmGetSFunctions
# define rtmGetSFunctions(rtm)         ((rtm)->childSfunctions)
#endif

#ifndef rtmSetSFunctions
# define rtmSetSFunctions(rtm, val)    ((rtm)->childSfunctions = (val))
#endif

#ifndef rtmGetSampleHitArray
# define rtmGetSampleHitArray(rtm)     ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmSetSampleHitArray
# define rtmSetSampleHitArray(rtm, val) ((rtm)->Timing.sampleHitArray = (val))
#endif

#ifndef rtmGetSampleHitPtr
# define rtmGetSampleHitPtr(rtm)       ((rtm)->Timing.sampleHits)
#endif

#ifndef rtmSetSampleHitPtr
# define rtmSetSampleHitPtr(rtm, val)  ((rtm)->Timing.sampleHits = (val))
#endif

#ifndef rtmGetSampleTimeArray
# define rtmGetSampleTimeArray(rtm)    ((rtm)->Timing.sampleTimesArray)
#endif

#ifndef rtmSetSampleTimeArray
# define rtmSetSampleTimeArray(rtm, val) ((rtm)->Timing.sampleTimesArray = (val))
#endif

#ifndef rtmGetSampleTimePtr
# define rtmGetSampleTimePtr(rtm)      ((rtm)->Timing.sampleTimes)
#endif

#ifndef rtmSetSampleTimePtr
# define rtmSetSampleTimePtr(rtm, val) ((rtm)->Timing.sampleTimes = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDArray
# define rtmGetSampleTimeTaskIDArray(rtm) ((rtm)->Timing.sampleTimeTaskIDArray)
#endif

#ifndef rtmSetSampleTimeTaskIDArray
# define rtmSetSampleTimeTaskIDArray(rtm, val) ((rtm)->Timing.sampleTimeTaskIDArray = (val))
#endif

#ifndef rtmGetSampleTimeTaskIDPtr
# define rtmGetSampleTimeTaskIDPtr(rtm) ((rtm)->Timing.sampleTimeTaskIDPtr)
#endif

#ifndef rtmSetSampleTimeTaskIDPtr
# define rtmSetSampleTimeTaskIDPtr(rtm, val) ((rtm)->Timing.sampleTimeTaskIDPtr = (val))
#endif

#ifndef rtmGetSimMode
# define rtmGetSimMode(rtm)            ((rtm)->simMode)
#endif

#ifndef rtmSetSimMode
# define rtmSetSimMode(rtm, val)       ((rtm)->simMode = (val))
#endif

#ifndef rtmGetSimTimeStep
# define rtmGetSimTimeStep(rtm)        ((rtm)->Timing.simTimeStep)
#endif

#ifndef rtmSetSimTimeStep
# define rtmSetSimTimeStep(rtm, val)   ((rtm)->Timing.simTimeStep = (val))
#endif

#ifndef rtmGetStartTime
# define rtmGetStartTime(rtm)          ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetStartTime
# define rtmSetStartTime(rtm, val)     ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetStepSize
# define rtmGetStepSize(rtm)           ((rtm)->Timing.stepSize)
#endif

#ifndef rtmSetStepSize
# define rtmSetStepSize(rtm, val)      ((rtm)->Timing.stepSize = (val))
#endif

#ifndef rtmGetStopRequestedFlag
# define rtmGetStopRequestedFlag(rtm)  ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequestedFlag
# define rtmSetStopRequestedFlag(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetTaskCounters
# define rtmGetTaskCounters(rtm)       ()
#endif

#ifndef rtmSetTaskCounters
# define rtmSetTaskCounters(rtm, val)  ()
#endif

#ifndef rtmGetTaskTimeArray
# define rtmGetTaskTimeArray(rtm)      ((rtm)->Timing.tArray)
#endif

#ifndef rtmSetTaskTimeArray
# define rtmSetTaskTimeArray(rtm, val) ((rtm)->Timing.tArray = (val))
#endif

#ifndef rtmGetTimePtr
# define rtmGetTimePtr(rtm)            ((rtm)->Timing.t)
#endif

#ifndef rtmSetTimePtr
# define rtmSetTimePtr(rtm, val)       ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTimingData
# define rtmGetTimingData(rtm)         ((rtm)->Timing.timingData)
#endif

#ifndef rtmSetTimingData
# define rtmSetTimingData(rtm, val)    ((rtm)->Timing.timingData = (val))
#endif

#ifndef rtmGetU
# define rtmGetU(rtm)                  ((rtm)->ModelData.inputs)
#endif

#ifndef rtmSetU
# define rtmSetU(rtm, val)             ((rtm)->ModelData.inputs = (val))
#endif

#ifndef rtmGetVarNextHitTimesListPtr
# define rtmGetVarNextHitTimesListPtr(rtm) ((rtm)->Timing.varNextHitTimesList)
#endif

#ifndef rtmSetVarNextHitTimesListPtr
# define rtmSetVarNextHitTimesListPtr(rtm, val) ((rtm)->Timing.varNextHitTimesList = (val))
#endif

#ifndef rtmGetY
# define rtmGetY(rtm)                  ((rtm)->ModelData.outputs)
#endif

#ifndef rtmSetY
# define rtmSetY(rtm, val)             ((rtm)->ModelData.outputs = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->ModelData.zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->ModelData.zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetZCSignalValues
# define rtmGetZCSignalValues(rtm)     ((rtm)->ModelData.zcSignalValues)
#endif

#ifndef rtmSetZCSignalValues
# define rtmSetZCSignalValues(rtm, val) ((rtm)->ModelData.zcSignalValues = (val))
#endif

#ifndef rtmGet_TimeOfLastOutput
# define rtmGet_TimeOfLastOutput(rtm)  ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmSet_TimeOfLastOutput
# define rtmSet_TimeOfLastOutput(rtm, val) ((rtm)->Timing.timeOfLastOutput = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->ModelData.derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->ModelData.derivs = (val))
#endif

#ifndef rtmGetChecksumVal
# define rtmGetChecksumVal(rtm, idx)   ((rtm)->Sizes.checksums[idx])
#endif

#ifndef rtmSetChecksumVal
# define rtmSetChecksumVal(rtm, idx, val) ((rtm)->Sizes.checksums[idx] = (val))
#endif

#ifndef rtmGetDWork
# define rtmGetDWork(rtm, idx)         ((rtm)->Work.dwork[idx])
#endif

#ifndef rtmSetDWork
# define rtmSetDWork(rtm, idx, val)    ((rtm)->Work.dwork[idx] = (val))
#endif

#ifndef rtmGetOffsetTime
# define rtmGetOffsetTime(rtm, idx)    ((rtm)->Timing.offsetTimes[idx])
#endif

#ifndef rtmSetOffsetTime
# define rtmSetOffsetTime(rtm, idx, val) ((rtm)->Timing.offsetTimes[idx] = (val))
#endif

#ifndef rtmGetSFunction
# define rtmGetSFunction(rtm, idx)     ((rtm)->childSfunctions[idx])
#endif

#ifndef rtmSetSFunction
# define rtmSetSFunction(rtm, idx, val) ((rtm)->childSfunctions[idx] = (val))
#endif

#ifndef rtmGetSampleTime
# define rtmGetSampleTime(rtm, idx)    ((rtm)->Timing.sampleTimes[idx])
#endif

#ifndef rtmSetSampleTime
# define rtmSetSampleTime(rtm, idx, val) ((rtm)->Timing.sampleTimes[idx] = (val))
#endif

#ifndef rtmGetSampleTimeTaskID
# define rtmGetSampleTimeTaskID(rtm, idx) ((rtm)->Timing.sampleTimeTaskIDPtr[idx])
#endif

#ifndef rtmSetSampleTimeTaskID
# define rtmSetSampleTimeTaskID(rtm, idx, val) ((rtm)->Timing.sampleTimeTaskIDPtr[idx] = (val))
#endif

#ifndef rtmGetVarNextHitTimeList
# define rtmGetVarNextHitTimeList(rtm, idx) ((rtm)->Timing.varNextHitTimesList[idx])
#endif

#ifndef rtmSetVarNextHitTimeList
# define rtmSetVarNextHitTimeList(rtm, idx, val) ((rtm)->Timing.varNextHitTimesList[idx] = (val))
#endif

#ifndef rtmIsContinuousTask
# define rtmIsContinuousTask(rtm, tid) ((tid) == 0)
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmIsSampleHit
# define rtmIsSampleHit(rtm, sti, tid) ((rtmIsMajorTimeStep((rtm)) && (rtm)->Timing.sampleHits[(rtm)->Timing.sampleTimeTaskIDPtr[sti]]))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmSetT
# define rtmSetT(rtm, val)                                       /* Do Nothing */
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmSetTFinal
# define rtmSetTFinal(rtm, val)        ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

#ifndef rtmGetTStart
# define rtmGetTStart(rtm)             ((rtm)->Timing.tStart)
#endif

#ifndef rtmSetTStart
# define rtmSetTStart(rtm, val)        ((rtm)->Timing.tStart = (val))
#endif

#ifndef rtmGetTaskTime
# define rtmGetTaskTime(rtm, sti)      (rtmGetTPtr((rtm))[(rtm)->Timing.sampleTimeTaskIDPtr[sti]])
#endif

#ifndef rtmSetTaskTime
# define rtmSetTaskTime(rtm, sti, val) (rtmGetTPtr((rtm))[sti] = (val))
#endif

#ifndef rtmGetTimeOfLastOutput
# define rtmGetTimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

#ifdef rtmGetRTWSolverInfo
#undef rtmGetRTWSolverInfo
#endif

#define rtmGetRTWSolverInfo(rtm)       &((rtm)->solverInfo)

/* Definition for use in the target main file */
#define Control_linear_rtModel         RT_MODEL_Control_linear

/* Block signals (auto storage) */
typedef struct {
  real_T Clock;                        /* '<S7>/Clock' */
  real_T FromWorkspace;                /* '<S11>/From Workspace' */
  real_T ManualSwitch5;                /* '<S7>/Manual Switch5' */
  real_T GetCurrentAxisPosition;       /* '<S10>/Get Current Axis' Position' */
  real_T Gain;                         /* '<S10>/Gain' */
  real_T q1;                           /* '<Root>/Sum1' */
  real_T radtodeg3;                    /* '<Root>/radtodeg3' */
  real_T Sum2;                         /* '<S5>/Sum2' */
  real_T FromWorkspace_b;              /* '<S12>/From Workspace' */
  real_T ManualSwitch4;                /* '<S7>/Manual Switch4' */
  real_T TransferFcn;                  /* '<S5>/Transfer Fcn' */
  real_T Sum7;                         /* '<S5>/Sum7' */
  real_T Clock_e;                      /* '<S9>/Clock' */
  real_T FromWorkspace_f;              /* '<S17>/From Workspace' */
  real_T ManualSwitch5_k;              /* '<S9>/Manual Switch5' */
  real_T GetCurrentAxisPosition2;      /* '<S10>/Get Current Axis' Position2' */
  real_T Gain2;                        /* '<S10>/Gain2' */
  real_T radtodeg7;                    /* '<Root>/radtodeg7' */
  real_T Sum1;                         /* '<S4>/Sum1' */
  real_T Fcn1;                         /* '<S9>/Fcn1' */
  real_T FromWorkspace_bg;             /* '<S18>/From Workspace' */
  real_T ManualSwitch4_h;              /* '<S9>/Manual Switch4' */
  real_T TransferFcn_n;                /* '<S4>/Transfer Fcn' */
  real_T Sum6;                         /* '<S4>/Sum6' */
  real_T radtodeg3_l;                  /* '<S4>/radtodeg3' */
  real_T GetCurrentAxisPosition1;      /* '<S10>/Get Current Axis' Position1' */
  real_T Gain1;                        /* '<S10>/Gain1' */
  real_T radtodeg1;                    /* '<Root>/radtodeg1' */
  real_T TransferFcn_k;                /* '<Root>/Transfer Fcn' */
  real_T TransferFcn_l;                /* '<S6>/Transfer Fcn' */
  real_T Sum8;                         /* '<S6>/Sum8' */
  real_T radtodeg2[2];                 /* '<Root>/radtodeg2' */
  real_T radtodeg5;                    /* '<Root>/radtodeg5' */
  real_T radtodeg4;                    /* '<Root>/radtodeg4' */
  real_T Sum9;                         /* '<Root>/Sum9' */
  real_T Sum10;                        /* '<Root>/Sum10' */
  real_T Limit1;                       /* '<Root>/Limit1' */
  real_T Limit2;                       /* '<Root>/Limit2' */
  real_T radtodeg6;                    /* '<Root>/radtodeg6' */
  real_T Limit1_a;                     /* '<S1>/Limit1' */
  real_T Limit2_o;                     /* '<S1>/Limit2' */
  real_T ManualSwitch2;                /* '<S10>/Manual Switch2' */
  real_T ManualSwitch1;                /* '<S10>/Manual Switch1' */
  real_T Sum1_p;                       /* '<S5>/Sum1' */
  real_T FromWorkspace_e;              /* '<S13>/From Workspace' */
  real_T ManualSwitch3;                /* '<S7>/Manual Switch3' */
  real_T Clock1;                       /* '<S7>/Clock1' */
  real_T FromWorkspace_c;              /* '<S19>/From Workspace' */
  real_T ManualSwitch3_o;              /* '<S9>/Manual Switch3' */
  real_T Clock1_i;                     /* '<S9>/Clock1' */
  real_T v1;                           /* '<Root>/MATLAB Function4' */
  real_T v2;                           /* '<Root>/MATLAB Function4' */
  real_T fs_d;                         /* '<Root>/MATLAB Function1' */
  real_T fd_d;                         /* '<Root>/MATLAB Function1' */
} BlockIO_Control_linear;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay1_RWORK;             /* '<S7>/Transport Delay1' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK;              /* '<S7>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay1_RWORK_l;           /* '<S9>/Transport Delay1' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay_RWORK_e;            /* '<S9>/Transport Delay' */

  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } TransportDelay2_RWORK;             /* '<S9>/Transport Delay2' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay1_PWORK;             /* '<S7>/Transport Delay1' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK;              /* '<S7>/Transport Delay' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace_PWORK;               /* '<S11>/From Workspace' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace_PWORK_c;             /* '<S12>/From Workspace' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay1_PWORK_i;           /* '<S9>/Transport Delay1' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay_PWORK_o;            /* '<S9>/Transport Delay' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace_PWORK_j;             /* '<S17>/From Workspace' */

  struct {
    void *TUbufferPtrs[2];
  } TransportDelay2_PWORK;             /* '<S9>/Transport Delay2' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace_PWORK_p;             /* '<S18>/From Workspace' */

  struct {
    void *LoggedData;
  } Scope1_PWORK;                      /* '<Root>/Scope1' */

  struct {
    void *LoggedData;
  } Scope2_PWORK;                      /* '<Root>/Scope2' */

  struct {
    void *LoggedData;
  } Scope3_PWORK;                      /* '<Root>/Scope3' */

  struct {
    void *LoggedData;
  } Scope4_PWORK;                      /* '<Root>/Scope4' */

  struct {
    void *LoggedData;
  } ele1_PWORK;                        /* '<Root>/ele1' */

  struct {
    void *LoggedData;
  } ele2_PWORK;                        /* '<Root>/ele2' */

  struct {
    void *LoggedData;
  } pitch1_PWORK;                      /* '<Root>/pitch1' */

  struct {
    void *LoggedData;
  } pitch2_PWORK;                      /* '<Root>/pitch2' */

  struct {
    void *LoggedData;
  } scope1_PWORK;                      /* '<Root>/scope1' */

  struct {
    void *LoggedData;
  } scope11_PWORK;                     /* '<Root>/scope11' */

  struct {
    void *LoggedData;
  } scope18_PWORK;                     /* '<Root>/scope18' */

  struct {
    void *LoggedData;
  } scope2_PWORK;                      /* '<Root>/scope2' */

  struct {
    void *LoggedData;
  } scope6_PWORK;                      /* '<Root>/scope6' */

  struct {
    void *LoggedData;
  } travel1_PWORK;                     /* '<Root>/travel1' */

  struct {
    void *LoggedData;
  } travel2_PWORK;                     /* '<Root>/travel2' */

  struct {
    void *LoggedData;
  } TravelAngle1_PWORK;                /* '<S1>/Travel Angle1' */

  struct {
    void *LoggedData;
  } TravelAngle2_PWORK;                /* '<S1>/Travel Angle2' */

  struct {
    void *LoggedData;
  } TravelAngle1_PWORK_j;              /* '<S10>/Travel Angle1' */

  struct {
    void *LoggedData;
  } Scope1_PWORK_p;                    /* '<S4>/Scope1' */

  struct {
    void *LoggedData;
  } Scope11_PWORK;                     /* '<S4>/Scope11' */

  struct {
    void *LoggedData;
  } Scope12_PWORK;                     /* '<S4>/Scope12' */

  struct {
    void *LoggedData;
  } Scope13_PWORK;                     /* '<S4>/Scope13' */

  struct {
    void *LoggedData;
  } Scope2_PWORK_k;                    /* '<S4>/Scope2' */

  struct {
    void *LoggedData;
  } Scope3_PWORK_f;                    /* '<S4>/Scope3' */

  struct {
    void *LoggedData;
  } Scope6_PWORK;                      /* '<S4>/Scope6' */

  struct {
    void *LoggedData;
  } Scope1_PWORK_n;                    /* '<S5>/Scope1' */

  struct {
    void *LoggedData;
  } Scope4_PWORK_d;                    /* '<S5>/Scope4' */

  struct {
    void *LoggedData;
  } Scope7_PWORK;                      /* '<S5>/Scope7' */

  struct {
    void *LoggedData;
  } Scope9_PWORK;                      /* '<S5>/Scope9' */

  struct {
    void *LoggedData;
  } scope4_PWORK;                      /* '<S5>/scope4' */

  struct {
    void *LoggedData;
  } Scope10_PWORK;                     /* '<S6>/Scope10' */

  struct {
    void *LoggedData;
  } Scope5_PWORK;                      /* '<S6>/Scope5' */

  struct {
    void *LoggedData;
  } Scope8_PWORK;                      /* '<S6>/Scope8' */

  struct {
    void *LoggedData;
  } scope5_PWORK;                      /* '<S6>/scope5' */

  struct {
    void *LoggedData;
  } Scope1_PWORK_k;                    /* '<S7>/Scope1' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace_PWORK_k;             /* '<S13>/From Workspace' */

  struct {
    void *LoggedData;
  } Scope2_PWORK_f;                    /* '<S7>/Scope2' */

  struct {
    void *LoggedData;
  } scope11_PWORK_j;                   /* '<S7>/scope11' */

  struct {
    void *LoggedData;
  } Scope1_PWORK_d;                    /* '<S9>/Scope1' */

  struct {
    void *TimePtr;
    void *DataPtr;
    void *RSimInfoPtr;
  } FromWorkspace_PWORK_g;             /* '<S19>/From Workspace' */

  struct {
    void *LoggedData;
  } Scope2_PWORK_g;                    /* '<S9>/Scope2' */

  struct {
    void *LoggedData;
  } scope11_PWORK_a;                   /* '<S9>/scope11' */

  int32_T sfEvent;                     /* '<Root>/MATLAB Function4' */
  int32_T sfEvent_b;                   /* '<Root>/MATLAB Function1' */
  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay1_IWORK;             /* '<S7>/Transport Delay1' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK;              /* '<S7>/Transport Delay' */

  struct {
    int_T PrevIndex;
  } FromWorkspace_IWORK;               /* '<S11>/From Workspace' */

  struct {
    int_T PrevIndex;
  } FromWorkspace_IWORK_o;             /* '<S12>/From Workspace' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay1_IWORK_e;           /* '<S9>/Transport Delay1' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay_IWORK_c;            /* '<S9>/Transport Delay' */

  struct {
    int_T PrevIndex;
  } FromWorkspace_IWORK_i;             /* '<S17>/From Workspace' */

  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } TransportDelay2_IWORK;             /* '<S9>/Transport Delay2' */

  struct {
    int_T PrevIndex;
  } FromWorkspace_IWORK_n;             /* '<S18>/From Workspace' */

  struct {
    int_T PrevIndex;
  } FromWorkspace_IWORK_d;             /* '<S13>/From Workspace' */

  struct {
    int_T PrevIndex;
  } FromWorkspace_IWORK_nh;            /* '<S19>/From Workspace' */

  uint8_T is_active_c2_Control_linear; /* '<Root>/MATLAB Function4' */
  uint8_T is_active_c1_Control_linear; /* '<Root>/MATLAB Function1' */
  boolean_T isStable;                  /* '<Root>/MATLAB Function4' */
  boolean_T doneDoubleBufferReInit;    /* '<Root>/MATLAB Function4' */
  boolean_T isStable_j;                /* '<Root>/MATLAB Function1' */
  boolean_T doneDoubleBufferReInit_m;  /* '<Root>/MATLAB Function1' */
} D_Work_Control_linear;

/* Continuous states (auto storage) */
typedef struct {
  real_T TransferFcn_CSTATE[2];        /* '<S5>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_l[2];      /* '<S4>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_h[2];      /* '<Root>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_e[2];      /* '<S6>/Transfer Fcn' */
  real_T TransferFcn2_CSTATE;          /* '<S10>/Transfer Fcn2' */
  real_T TransferFcn1_CSTATE;          /* '<S10>/Transfer Fcn1' */
} ContinuousStates_Control_linear;

/* State derivatives (auto storage) */
typedef struct {
  real_T TransferFcn_CSTATE[2];        /* '<S5>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_l[2];      /* '<S4>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_h[2];      /* '<Root>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_e[2];      /* '<S6>/Transfer Fcn' */
  real_T TransferFcn2_CSTATE;          /* '<S10>/Transfer Fcn2' */
  real_T TransferFcn1_CSTATE;          /* '<S10>/Transfer Fcn1' */
} StateDerivatives_Control_linear;

/* State disabled  */
typedef struct {
  boolean_T TransferFcn_CSTATE[2];     /* '<S5>/Transfer Fcn' */
  boolean_T TransferFcn_CSTATE_l[2];   /* '<S4>/Transfer Fcn' */
  boolean_T TransferFcn_CSTATE_h[2];   /* '<Root>/Transfer Fcn' */
  boolean_T TransferFcn_CSTATE_e[2];   /* '<S6>/Transfer Fcn' */
  boolean_T TransferFcn2_CSTATE;       /* '<S10>/Transfer Fcn2' */
  boolean_T TransferFcn1_CSTATE;       /* '<S10>/Transfer Fcn1' */
} StateDisabled_Control_linear;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Backward compatible GRT Identifiers */
#define rtB                            Control_linear_B
#define BlockIO                        BlockIO_Control_linear
#define rtX                            Control_linear_X
#define ContinuousStates               ContinuousStates_Control_linear
#define rtXdot                         Control_linear_Xdot
#define StateDerivatives               StateDerivatives_Control_linear
#define tXdis                          Control_linear_Xdis
#define StateDisabled                  StateDisabled_Control_linear
#define rtP                            Control_linear_P
#define Parameters                     Parameters_Control_linear
#define rtDWork                        Control_linear_DWork
#define D_Work                         D_Work_Control_linear

/* Parameters (auto storage) */
struct Parameters_Control_linear_ {
  real_T GT400SVInitialization1_P1;    /* Expression: openloop
                                        * Referenced by: '<Root>/GT400-SV Initialization1'
                                        */
  real_T TransportDelay1_Delay;        /* Expression: 0
                                        * Referenced by: '<S7>/Transport Delay1'
                                        */
  real_T TransportDelay1_InitOutput;   /* Expression: 0
                                        * Referenced by: '<S7>/Transport Delay1'
                                        */
  real_T TransportDelay_Delay;         /* Expression: 12
                                        * Referenced by: '<S7>/Transport Delay'
                                        */
  real_T TransportDelay_InitOutput;    /* Expression: 0
                                        * Referenced by: '<S7>/Transport Delay'
                                        */
  real_T radtodeg1_Gain;               /* Expression: (4*pi/180)/15
                                        * Referenced by: '<S7>/radtodeg1'
                                        */
  real_T Saturation_UpperSat;          /* Expression: (4*pi/180)
                                        * Referenced by: '<S7>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: -(4*pi/180)
                                        * Referenced by: '<S7>/Saturation'
                                        */
  real_T GetCurrentAxisPosition_P1;    /* Expression: axis
                                        * Referenced by: '<S10>/Get Current Axis' Position'
                                        */
  real_T Gain_Gain;                    /* Expression: 2*pi / 2400
                                        * Referenced by: '<S10>/Gain'
                                        */
  real_T Constant7_Value;              /* Expression: -0.545
                                        * Referenced by: '<Root>/Constant7'
                                        */
  real_T radtodeg3_Gain;               /* Expression: 1
                                        * Referenced by: '<Root>/radtodeg3'
                                        */
  real_T Constant_Value;               /* Expression: (4*pi/180)/15
                                        * Referenced by: '<S7>/Constant'
                                        */
  real_T Step1_Time;                   /* Expression: 12
                                        * Referenced by: '<S7>/Step1'
                                        */
  real_T Step1_Y0;                     /* Expression: 0
                                        * Referenced by: '<S7>/Step1'
                                        */
  real_T Step1_YFinal;                 /* Expression: 1
                                        * Referenced by: '<S7>/Step1'
                                        */
  real_T TransferFcn_A[2];             /* Computed Parameter: TransferFcn_A
                                        * Referenced by: '<S5>/Transfer Fcn'
                                        */
  real_T TransferFcn_C[2];             /* Computed Parameter: TransferFcn_C
                                        * Referenced by: '<S5>/Transfer Fcn'
                                        */
  real_T TransportDelay1_Delay_m;      /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay1'
                                        */
  real_T TransportDelay1_InitOutput_g; /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay1'
                                        */
  real_T TransportDelay_Delay_b;       /* Expression: 12
                                        * Referenced by: '<S9>/Transport Delay'
                                        */
  real_T TransportDelay_InitOutput_m;  /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay'
                                        */
  real_T radtodeg1_Gain_b;             /* Expression: (4*pi/180)/15
                                        * Referenced by: '<S9>/radtodeg1'
                                        */
  real_T Saturation_UpperSat_b;        /* Expression: (4*pi/180)
                                        * Referenced by: '<S9>/Saturation'
                                        */
  real_T Saturation_LowerSat_m;        /* Expression: -(4*pi/180)
                                        * Referenced by: '<S9>/Saturation'
                                        */
  real_T GetCurrentAxisPosition2_P1;   /* Expression: axis
                                        * Referenced by: '<S10>/Get Current Axis' Position2'
                                        */
  real_T Gain2_Gain;                   /* Expression: -0.5 * 2*pi / 2400
                                        * Referenced by: '<S10>/Gain2'
                                        */
  real_T radtodeg7_Gain;               /* Expression: 1
                                        * Referenced by: '<Root>/radtodeg7'
                                        */
  real_T TransportDelay2_Delay;        /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay2'
                                        */
  real_T TransportDelay2_InitOutput;   /* Expression: 0
                                        * Referenced by: '<S9>/Transport Delay2'
                                        */
  real_T Constant_Value_n;             /* Expression: (4*pi/180)/15
                                        * Referenced by: '<S9>/Constant'
                                        */
  real_T Step1_Time_n;                 /* Expression: 12
                                        * Referenced by: '<S9>/Step1'
                                        */
  real_T Step1_Y0_i;                   /* Expression: 0
                                        * Referenced by: '<S9>/Step1'
                                        */
  real_T Step1_YFinal_f;               /* Expression: 1
                                        * Referenced by: '<S9>/Step1'
                                        */
  real_T TransferFcn_A_b[2];           /* Computed Parameter: TransferFcn_A_b
                                        * Referenced by: '<S4>/Transfer Fcn'
                                        */
  real_T TransferFcn_C_g[2];           /* Computed Parameter: TransferFcn_C_g
                                        * Referenced by: '<S4>/Transfer Fcn'
                                        */
  real_T k_Gain[2];                    /* Expression: K_psi2theta
                                        * Referenced by: '<S4>/k'
                                        */
  real_T radtodeg3_Gain_l;             /* Expression: -1
                                        * Referenced by: '<S4>/radtodeg3'
                                        */
  real_T GetCurrentAxisPosition1_P1;   /* Expression: axis
                                        * Referenced by: '<S10>/Get Current Axis' Position1'
                                        */
  real_T Gain1_Gain;                   /* Expression: -2*pi / 2400
                                        * Referenced by: '<S10>/Gain1'
                                        */
  real_T radtodeg1_Gain_f;             /* Expression: 1
                                        * Referenced by: '<Root>/radtodeg1'
                                        */
  real_T TransferFcn_A_j[2];           /* Computed Parameter: TransferFcn_A_j
                                        * Referenced by: '<Root>/Transfer Fcn'
                                        */
  real_T TransferFcn_C_o[2];           /* Computed Parameter: TransferFcn_C_o
                                        * Referenced by: '<Root>/Transfer Fcn'
                                        */
  real_T TransferFcn_A_m[2];           /* Computed Parameter: TransferFcn_A_m
                                        * Referenced by: '<S6>/Transfer Fcn'
                                        */
  real_T TransferFcn_C_f[2];           /* Computed Parameter: TransferFcn_C_f
                                        * Referenced by: '<S6>/Transfer Fcn'
                                        */
  real_T radtodeg2_Gain[8];            /* Expression: K
                                        * Referenced by: '<Root>/radtodeg2'
                                        */
  real_T radtodeg5_Gain;               /* Expression: 180/pi
                                        * Referenced by: '<Root>/radtodeg5'
                                        */
  real_T radtodeg4_Gain;               /* Expression: 180/pi
                                        * Referenced by: '<Root>/radtodeg4'
                                        */
  real_T Constant1_Value;              /* Expression: a1
                                        * Referenced by: '<Root>/Constant1'
                                        */
  real_T Constant2_Value;              /* Expression: a2
                                        * Referenced by: '<Root>/Constant2'
                                        */
  real_T Constant3_Value;              /* Expression: b1
                                        * Referenced by: '<Root>/Constant3'
                                        */
  real_T Constant4_Value;              /* Expression: a3
                                        * Referenced by: '<Root>/Constant4'
                                        */
  real_T Constant5_Value;              /* Expression: b2
                                        * Referenced by: '<Root>/Constant5'
                                        */
  real_T Gain_Gain_j;                  /* Expression: 0.5
                                        * Referenced by: '<Root>/Gain'
                                        */
  real_T Limit1_UpperSat;              /* Expression: 4
                                        * Referenced by: '<Root>/Limit1'
                                        */
  real_T Limit1_LowerSat;              /* Expression: 0.05
                                        * Referenced by: '<Root>/Limit1'
                                        */
  real_T Gain1_Gain_h;                 /* Expression: 0.5
                                        * Referenced by: '<Root>/Gain1'
                                        */
  real_T Limit2_UpperSat;              /* Expression: 4
                                        * Referenced by: '<Root>/Limit2'
                                        */
  real_T Limit2_LowerSat;              /* Expression: 0.05
                                        * Referenced by: '<Root>/Limit2'
                                        */
  real_T radtodeg6_Gain;               /* Expression: 180/pi
                                        * Referenced by: '<Root>/radtodeg6'
                                        */
  real_T Gain2_Gain_l;                 /* Expression: 32767/10
                                        * Referenced by: '<S1>/Gain2'
                                        */
  real_T Limit1_UpperSat_j;            /* Expression: 15000
                                        * Referenced by: '<S1>/Limit1'
                                        */
  real_T Limit1_LowerSat_m;            /* Expression: -15000
                                        * Referenced by: '<S1>/Limit1'
                                        */
  real_T Gain1_Gain_d;                 /* Expression: 32767/10
                                        * Referenced by: '<S1>/Gain1'
                                        */
  real_T Limit2_UpperSat_j;            /* Expression: 15000
                                        * Referenced by: '<S1>/Limit2'
                                        */
  real_T Limit2_LowerSat_m;            /* Expression: -15000
                                        * Referenced by: '<S1>/Limit2'
                                        */
  real_T TransferFcn2_A;               /* Computed Parameter: TransferFcn2_A
                                        * Referenced by: '<S10>/Transfer Fcn2'
                                        */
  real_T TransferFcn2_C;               /* Computed Parameter: TransferFcn2_C
                                        * Referenced by: '<S10>/Transfer Fcn2'
                                        */
  real_T SetCurrentAxisCommand1_P1;    /* Expression: axis
                                        * Referenced by: '<S10>/Set Current Axis' Command1'
                                        */
  real_T TransferFcn1_A;               /* Computed Parameter: TransferFcn1_A
                                        * Referenced by: '<S10>/Transfer Fcn1'
                                        */
  real_T TransferFcn1_C;               /* Computed Parameter: TransferFcn1_C
                                        * Referenced by: '<S10>/Transfer Fcn1'
                                        */
  real_T SetCurrentAxisCommand2_P1;    /* Expression: axis
                                        * Referenced by: '<S10>/Set Current Axis' Command2'
                                        */
  real_T Constant7_Value_l;            /* Expression: -0.545
                                        * Referenced by: '<S5>/Constant7'
                                        */
  real_T Constant1_Value_l;            /* Expression: 0
                                        * Referenced by: '<S7>/Constant1'
                                        */
  real_T Constant1_Value_d;            /* Expression: 0
                                        * Referenced by: '<S9>/Constant1'
                                        */
  uint8_T ManualSwitch_CurrentSetting; /* Computed Parameter: ManualSwitch_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch'
                                        */
  uint8_T ManualSwitch5_CurrentSetting;/* Computed Parameter: ManualSwitch5_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch5'
                                        */
  uint8_T ManualSwitch1_CurrentSetting;/* Computed Parameter: ManualSwitch1_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch1'
                                        */
  uint8_T ManualSwitch4_CurrentSetting;/* Computed Parameter: ManualSwitch4_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch4'
                                        */
  uint8_T ManualSwitch_CurrentSetting_e;/* Computed Parameter: ManualSwitch_CurrentSetting_e
                                         * Referenced by: '<S9>/Manual Switch'
                                         */
  uint8_T ManualSwitch5_CurrentSetting_d;/* Computed Parameter: ManualSwitch5_CurrentSetting_d
                                          * Referenced by: '<S9>/Manual Switch5'
                                          */
  uint8_T ManualSwitch1_CurrentSetting_a;/* Computed Parameter: ManualSwitch1_CurrentSetting_a
                                          * Referenced by: '<S9>/Manual Switch1'
                                          */
  uint8_T ManualSwitch4_CurrentSetting_k;/* Computed Parameter: ManualSwitch4_CurrentSetting_k
                                          * Referenced by: '<S9>/Manual Switch4'
                                          */
  uint8_T ManualSwitch2_CurrentSetting;/* Computed Parameter: ManualSwitch2_CurrentSetting
                                        * Referenced by: '<S10>/Manual Switch2'
                                        */
  uint8_T ManualSwitch1_CurrentSetting_f;/* Computed Parameter: ManualSwitch1_CurrentSetting_f
                                          * Referenced by: '<S10>/Manual Switch1'
                                          */
  uint8_T ManualSwitch2_CurrentSetting_l;/* Computed Parameter: ManualSwitch2_CurrentSetting_l
                                          * Referenced by: '<S7>/Manual Switch2'
                                          */
  uint8_T ManualSwitch3_CurrentSetting;/* Computed Parameter: ManualSwitch3_CurrentSetting
                                        * Referenced by: '<S7>/Manual Switch3'
                                        */
  uint8_T ManualSwitch2_CurrentSetting_m;/* Computed Parameter: ManualSwitch2_CurrentSetting_m
                                          * Referenced by: '<S9>/Manual Switch2'
                                          */
  uint8_T ManualSwitch3_CurrentSetting_d;/* Computed Parameter: ManualSwitch3_CurrentSetting_d
                                          * Referenced by: '<S9>/Manual Switch3'
                                          */
};

/* Real-time Model Data Structure */
struct tag_RTM_Control_linear {
  const char_T *path;
  const char_T *modelName;
  struct SimStruct_tag * *childSfunctions;
  const char_T *errorStatus;
  SS_SimMode simMode;
  RTWExtModeInfo *extModeInfo;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * ModelData:
   * The following substructure contains information regarding
   * the data used in the model.
   */
  struct {
    void *blockIO;
    const void *constBlockIO;
    void *defaultParam;
    ZCSigState *prevZCSigState;
    real_T *contStates;
    real_T *derivs;
    void *zcSignalValues;
    void *inputs;
    void *outputs;
    boolean_T *contStateDisabled;
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T blkStateChange;
    real_T odeY[10];
    real_T odeF[3][10];
    ODE3_IntgData intgData;
  } ModelData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T checksums[4];
    uint32_T options;
    int_T numContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * SpecialInfo:
   * The following substructure contains special information
   * related to other components that are dependent on RTW.
   */
  struct {
    const void *mappingInfo;
    void *xpcData;
  } SpecialInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    time_T stepSize1;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    void *timingData;
    real_T *varNextHitTimesList;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;

  /*
   * Work:
   * The following substructure contains information regarding
   * the work vectors in the model.
   */
  struct {
    void *dwork;
  } Work;
};

/* Block parameters (auto storage) */
extern Parameters_Control_linear Control_linear_P;

/* Block signals (auto storage) */
extern BlockIO_Control_linear Control_linear_B;

/* Continuous states (auto storage) */
extern ContinuousStates_Control_linear Control_linear_X;

/* Block states (auto storage) */
extern D_Work_Control_linear Control_linear_DWork;

/* Model entry point functions */
extern void Control_linear_initialize(void);
extern void Control_linear_output(void);
extern void Control_linear_update(void);
extern void Control_linear_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Control_linear *const Control_linear_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Control_linear'
 * '<S1>'   : 'Control_linear/Heli-model with sat and disturbance 1'
 * '<S2>'   : 'Control_linear/MATLAB Function1'
 * '<S3>'   : 'Control_linear/MATLAB Function4'
 * '<S4>'   : 'Control_linear/Subsystem'
 * '<S5>'   : 'Control_linear/Subsystem1'
 * '<S6>'   : 'Control_linear/Subsystem2'
 * '<S7>'   : 'Control_linear/goal_ele'
 * '<S8>'   : 'Control_linear/goal_pitch'
 * '<S9>'   : 'Control_linear/goal_travel'
 * '<S10>'  : 'Control_linear/Heli-model with sat and disturbance 1/Subsystem'
 * '<S11>'  : 'Control_linear/goal_ele/Signal From Workspace'
 * '<S12>'  : 'Control_linear/goal_ele/Signal From Workspace1'
 * '<S13>'  : 'Control_linear/goal_ele/Signal From Workspace2'
 * '<S14>'  : 'Control_linear/goal_pitch/Signal From Workspace'
 * '<S15>'  : 'Control_linear/goal_pitch/Signal From Workspace1'
 * '<S16>'  : 'Control_linear/goal_pitch/Signal From Workspace2'
 * '<S17>'  : 'Control_linear/goal_travel/Signal From Workspace'
 * '<S18>'  : 'Control_linear/goal_travel/Signal From Workspace1'
 * '<S19>'  : 'Control_linear/goal_travel/Signal From Workspace2'
 */
#endif                                 /* RTW_HEADER_Control_linear_h_ */
